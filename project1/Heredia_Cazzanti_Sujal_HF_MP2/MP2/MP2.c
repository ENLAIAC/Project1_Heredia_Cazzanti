#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <trexio.h>  //include the 'exit' function to terminate the program if any read goes wrong


///////////////////////////////////// AUXILIARY FUNCTIONS //////////////////////////
//Two-electron integrals obey 8-fold permutational symmetry
//Since TREXIO stores only one permutation for each quartet, we need a way to retrieve <pq|rs>
//even if the requested permutation is not the one stored.
//Here we canonicalize the 4 indexes (p,q,r,s) by taking the lexicographically smallest
//tuple among the 8 equivalent permutations, and pack it into a 64-bit integer key.

static inline uint64_t pack4_u16(uint16_t a, uint16_t b, uint16_t c, uint16_t d){
	return ((uint64_t)a << 48) | ((uint64_t)b << 32) | ((uint64_t)c << 16) | (uint64_t)d;
}

static inline int tuple_lt(uint16_t a1,uint16_t b1,uint16_t c1,uint16_t d1,
                           uint16_t a2,uint16_t b2,uint16_t c2,uint16_t d2){
	if (a1 != a2) return a1 < a2;
	if (b1 != b2) return b1 < b2;
	if (c1 != c2) return c1 < c2;
	return d1 < d2;
}

static uint64_t canonical_key_8fold(int p, int q, int r, int s){
	//<ij|kl> = <il|kj> = <kl|ij> = <kj|il> = <ji|lk> = <li|jk> = <lk|ji> = <jk|li>
	uint16_t t[8][4] = {
		{(uint16_t)p,(uint16_t)q,(uint16_t)r,(uint16_t)s}, // <pq|rs>
		{(uint16_t)p,(uint16_t)s,(uint16_t)r,(uint16_t)q}, // <ps|rq>
		{(uint16_t)r,(uint16_t)s,(uint16_t)p,(uint16_t)q}, // <rs|pq>
		{(uint16_t)r,(uint16_t)q,(uint16_t)p,(uint16_t)s}, // <rq|ps>
		{(uint16_t)q,(uint16_t)p,(uint16_t)s,(uint16_t)r}, // <qp|sr>
		{(uint16_t)s,(uint16_t)p,(uint16_t)q,(uint16_t)r}, // <sp|qr>
		{(uint16_t)s,(uint16_t)r,(uint16_t)q,(uint16_t)p}, // <sr|qp>
		{(uint16_t)q,(uint16_t)r,(uint16_t)s,(uint16_t)p}  // <qr|sp>
	};

	uint16_t b0=t[0][0], b1=t[0][1], b2=t[0][2], b3=t[0][3];
	for (int m=1; m<8; m++){
		if (tuple_lt(t[m][0],t[m][1],t[m][2],t[m][3], b0,b1,b2,b3)){
			b0=t[m][0]; b1=t[m][1]; b2=t[m][2]; b3=t[m][3];
		}
	}
	return pack4_u16(b0,b1,b2,b3);
}

//We store in an array and sort it. Then we can retrieve any <pq|rs> with bsearch.
typedef struct {
	uint64_t key;
	double val;
} eri_kv_t;

static int cmp_eri_kv(const void* a, const void* b){
	const eri_kv_t* x = (const eri_kv_t*)a;
	const eri_kv_t* y = (const eri_kv_t*)b;
	if (x->key < y->key) return -1;
	if (x->key > y->key) return  1;
	return 0;
}

static double eri_get(const eri_kv_t* arr, int64_t n, int p, int q, int r, int s){
	eri_kv_t needle;
	needle.key = canonical_key_8fold(p,q,r,s);

	const eri_kv_t* found = (const eri_kv_t*) bsearch(
		&needle, arr, (size_t)n, sizeof(eri_kv_t), cmp_eri_kv
	);

	if (found == NULL) return 0.0;
	return found->val;
}


int main(void){
	//////////////////////////////////// TREXIO VARIABLES INITIALIZATION ///////////////////////////////
	
	trexio_exit_code rc; //This variable stores a message about the status of the trexio.h function. If is succesfully called and ended it stores a 'TREXIO SUCCESS', otherwise it sotres the error arised
	trexio_t* trexio_file=trexio_open("h2o.h5", 'r', TREXIO_AUTO, &rc);

	///////////////////////////////////// VARIABLES DECLARATION PART /////////////////////////////////////
	
	int num_elec; //Variable devoted to store the number of electron: useful to track the number of occupied 
		      //orbitals, i.e., those that contribute to the HF energy
	int mo; //Stores the number of molecular orbital: used to define the dimension of the two electron integral 
		//matr	ix. Here, mo includes both virtual and occupied orbitals
	int64_t integrals; //amount of 2e non-zero integrals. Here, occupied-occupied, occupied-virtual and 
			   //virtual-virtual 2-electron integrals are included. The ones that contribute to the energy 
			   //are the former ones.
	double* mo_energy; //Array devoted to store the MO energies eps_p
	int* indexes; //Array storing the 4 indexes associated to each 2e integral
	double* two_el_int; //Array storing the values <pq|rs> corresponding to the indexes above
	eri_kv_t* eri_table; //Canonicalized (key,value) array for ERIs
	double emp2=0.0; //MP2 correlation energy

	///////////////////////////////////////////// PROGRAM STARTS ////////////////////////////////////////
	//Reading from file phase:
	//Number of occupied orbitals
	rc = trexio_read_electron_up_num(trexio_file, &num_elec);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the number of electrons: %s\n", trexio_string_of_error(rc));
		exit(1);
	}

	//Number of MOs
	rc = trexio_read_mo_num(trexio_file, &mo);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the number of molecular orbitals: %s\n", trexio_string_of_error(rc));
		exit(1);
	}

	//MO energies 
	mo_energy = malloc(mo*sizeof(double));
	if ( mo_energy == NULL ){
		printf("Memory allocation went wrong");
		exit(1);
	}
	rc = trexio_read_mo_energy(trexio_file, mo_energy);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the MO energies: %s\n", trexio_string_of_error(rc));
		exit(1);
	}

	//Number of non-zero 2-electron integrals stored in sparse format
	rc = trexio_read_mo_2e_int_eri_size(trexio_file, &integrals);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the number of 2-electron integrals: %s\n", trexio_string_of_error(rc));
		exit(1);
	}

	//Allocating (reserving) a finite amount of memory to store the indexes of the orbitals in the two-electron
	//integrals. Each integral is associated with 4 integers (i.e., the indexes).
	indexes = malloc(integrals*4*sizeof(int)); 
	//Memory allocation success verification
	if ( indexes == NULL ){
		printf("Memory allocation went wrong");
		exit(1);
	}

	//Same as 'indexes' but here the variable is devoted to store the integral values
	two_el_int = malloc(integrals*sizeof(double));
	if ( two_el_int == NULL ){
		printf("Memory allocation went wrong");
		exit(1);
	}

	//Finally reading the two-electrons integrals. Bear in mind 'indexes' and 'two_el_int' are pointers.
	rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &integrals, indexes, two_el_int);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the 2-electron integrals: %s\n", trexio_string_of_error(rc));
		exit(1);
	}

	//////////////////////////////////////// SYMMETRY HANDLING /////////////////////////////
	//We transform the sparse (indexes,value) storage into a sorted (key,value) table where the key is canonical
	//with respect to the 8-fold symmetry. Then <pq|rs> can be retrieved with eri_get(...).

	eri_table = malloc(integrals*sizeof(eri_kv_t));
	if ( eri_table == NULL ){
		printf("Memory allocation went wrong");
		exit(1);
	}

	for (int64_t n=0; n<integrals; n++){
		int p = indexes[4*n + 0];
		int q = indexes[4*n + 1];
		int r = indexes[4*n + 2];
		int s = indexes[4*n + 3];

		eri_table[n].key = canonical_key_8fold(p,q,r,s);
		eri_table[n].val = two_el_int[n];
	}

	qsort(eri_table, (size_t)integrals, sizeof(eri_kv_t), cmp_eri_kv);

	//////////////////////////////////////// MP2 ENERGY CALCULATION //////////////////////////////////

	for (int i=0; i<num_elec; i++){
		for (int j=0; j<num_elec; j++){
			for (int a=num_elec; a<mo; a++){
				for (int b=num_elec; b<mo; b++){

					double ijab = eri_get(eri_table, integrals, i, j, a, b); // <ij|ab>
					if (ijab == 0.0) continue; 

					double ijba = eri_get(eri_table, integrals, i, j, b, a); // <ij|ba>
					double denom = mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b];

					emp2 += ijab * ( (2.0*ijab) - ijba ) / denom;
				}
			}
		}
	}

	printf("MP2 correlation energy: %f \n", emp2);

	/////////////////////////////////////// MEMORY DEALLOCATION PHASE ///////////////////////////////////
	trexio_close(trexio_file);

	free(mo_energy);
	mo_energy=NULL;

	free(indexes);
	indexes=NULL;

	free(two_el_int);
	two_el_int=NULL;

	free(eri_table);
	eri_table=NULL;
}
