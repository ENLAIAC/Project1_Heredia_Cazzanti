#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>

int main(void){
//////////////////////////////////////////////////////////TREXIO VARIABLES INITIALIZATION///////////////////////////////
	
	trexio_exit_code rc; //This variable stores a message about the status of the trexio.h function. If is succesfully called and ended it stores a 'TREXIO SUCCESS', otherwise it sotres the error arised
	trexio_t* trexio_file=trexio_open("h2o.h5", 'r', TREXIO_AUTO, &rc);

/////////////////////////////////////////VARIABLES DECLARATION PART/////////////////////////////////////////////////////
	
	double Vnn; //This variable stores the nuclear repulsion
	int num_elec; //Variable devoted to store the number of electron
	int mo; //Stores the number of molecular orbital: used to define the dimension of the two electron integral matrix 
	int64_t integrals; //amount of 4e non-zero integrals
	double energy;

	/////////////////////////////////////////////PROGRAM STARTS//////////////////////////////////////////////////
	rc=trexio_read_nucleus_repulsion(trexio_file, &Vnn);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the nucleus repulsion: %s\n", trexio_string_of_error(rc));
	}
	rc=trexio_read_electron_up_num(trexio_file, &num_elec);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the number of electrons: %s\n", trexio_string_of_error(rc));
	}
	rc=trexio_read_mo_num(trexio_file, &mo);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the number of orbitals: %s\n", trexio_string_of_error(rc));
	}
	double Ven[mo][mo]; //malloc(mo*mo*sizeof(double)); //Pointer that stores the electron-nucleus repulsion potential energy
	rc=trexio_read_mo_1e_int_core_hamiltonian(trexio_file, &Ven[0][0]);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the 1-electron orbitals: %s\n", trexio_string_of_error(rc));
	}
	rc=trexio_read_mo_2e_int_eri_size(trexio_file, &integrals);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the number of non vanishing 2-electron orbitals: %s\n", trexio_string_of_error(rc));
	}
	int* indexes = malloc(integrals*4*sizeof(int)); //Is a pointer: stores the memory label that was reserved to store the indexes of the two-electron integrals
	//memory allocation success verification
	if ( indexes == NULL ){
		printf("Memory allocation went wrong");
		exit(1);
	}
	double* two_el_int = malloc(integrals*sizeof(double)); //Is a pointer: it stores the memory reference that will host the 2e integrals values
	double E_coul=0; //This variable is designed to cumulate the Coulomb repulsion contributions, hence integrals whose indexes are i=j/k=l (ii|jj)
	double E_xc=0; //This variable yield the exchange energy. Its contributions are the integrals whose indexes satisfy i=l/j=k (ij|ji)
	//memory allocation success verification
	if ( two_el_int == NULL ){
		printf("Memory allocation went wrong");
		exit(1);
	}
	rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &integrals, indexes, two_el_int);
	int i, j, k, l; //2-electrons integral indexes
	//i=indexes[5*4+2];
	energy=Vnn; //is summed one time onlyi
	printf("Nuclear repulsion energy:%f", energy);

	//Follows the calculation of one electron integral sum: The integrals store by the read of one electrons integrals return the integrals of all orbitals both occupied and virtual, thus a first filter should be set. Secondly, among the occupied orbitals only the diagonal terms (i.e., hii)n contributes to the 1-el energy. Storing the amount of mo was required to provide the right size to the array that stores the one-electron integrals, but the sum (i.e., the loop) doesn't run over the whole set of MOs
	
	int one_el_sum=0;
	for (int n=0; n<=num_elec; n++){	
		one_el_sum+=2*Ven[n][n];
	}

	energy+=one_el_sum;
	printf("Nuclear repulsion energy:%f", energy);
	double two_el_sum=0; //Stores the sum of the two electron integrals.
	for (int n=0; n<=num_elec; n++){
		i = indexes[n*4+0];
		j = indexes[n*4+1];
		k = indexes[n*4+2];
		l = indexes[n*4+3];
		if (i==k && j==l){
			printf("Indici:%d %d %d %d \n", i, j, k, l);
			two_el_sum+=(2*two_el_int[i][j][k][l]-two_el_int[i][k][j][l]);
		}
	}
	energy+=two_el_sum;
	printf("Total energy of water:%f ",energy);
	trexio_close(trexio_file);
	//free(Ven);
	//Ven=NULL;
	free(indexes);
	indexes=NULL;
	free(two_el_int);
	two_el_int=NULL;
}
