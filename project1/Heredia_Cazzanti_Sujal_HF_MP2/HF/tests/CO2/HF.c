#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>  //include the 'exit' function to terminate the program if any read goes wrong

int main(void){
	//////////////////////////////////// TREXIO VARIABLES INITIALIZATION ///////////////////////////////
	
	trexio_exit_code rc; //This variable stores a message about the status of the trexio.h function. If is succesfully called and ended it stores a 'TREXIO SUCCESS', otherwise it sotres the error arised
	trexio_t* trexio_file=trexio_open("co2.h5", 'r', TREXIO_AUTO, &rc);

	///////////////////////////////////// VARIABLES DECLARATION PART /////////////////////////////////////
	
	double Vnn; //This variable stores the nuclear repulsion
	int num_elec; //Variable devoted to store the number of electron: useful to track the number of occupied 
		      //orbitals, i.e., those that contribute to the HF energy
	int mo; //Stores the number of molecular orbital: used to define the dimension of the two electron integral 
		//matr	ix. Here, mo includes both virtual and occupied orbitals
	int64_t integrals; //amount of 2e non-zero integrals. Here, occupied-occupied, occupied-virtual and 
			   //virtual-virtual 2-electron integrals are included. The ones that contribute to the energy 
			   //are the former ones.
	double energy; //This variable will store the final energy

	///////////////////////////////////////////// PROGRAM STARTS ////////////////////////////////////////
	//Reading from file phase:
	//- Nuclear-Nuclear repulsion (Vnn)
	rc=trexio_read_nucleus_repulsion(trexio_file, &Vnn); //Trexio requires to pass a pointer, i.e., '&Vnn' 
	if ( rc != TREXIO_SUCCESS){ 
		printf ("Error reading the nucleus repulsion: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	printf("Nuclear-Nuclear repulsion energy: %f \n", Vnn);
	//- Number of spin-up electrons: stands for the number of occupied spatial orbitals
	rc=trexio_read_electron_up_num(trexio_file, &num_elec);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the number of electrons: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	//- Reading the number of molecular orbital (both virtual and occupied). To note that only occupied contribute
	//to 1electron and two electron energy
	rc=trexio_read_mo_num(trexio_file, &mo); 
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the number of orbitals: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	//- Reading the electorn-nuclear repulsion (i.e., one electron orbitals). The terms that contribute to the 
	//energy are the ones that couples the same occupied orbital= <i|h|i> where i \in {occupied}
	double Ven[mo][mo]; //According to the amount of MOs in the molecule (still both virtual and occupeid). 
			    //The dimension of the one-electron operator is mo x mo (see the previous reading)
			    //One electron operator computes both kinetic energy of i-th electron and electron-nucleus
			    //potential energy
	rc=trexio_read_mo_1e_int_core_hamiltonian(trexio_file, &Ven[0][0]);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the 1-electron orbitals: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	//- Reading the amount of non-vanishing integrals within the set of two-electron integrals
	rc=trexio_read_mo_2e_int_eri_size(trexio_file, &integrals);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the number of 2-electron orbitals: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	//- Allocating (reserving) a finite amount of memory to store the indexes of the orbitals in the two-electron
	//integrals. Trexio stores in a sparse matrix 2-el integrals regradless the orbitals involved are occupied or
	//virtual and whether that term actually contribute to the two electron energy. Each integral is associated 
	//with 4 integers (i.e., the indexes).
	int* indexes = malloc(integrals*4*sizeof(int)); 
	//Memory allocation success verification
	if ( indexes == NULL ){
		printf("Memory allocation went wrong");
		exit(1);
	}
	//- Same as 'indexes' but here the variable is devoted to store the two-electron integrals (i.e., <ij|kl>)
	//While 'integrals' stores the amount of integrals, 'two_el_int' stores the evaluation of these integrals. 
	//Again, we're interested in those 2-el integrals of the form <ij|ij> for Coulomb energy contributions, 
	//and to those <ij|ji> for exchange energy.  
	double* two_el_int = malloc(integrals*sizeof(double));
	if ( two_el_int == NULL ){
		printf("Memory allocation went wrong");
		exit(1);
	}
	//- Finally reading the two-electrons integrals. Bear in mind 'indexes' and 'two_el_int' are pointers
	
	rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &integrals, indexes, two_el_int);
	if ( rc != TREXIO_SUCCESS){
		printf ("Error reading the 2-electron orbitals: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	
	//////////////////////////////////////// ENERGY CALCULATION //////////////////////////////////
	
	
	//Contribution of the nuclear-nuclear repulsion (only once)
	energy=Vnn; 	
	printf("Nuclear repulsion energy: %f \n", energy); //To verify the energy adds up to the nuclear-nuclear repulsion

	//1-el energy computation
	double one_el_en=0; //One electron energy. REMINDER: Only <i|h|i> terms contribute to it
	for (int n=0; n<num_elec; n++){	
		one_el_en=one_el_en+(2*Ven[n][n]);
	}
	printf("One electron energy: %f \n", one_el_en);
	energy+=one_el_en;
	printf("Nuclear-nuclear repulsion + one electron energy: %f \n", energy);

	//2-el energy computation
	double E_coul=0; //This variable is designed to cumulate the Coulomb repulsion contributions, hence integrals 
			 //whose indexes are i=k/j=l <ij|ij>
	double E_xc=0; //This variable yield the exchange energy. Its contributions are the integrals whose indexes 
		       //satisfy i=l/j=k <ij|ji>
	int i, j, k, l; //2-electrons integral indexes
	double two_el_en=0; //Stores the sum of the two electron integrals.
	int nJ=0, nK=0;
	
	
	
	for (int n=0; n<integrals; n++){
		i = indexes[n*4+0];
		j = indexes[n*4+1];
		k = indexes[n*4+2];
		l = indexes[n*4+3];
		if (i<num_elec && j<num_elec && k<num_elec && l<num_elec){
			printf("Indexes:%d %d %d %d \n", i, j, k, l);
			if (i==k && j==l){
				nJ++;
				printf("COULOMB YES");
				if (i==j){
					nK++;
					//printf("Indexes:%d %d %d %d \n", i, j, k, l);
					E_coul+=two_el_int[n];
					printf("Coulomb energy: %f \n", E_coul);
				}
				else{
					//printf("Indexes:%d %d %d %d \n", i, j, k, l);
					E_coul+=4*two_el_int[n];
					printf("Coulomb energy: %f \n", E_coul);
				}
			}
			else if ((i==l && j==k) || (i==j && k==l)) {
				nK++;
				printf("EXCHANGE YES");
				//printf("Indexes:%d %d %d %d \n", i, j, k, l);
				E_xc-=2*two_el_int[n];
				printf("Exchange energy: %f \n", E_xc);
			}
		}
	}
	printf("Amount of Coulomb contributions: %d \n", nJ);
	printf("Amount of Exchange contributions: %d \n", nK);
	two_el_en=E_coul+E_xc; //Here the contributions are summed because 'E_xc' already carries the minus sign
	printf("Two electron_energy: %f \n", two_el_en);			     
	energy+=two_el_en;
	printf("Final energy: %f \n", energy);
	
	
	/////////////////////////////////////// MEMORY DEALLOCATION PHASE ///////////////////////////////////
	trexio_close(trexio_file);
	free(indexes);
	indexes=NULL;
	free(two_el_int);
	two_el_int=NULL;
}
