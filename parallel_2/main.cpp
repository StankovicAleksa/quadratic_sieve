#include <iostream> 
#include <string> 
#include "init.hpp"
#include "parameters.hpp"
#include "sieving.hpp"
#include "sieving_worker.hpp"
#include "sieving_coord.hpp"
#include "candidates.hpp"
#include <mpi.h>

/*** BE CAREFUL WITH VALUES 2 AND -1 ***/

using namespace quadratic_sieve;
namespace quadratic_sieve{
	int mpi_rank;
	int mpi_size;
} 
boost::multiprecision::cpp_int n; // the number to be factorized 
int main(int argc,char *argv[]) {
	
	/*** MPI stuff ***/	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

	/*** coordinator process ***/
	if ( mpi_rank == MY_ROOT){
		input_number(); // input the number
		broadcast_number(); // broadcast the number to the all other processes
		setup_parameters();  // setup the environment parameters, we need this only for coordinator
		generate_factor_base();
		printf("Sieving\n");
		sieving_coord();
		printf("Searching\n");
		//find_choice_vector();
	}
	/*** Worker processes ***/
	else {
		recieve_number();
		setup_parameters();
		generate_factor_base();
		sieving_work();
		// answer multiplications
	}
	MPI_Finalize();	
	return 0;
}
