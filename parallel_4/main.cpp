#include <iostream> 
#include <string> 
#include "init.hpp"
#include "parameters.hpp"
#include "sieving.hpp"
#include "sieving_worker.hpp"
#include "sieving_coord.hpp"
#include "candidates_main.hpp"
#include "candidates_helper.hpp"
#include <mpi.h>

/*** BE CAREFUL WITH VALUES 2 AND -1 ***/

using namespace quadratic_sieve;
namespace quadratic_sieve{
	int mpi_rank;
	int mpi_size;
} 
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
		generate_factor_base(); // generate base of all factors in the sieve
		sieving_coord();
		MPI_Barrier(MPI_COMM_WORLD); /*** Not really required, included for synchronization of processes before going into next phase ***/

		find_choice_vector();
	}
	/*** Worker processes ***/
	else {
		recieve_number();
		setup_parameters();
		generate_factor_base();
		sieving_work();				// answer sieving requests
		MPI_Barrier(MPI_COMM_WORLD);   /*** Not really required, included for synchronization of processes before going into next phase ***/
		answer_multiplication_requests(); // answer multiplications
	}
	MPI_Finalize();	
	return 0;
}
