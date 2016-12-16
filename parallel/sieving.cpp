#include <mpi.h>
#include "sieving.hpp"
#include "parameters.hpp"
#include <mpi.h>
#include <boost/multiprecision/cpp_int.hpp>
#include "mod_arithm.hpp"
#include "sparse_binary_vector.hpp"

namespace quadratic_sieve{

	std::vector<cpp_int> smooth_numbers;
	std::vector<vec1> exponents;	/*** matrix containing vectors found so far ***/
	std::vector<vec> exponents_binary;
	int found_vectors_num;

	int total_number=-1;	/*** Total_number of sieved vectors the matrix should have in the end.
													* When all processes finish sieving this information is used to ensure 
													* they all recieved all vectors
													***/ 
	MPI_Request	*broad_req;
	b_info * info_array;
	b_data * data_array;
	int *dumm;
	int *phases;   /*** 0 - waiting for info from node i , 1 - waiting for data ***/
	/*** processes vectors recieved from sieving ***/
	void init_matrix_communication(){
		broad_req=new MPI_Request[mpi_size];
		phases= new int[mpi_size];	
		info_array=new b_info[mpi_size];
		data_array=new b_data[mpi_size];
		dumm=new int[mpi_size];
		for ( int i=0;i<mpi_size;i++)
			phases[i]=0;
		for ( int i=0;i<mpi_size;i++){
			if ( i==mpi_rank || i == MY_ROOT ) continue;
			else{
			//MPI_Ibcast(info_array[i].info,4,MPI_INT,i,MPI_COMM_WORLD,&broad_req[i]);
			//if ( i== 2)
			MPI_Irecv(dumm+i,1,MPI_INT,i,MPI_COMM_WORLD,&broad_req[i]);
			}
		}
	}

	void step_matrix_construction_protocol(){
		for ( int i=0;i<mpi_size;i++){
			if (  i==MY_ROOT ||  i == mpi_rank ) continue;
			int flag;
			MPI_Status status;
			MPI_Test(&broad_req[i],&flag,&status);
			if ( flag == true ){
				if ( phases[i]==0 ){
					printf("dest: %d %d\n",i,mpi_rank);
					int sz=0; for ( int j=0;j<3;j++) sz+=info_array[i].info[j];
					printf("SZ(dest): %d\n",sz);std::cout.flush();
					//data_array[i].data=new int[sz];
					//MPI_Ibcast(data_array[i].data,sz,MPI_INT,i,MPI_COMM_WORLD,broad_req+i);
					phases[i]=1;
				}	
				else if ( phases[i]==1){
					//process_data(	info_array[i],data_array[i]);
					//phases[i]=0;
				}
			}
		}
	}

	/*** major function ***/
	void process_data(b_info &info, b_data &data){
		printf("processing data");
	}
	/*** returns true if we formed completely matrix ***/
	bool all_sieving_recieved(){
		return false;
	}
	void finalize_matrix_communication(){
		delete [] info_array;
		delete []data_array;
		for ( int i=0;i<mpi_size;i++){
			MPI_Request_free(broad_req+i);
		}
		delete []broad_req;
	}
}
