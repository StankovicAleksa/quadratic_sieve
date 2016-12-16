#include "candidates_helper.hpp"
#include <mpi.h>
#include "parameters.hpp"
#include "mod_arithm.hpp"
namespace quadratic_sieve{
	void answer_multiplication_requests(){
		int *dummy_buff=0;
		int *recv_data;
		int chunk_size; // if sz is -1, stop the operation
		
		while (true ){
			MPI_Bcast(&chunk_size,1,MPI_INT,MY_ROOT,MPI_COMM_WORLD);
			if ( chunk_size == -1 ) {
				break;	 // we have finished! :)
			}
			recv_data=new int[chunk_size];
			MPI_Scatter(dummy_buff,chunk_size,MPI_INT,recv_data,chunk_size,MPI_INT,MY_ROOT,MPI_COMM_WORLD);
			// Essentialy same code as in candidates_main
			vec x_temp;x_temp.resize(base_size+1);
			for ( int i=0;i<chunk_size;i++){
				if ( recv_data[i]==-1 ) break;
				else x_temp.push_back(recv_data[i]);
			}
			vec y_temp; y_temp.resize(base_size+1);
			matrix_multiply(y_temp,exponents_binary,x_temp);
			int final_data[base_size+1];
			for ( int i=0;i<base_size+1;i++){
				final_data[i]=0;
			}
			std::vector<int> vec_data=y_temp.getData();
			for ( std::vector<int>::iterator it=vec_data.begin();it!=vec_data.end();it++){
				final_data[*it]=1;
			}

			MPI_Reduce(final_data,dummy_buff,base_size+1,MPI_INT,MPI_BXOR,MY_ROOT,MPI_COMM_WORLD);
			delete recv_data;
		}
	}	
}
