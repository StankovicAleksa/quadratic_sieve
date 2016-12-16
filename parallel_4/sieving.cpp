#include <mpi.h>
#include "sieving.hpp"
#include "parameters.hpp"
#include <mpi.h>
#include <boost/multiprecision/cpp_int.hpp>
#include "mod_arithm.hpp"
#include "sparse_binary_vector.hpp"
#include <set>

struct comp_vector{
	int index;
	b_info info;
	b_data data;
};
bool operator < (comp_vector v1,comp_vector v2){
	return v1.index<v2.index;
}
namespace quadratic_sieve{


	/*** Global variables ***/
	std::vector<cpp_int> smooth_numbers;
	std::vector<vec1> exponents;	/*** matrix containing vectors found so far ***/
	std::vector<vec> exponents_binary;


	/*** variables for synchronization of matrix construction ***/
	int processed_columns=0;
	std::multiset<comp_vector> vec_list;
	MPI_Request final_bcast;
	bool final_processed=false;
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
			if ( i == MY_ROOT ) continue;
			else{
				MPI_Irecv(info_array[i].info,4,MPI_INT,i,4,MPI_COMM_WORLD,&broad_req[i]);
			}
		}
	 if ( mpi_rank != MY_ROOT)	
		MPI_Ibcast(&total_number,1,MPI_INT,MY_ROOT,MPI_COMM_WORLD,&final_bcast);  // info from coordinator that all processes have complete matrix
	} 

	void step_matrix_construction_protocol(){
		for ( int i=0;i<mpi_size;i++){
			if (  i==MY_ROOT ) continue;
			int flag;
			MPI_Status status;
			MPI_Test(&broad_req[i],&flag,&status);
			if ( flag == true ){
				if ( phases[i]==0 ){
					int sz=0; for ( int j=0;j<3;j++) sz+=info_array[i].info[j];
					data_array[i].data=new int[sz];
					data_array[i].sz=sz;
					MPI_Irecv(data_array[i].data,sz,MPI_INT,i,5,MPI_COMM_WORLD,broad_req+i);
					phases[i]=1;
				}	
				else if ( phases[i]==1){
					process_data(	info_array[i],data_array[i]);
					phases[i]=0;
					MPI_Irecv(info_array[i].info,4,MPI_INT,i,4,MPI_COMM_WORLD,&broad_req[i]);
				}
			}
		}
		if(!final_processed && mpi_rank != MY_ROOT){
			int flag; MPI_Status status;
			MPI_Test(&final_bcast,&flag,&status);
			if ( flag == true){
				final_processed=true;
			}
		}
	}
	void set_final_matrix_width(int width){
		total_number=width;	
	}
	/*** major function ***/
	void process_data(b_info &info, b_data &data){
		comp_vector v;
		v.info=info;v.data=data;
		v.index=info.info[3];
		vec_list.insert(v);
		while ( !vec_list.empty() ){
			v=*vec_list.begin();
			if ( v.index == processed_columns ){
				vec_list.erase(vec_list.begin());
				/*** unpack data ***/
				int bin_sz=v.info.info[0];
				int full_sz=v.info.info[1];
				int smooth_sz=v.info.info[2];
				vec exp_bin; exp_bin.resize(base_size+1);
				vec1 exp; exp.resize(base_size+1); exp.clear();
				int index=0;
				while (index<bin_sz){
					if ( v.data.data[index]==-1 ){
						exponents_binary.push_back(exp_bin);
						exp_bin=vec(); exp_bin.resize(base_size+1);
					} else{
						exp_bin.push_back(v.data.data[index]);
					}
					index++;
				}
				while( index<bin_sz+full_sz ){
					if ( v.data.data[index]==-1 ){
						exponents.push_back(exp);
						exp=vec1(); exp.resize(base_size+1);
						index++;
					} else {
						exp(v.data.data[index])=v.data.data[index+1];
						index+=2;
					}
				}
				// finally add smooth numbers and update processed_columns
				for ( int i=0;i<smooth_sz;i++){
					smooth_numbers.push_back(sqrt_n+v.data.data[index]);
					index++;
					processed_columns++;
				}	
			}
			else {
				break;
			}
		}
	}
	/*** returns true if we formed completely matrix ***/
	bool all_sieving_recieved(){
		return total_number==processed_columns;
	}
	void finalize_matrix_communication(){
		delete [] info_array;
		delete []data_array;
		for ( int i=0;i<mpi_size;i++){
			//MPI_Cancel(broad_req+i);
			//MPI_Request_free(broad_req+i);
		}
		delete []broad_req;
	}
}
