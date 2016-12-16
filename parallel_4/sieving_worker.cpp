#include "sieving_worker.hpp"
#include "sieving.hpp"
#include "parameters.hpp"
#include <mpi.h>
#include "mod_arithm.hpp"
#include <queue>

using namespace quadratic_sieve;
namespace quadratic_sieve{

	int phase=0;  // the phase process is in: 0 - free, 1 - sieving
	int info[2]; // contains start and end for the sieving process
	int matrix_index;     // index in matrix given for our computation
	int sieved_vectors=0;

	MPI_Request req; // general request 



	/*** SIEVING INFO ***/
	/*** place where we store temporary results after sieving ***/
	std::vector<int> found_smooth_numbers;
	std::vector<vec1> found_exponents;	/*** matrix containing vectors found so far ***/
	std::vector<vec> found_exponents_binary;

	/*** BROADCAST INFO ***/
	/*** Queue for broadcasting ***/
	std::queue<b_data> q_b_data;
	std::queue<b_info> q_i_data;


	void sieving_work(){
		init_matrix_communication();  // we initialize communication variables required for construction of matrix
		// initialization
		MPI_Request req; // general request 
		MPI_Irecv(info,2,MPI_INT,MY_ROOT,0,MPI_COMM_WORLD,&req);

		int flag; MPI_Status status;
		while (!all_sieving_recieved()){
			MPI_Test(&req,&flag,&status);
			if ( flag == true ) { // we have something to do
				if (phase == 0 ) {
					//printf("%d %d\n",info[0],info[1]);
					phase=1;
					sieve_interval(info[0],info[1]);
					//printf("interval sieved!");std::cout.flush();
					//printf("%d\n",sieved_vectors);
					MPI_Request tmp_req;
					MPI_Isend(&sieved_vectors,1,MPI_INT,MY_ROOT,1,MPI_COMM_WORLD,&tmp_req); /*** inform root we finished the work ***/
					MPI_Irecv(&matrix_index,1,MPI_INT,MY_ROOT,2,MPI_COMM_WORLD,&req);
				} else if (phase == 1 ){
					MPI_Request tmp_req;
					MPI_Test(&req,&flag,&status);
					if ( flag == true){
						do_broadcast(matrix_index);
						int dummy=0;
						MPI_Isend(&dummy,1,MPI_INT,MY_ROOT,3,MPI_COMM_WORLD,&tmp_req);
						MPI_Irecv(info,2,MPI_INT,MY_ROOT,0,MPI_COMM_WORLD,&req);
						phase=0;
					}
				}
			}
			step_matrix_construction_protocol();	
		}
		finalize_matrix_communication();
		/*** FREE MEMORY ***/
		MPI_Cancel(&req);
		MPI_Request_free(&req);  // free the request
		while ( !q_b_data.empty()){
			b_data d=q_b_data.front();
			q_b_data.pop();
			delete [] d.data;
		}
	}

	void do_broadcast(int matrix_index){
		b_data b;
		b_info b_i;
		
		for ( unsigned int i=0;i<found_exponents_binary.size();i++)
			b_i.info[0]+=found_exponents_binary[i].non_zero_elements();
		for ( unsigned int i=0;i<found_exponents.size();i++){
			b_i.info[1]+=2*found_exponents[i].nnz(); // adding number of non-zero entries
		}
		b_i.info[1]+=found_exponents.size();
		b_i.info[0]+=found_exponents_binary.size();
		b_i.info[2]=found_smooth_numbers.size();
		b_i.info[3]=matrix_index;
		q_i_data.push(b_i);
		
		int sz=0; for ( int i=0;i<3;i++) sz+=b_i.info[i];
		sz+=b_i.info[1];
		b.data=new int[sz];
		int index=0;
		for ( unsigned int i=0;i<found_exponents_binary.size();i++){
			std::vector<int> dat=found_exponents_binary[i].getData();
			for ( std::vector<int> ::iterator it=dat.begin();it!=dat.end();it++){
				b.data[index++]=*it;
			}
			b.data[index++]=-1;
		}
		for ( unsigned int i=0;i<found_exponents.size();i++){
			for (unsigned int j=0;j<found_exponents[i].size();j++){
				if ( found_exponents[i][j]!=0 ){
					b.data[index++]=j;
					b.data[index++]=found_exponents[i][j];
				}
			}
			b.data[index++]=-1;
		}
		std::vector<int>::iterator it;
		for ( it=found_smooth_numbers.begin();it!=found_smooth_numbers.end();it++)
			b.data[index++]=*it;
		q_b_data.push(b);	
		b.sz=index ;	
		// sending the data
		MPI_Request tmp_req1,tmp_req2;
		for ( int i=0;i<mpi_size;i++){
				MPI_Isend(b_i.info,4,MPI_INT,i,4,MPI_COMM_WORLD,&tmp_req1);
				MPI_Isend(b.data,b.sz,MPI_INT,i,5,MPI_COMM_WORLD,&tmp_req2);
		}


		// clear vectors, so that they are ready for next sieving	
		found_exponents_binary.clear();
		found_exponents.clear();
		found_smooth_numbers.clear();	
	}


	void sieve_interval(int start, int end){
		cpp_int *U;			/*** Here I use more memory than needed, I can reduce it by subracting all the elements by sqrt(n) ***/
		cpp_int *V;			/*** Here I can not use less memory than the amount I am using now ***/
		/*** Additional remark: the elements size of array V will grow very fast(quadratic) ***/
		int signs[2]={-1,1};
		int sz=end-start;
		sieved_vectors=0;
		vec1 *exp_vectors= new vec1[sz];
		for ( int i=0;i<sz;i++)
			exp_vectors[i].resize(base_size+1);
		U=new cpp_int[sz];
		V=new cpp_int[sz];
		// I should change this part of code, because cpp_int is not of fixed size!!!
		for ( int cnt=1;cnt<2;cnt++){ // taking care of both signs
			for ( int i=0;i<sz;i++){
				U[i]=sqrt_n+signs[cnt]*(start+i);
				V[i]=U[i]*U[i]-n;
			}
			/*** allocating memory for results ***/
			/*** sieving ***/
			for ( int i=0;i<base_size;i++){
				if (base[i]==2 ) 
					continue;
				std::pair<int,int> residuals=get_quadratic_residuals(base[i]);
				int p=base[i];
				int r1=residuals.first;
				int r2=residuals.second;
				//process r1 first
				//cpp_int tmp_offset=(sqrt_n-u)%p;
				cpp_int tmp_offset=U[0]%p;
				int offset=tmp_offset.convert_to<int>();
				if (offset < 0 ) offset+=p;   // make sure offset is positive
				int first_index=r1-offset;
				if ( first_index < 0 ) 
					first_index+=p;
				for(int ind=first_index; ind<sz;ind+=p ){
					while ( V[ind]%p==0 ){
						V[ind]/=p;
						exp_vectors[ind](i)+=1;
					}	
				}
				// process r2				
				tmp_offset=U[0]%p;
				offset=tmp_offset.convert_to<int>();		 
				if (offset < 0 ) offset+=p;   // make sure offset is positive
				first_index=r2-offset;
				if ( first_index < 0 ) 
					first_index+=p;
				for(int ind=first_index; ind<sz;ind+=p ){
					while ( V[ind]%p==0 ){
						V[ind]/=p;
						exp_vectors[ind](i)+=1;
					}	
				}
			}

			// process 2 and -1
			int offset=0;
			if ( (V[0]&1)==1)
				offset=1;
			for ( int ind=offset;ind<sz;ind+=2){
				while ( (V[ind]&1 ) !=0 ) {
					V[ind]=(V[ind] >> 1);
					exp_vectors[ind](0)+=1;
				}
			}
			for ( int ind=0;ind<sz;ind++){
				if ( V[ind]<0 ){
					V[ind]=-V[ind];
					(exp_vectors[ind])(base_size)+=1;
				} else{
					break; // optimization
				}
			}
			for ( int i=0;i<sz;i++)
				if ( V[i] == 1 ){
					sieved_vectors++;
					//std::cout << "ORIGINALY SIEVED: " << U[i] <<std::endl;
					/*** inform about the number! ***/
					found_smooth_numbers.push_back(signs[cnt]*(start+i));
					found_exponents.push_back(exp_vectors[i]);
					found_exponents_binary.push_back(vec(exp_vectors[i]));
				}
			// free the memory 
		}	
		delete[] U;
		delete[] V;
		delete[] exp_vectors;
	}	

}
