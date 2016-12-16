#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include "sparse_binary_vector.hpp"
#include "parameters.hpp"
#include "candidates.hpp"
#include "mod_arithm.hpp"
using namespace quadratic_sieve;
namespace quadratic_sieve{

	/*** helper function for switching between mapped_vector and compressed_vector ***/
	//inline void push_back(vec &a, int index, int el){
	//	a[index]=el;
	//}

	void init_seed(){

		/*** MPI stuff ***/	
		int rank,size,nameLen;
		char processorName[MPI_MAX_PROCESSOR_NAME];
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Get_processor_name(processorName,&nameLen);
		/*** take different seed for each process ***/
		srand((unsigned)time(NULL)+rank*size + rank+ nameLen);
	}
	vec lin_solution;
	int lin_sol_size=0;

	inline void random_vec(vec &x, int size){
		for ( int i=0;i<size;i++){
			if ( (rand()&1) ) 
				x.push_back(i);	
		}
	}

	inline void calc_sequence(vec &y, std::vector<vec> &matrix, vec &x, int len,int iterations,vec &sequence){
		vec tmp_vec; tmp_vec.resize(len);
		vec t_vec;   t_vec.resize(len);
		tmp_vec = x;
		int tmp_num;
		sequence.clear();
		for ( int i=0;i<iterations;i++){
			tmp_num=y*tmp_vec;
			if ( tmp_num == 1)
				sequence.push_back(i);
			t_vec.clear();
			t_vec=tmp_vec;
			matrix_multiply(tmp_vec,exponents_binary,t_vec);
		}
	}
	void find_choice_vector(){
		init_seed();
		if ( exponents_binary.size()==0) {
			printf("More sieving is required(0 numbers in the sieve)!!!\n");
			return;
		}
		int height=exponents_binary[0].size();
		int width=exponents_binary.size();

		if ( width <= height) {
			printf("More sieving is required!!!\n");
			return;
		}
		vec sequence;    sequence.resize(2*width+1);
		vec tmp_vec; 		 tmp_vec.resize(width);
		vec t_vec; 		   t_vec.resize(width);
		vec c; 					 c.resize(2*width+1);
		vec x_base;		   x_base.resize(width);
		vec x;      		 x.resize(width);
		vec y;					 y.resize(width);
		std::cout.flush();
		while (true){
			x+=y;
			random_vec(x_base,width);
			/** define x=M*x_base */
			matrix_multiply(x,exponents_binary,x_base);
			/** generate x_base vector randomly **/
			printf("MATRIX CALCULATED");
			std::cout.flush();
			/** generate vector y randomly **/
			random_vec(y,width);

			/*** tries to find non-zero kernel ***/	

			/*** calculate sequence ***/
			calc_sequence(y,exponents_binary,x,width,2*width,sequence);
			int l=berlekamp_massey(sequence,2*width,c);
			/*** try to find non-zero vector ***/
			vec res; res.resize(width);
			tmp_vec=x_base;
			std::vector < int > data=c.getData();

			// *** very careful here *** 
			std::vector<int>::reverse_iterator it=data.rbegin();
			for ( int i=l;i>=0;i--){ 
				while ( it!= data.rend() && *it>i ) 
					it++;
				if ( it!=data.rend() && *it==i ){
					//if ( c[i]!=0){
					res+=tmp_vec;
				}
				//}
				t_vec=tmp_vec;
				matrix_multiply(tmp_vec,exponents_binary,t_vec);
			}

			matrix_multiply(tmp_vec,exponents_binary,res);	
			if ( tmp_vec.non_zero_elements() ==0 ){
				// try factorization!!!
				lin_solution=res;		
				lin_sol_size=width;
				if (try_factorization()){
					break;
				}
			}
		}	
	}
	void p_find_choice_vector(){
	}
	bool try_factorization(){
		//int width=exponents_binary.size();
		//int height=exponents_binary[0].size();
		cpp_int a=1;
		cpp_int b=1;
		std::vector <int > linData=lin_solution.getData();
		std::vector < int > :: iterator it=linData.begin();
		for (;it!=linData.end();it++){
			a=a*smooth_numbers[*it]%n;
		}
		int *coeff=new int[base_size];
		for ( int i=0;i<base_size;i++)
			coeff[i]=0;

		it=linData.begin();
		for (;it!=linData.end();it++){
			for ( int j=0;j<base_size;j++){
				coeff[j]+=exponents[*it][j];
			}
		}
		for ( int i=0;i<base_size;i++){
			coeff[i]/=2;
			cpp_int p=base[i];
			b=b*power(p,coeff[i],n)%n;
		}
		cpp_int try_1=gcd(a+b,n);
		cpp_int try_2=gcd(a-b,n);
		//std::cout<<n<<std::endl;
		std::cout<< a << " "<< b<< std::endl;
		//std::cout<< a*a%n << " " << b*b%n << std::endl;
		std::cout << try_1 << " "<< try_2<<std::endl<<std::endl;
		if ( try_1!=1 && try_1!=n )
			return true;
		if ( try_2!=1 && try_2!=n )
			return true;
		return false;
	}


	/*** Multiply matrix with a vector x, save result in y ***/
	void matrix_multiply(vec &y,std::vector<vec> &matrix,vec &x){
		//int width=matrix.size();
		//int height=matrix[0].size();
		std::vector<int>::iterator it_x;
		std::vector<int> &x_data=x.getData();
		y.clear();  // make it 0

		for ( it_x=x_data.begin();it_x!=x_data.end();it_x++){
			int j=*it_x;
			y+=matrix[j];
		}
	}	

	void arraycopy(int *src,int src_pos,int *dest, int dest_pos, int length){
		for ( int i=0;i<length;i++)
			dest[dest_pos+i]=src[src_pos+i];
	}
	int berlekamp_massey(vec &v_array, int sz, vec &v_c){
		int *c=new int[sz];
		int *array=new int[sz];
		int *b=new int[sz];
		int *t=new int[sz];


		//vec b; b.resize(sz);
		//vec t; t.resize(sz);
		for ( int i=0;i<sz;i++)
			array[i]=b[i]=t[i]=c[i]=0;
		std::vector<int> &array_data=v_array.getData();
		for (	std::vector<int>::iterator it=array_data.begin();it!=array_data.end();it++)
			array[*it]=1;
		b[0] = 1;
		c[0] = 1;
		int l = 0;
		int m = -1;
		for (int n = 0; n < sz; n++) {
			int d = 0;
			for (int i = 0; i <= l; i++) {
				d = (d^(c[i] * array[n - i]));
			}
			if (d == 1) {
				arraycopy(c, 0, t, 0, sz);
				int sz_M = n - m;
				for (int j = 0; j < sz - sz_M; j++) {
					c[sz_M + j] = (c[sz_M+j]^b[j]);
				}
				if (l <= n / 2) {
					l = n + 1 - l;
					m = n;
					arraycopy(t, 0, b, 0, sz);
				}
			}
		}
		v_c.clear();
		for ( int i=0;i<=l;i++)
			if ( c[i]==1 ){
				v_c.push_back(i);
			}
		return l;
	}
}
