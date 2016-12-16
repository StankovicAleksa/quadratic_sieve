#include "candidates.hpp"
#include "parameters.hpp"
#include "mod_arithm.hpp"
#include <cstdio>
#include <cstdlib>
#include <ctime>
namespace quadratic_sieve{
	int *lin_solution=0;
	int lin_sol_size=0;

	inline void random_vec(int *vec, int size){
		for ( int i=0;i<size;i++){
			if ( (rand()&1) ) vec[i]=0;
			else vec[i]=1;
		}
	}

	inline void calc_sequence(int * y, std::vector<vec> &matrix, int *x, int len,int iterations,int *sequence){
		int *tmp_vec	=(int*) malloc( len*sizeof(int));
		int *t_vec=(int*) malloc(len*sizeof(int));

		for ( int i=0;i<len;i++)
			tmp_vec[i]=x[i];

		for ( int i=0;i<iterations;i++){
			vectors_multiply(sequence[i],y,tmp_vec,len);
			for ( int j=0;j<len;j++)
				t_vec[j]=tmp_vec[j];
			matrix_multiply(tmp_vec,exponents_binary,t_vec);
		}
		free(tmp_vec);
		free(t_vec);
	}
	void find_choice_vector(){
		int height=exponents_binary[0].size();
		int width=exponents_binary.size();

		if ( width <= height) {
			printf("More sieving is required!!!\n");
			return;
		}
		int * sequence=(int*) malloc((2*width+1)*sizeof(int));
		int *tmp_vec	=(int*) malloc(   width   *sizeof(int));
		int *t_vec    =(int*) malloc(   width   *sizeof(int));
		int *c        =(int*) malloc((2*width+1)*sizeof(int));
		int* x_base   =(int*) malloc(   width   *sizeof(int));
		int* x			  =(int*) malloc(   width   *sizeof(int));
		int* y        =(int*) malloc(   width   *sizeof(int));

		while (true){
			/** generate x_base vector randomly **/
			random_vec(x_base,width);
			/** define x=M*x_base */
			matrix_multiply(x,exponents_binary,x_base);

			/** generate vector y randomly **/
			random_vec(y,width);

			/*** tries to find non-zero kernel ***/	

			/*** calculate sequence ***/
			clock_t start=clock();
			calc_sequence(y,exponents_binary,x,width,2*width,sequence);
			clock_t end=clock();
      printf("TIME is : %lf\n",(float)(end-start)/CLOCKS_PER_SEC); 

			int l=berlekamp_massey(sequence,2*width,c);
			/*** try to find non-zero vector ***/
			int *res	=(int*) malloc( width*sizeof(int));

			for ( int i=0;i<width;i++){
				tmp_vec[i]=x_base[i];
				res[i]=0;
			}
			for ( int i=l;i>=0;i--){ // CHANGED!!!
				if ( c[i]!=0){
					add_vectors(res,tmp_vec,width);	
				}
				for ( int j=0;j<width;j++)
					t_vec[j]=tmp_vec[j];
				matrix_multiply(tmp_vec,exponents_binary,t_vec);
			}

			matrix_multiply(tmp_vec,exponents_binary,res);	
			bool is_kernel=true;
			for ( int i=0;i<width;i++){
				if ( tmp_vec[i] !=0 ){
					is_kernel=false;
					break;
				}
			}
			if ( is_kernel){
				// try factorization!!!
				lin_solution=res;		
				lin_sol_size=width;
				if (try_factorization()){
					break;
				}
			}
		}	
		free(sequence);
		free(tmp_vec);
		free(t_vec);
		free(c);
		free(x_base);
		free(x);
		free(y);

		}
	/*** a=a+b  ***/
	void add_vectors(int * a, int * b, int sz){
		for ( int i=0;i<sz;i++)
			a[i]^=b[i];
	}

	void p_find_choice_vector(){
	}
	bool try_factorization(){
		//int width=exponents_binary.size();
		//int height=exponents_binary[0].size();
		cpp_int a=1;
		cpp_int b=1;

		for ( int i=0;i<lin_sol_size;i++){
			if ( lin_solution[i] == 1 )
			a=a*smooth_numbers[i]%n;
		}
		int *coeff=(int*)malloc(base_size*sizeof(int));
		for ( int i=0;i<base_size;i++)
			coeff[i]=0;

		for ( int i=0;i<lin_sol_size;i++){
			if (lin_solution[i]!=0){
				for ( int j=0;j<base_size;j++){
					coeff[j]+=exponents[i][j];
				}
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
	void matrix_multiply(int *y,std::vector<vec> &matrix, int *x){
		int width=matrix.size();
		int height=matrix[0].size();
		for ( int i=0;i<width;i++)
			y[i]=0;
		for ( int i=0;i<height;i++){
			for ( int j=0;j<width;j++){
				if ( x[j] !=0 ){
					y[i]+=matrix[j][i];
					y[i]%=2;
				}
			}
		}
	}	

	/*** implement multiplication of vectors, c=a*b;  ***/
	void vectors_multiply(int &c,int *a, int *b, int sz){
		c=0;
		for ( int i=0;i<sz;i++)
			c^=(a[i]&b[i]);
	}
	void arraycopy(int *src,int src_pos,int *dest, int dest_pos, int length){
		for ( int i=0;i<length;i++)
			dest[dest_pos+i]=src[src_pos+i];
	}

	int berlekamp_massey(int *array, int sz, int *c){
		int b[sz];
		int t[sz];
		for ( int i=0;i<sz;i++)
			b[i]=t[i]=c[i]=0;
		b[0] = 1;
		c[0] = 1;
		int l = 0;
		int m = -1;
		for (int n = 0; n < sz; n++) {
		    int d = 0;
		    for (int i = 0; i <= l; i++) {
		        d ^= (c[i] * array[n - i]);
		    }
		    if (d == 1) {
		        arraycopy(c, 0, t, 0, sz);
		        int sz_M = n - m;
		        for (int j = 0; j < sz - sz_M; j++) {
		            c[sz_M + j] ^= b[j];
		        }
		        if (l <= n / 2) {
		            l = n + 1 - l;
		            m = n;
		            arraycopy(t, 0, b, 0, sz);
		        }
		    }
		}
		return l;
  }
}
