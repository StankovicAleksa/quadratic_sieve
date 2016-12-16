#include "sieving.hpp"
#include "parameters.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include "mod_arithm.hpp"

namespace quadratic_sieve{

	std::vector<cpp_int> smooth_numbers;
	std::vector<vec> exponents;	/*** matrix containing vectors found so far ***/
	std::vector<vec> exponents_binary;
	vec **exp_vectors=0;
	int found_vectors_num;
	int max_vectors_size=base_size+10;  /*** we expect that after most 10 unsucessful trials of candidates a^2=b^2 (mod n)
																				* we will surely factorize the number 
																				***/


	/*** WARNING: here I will probably need distributed memory!!! ***/
	cpp_int **U;			/*** Here I use more memory than needed, I can reduce it by subracting all the elements by sqrt(n) ***/
	cpp_int **V;			/*** Here I can not use less memory than the amount I am using now ***/
									/*** Additional remark: the elements size of array V will grow very fast(quadratic) ***/
	void sieving_phase(){
		/*** constructing the arrays ***/
		cpp_int sqrt_n=sqrt(n);
		//std::cout << sqrt_n;
		U=(cpp_int**)malloc((2*u+1)*sizeof(cpp_int*));
		V=(cpp_int**)malloc((2*u+1)*2*sizeof(cpp_int*));
		// I should change this part of code, because cpp_int is not of fixed size!!!
		for ( int i=-u;i<=u;i++){
			cpp_int* a=(cpp_int*) malloc(sizeof(cpp_int));
			//V[i+u]=(cpp_int*) malloc(sizeof(cpp_int));
			//*a=sqrt_n+i;
			//(*(U[i+u]))=sqrt_n+i;
		//	*V[i+u]=(*U[i+u])*(*U[i+u]);
		//	*V[i+u]-=n;
		}
		return ;
		/*** allocating memory for results ***/
		exp_vectors=(vec**) malloc((2*u+1)*sizeof(vec*));
		for ( int i=0;i<2*u+1;i++){
			exp_vectors[i]=(vec*) malloc(sizeof(vec));
			exp_vectors[i]->resize(base_size+1);
		}
		/*** sieving ***/
		for ( int i=0;i<base_size;i++){
			if (base[i]==2 ) 
				continue;
			std::pair<int,int> residuals=get_quadratic_residuals(base[i]);
			int p=base[i];
			int r1=residuals.first;
			int r2=residuals.second;

			//process r1 first
			cpp_int tmp_offset=(sqrt_n-u)%p;
			int offset=tmp_offset.convert_to<int>();		 
			if (offset < 0 ) offset+=p;   // make sure offset is positive
			int first_index=r1-offset;
			if ( first_index < 0 ) 
				first_index+=p;
			for(int ind=first_index; ind<=2*u;ind+=p ){
				while ( *V[ind]%p==0 ){
					*V[ind]/=p;
					(*exp_vectors[ind])(i)+=1;
				}	
			}
			// process r2				
			tmp_offset=(sqrt_n-u)%p;
			offset=tmp_offset.convert_to<int>();		 
			if (offset < 0 ) offset+=p;   // make sure offset is positive
			first_index=r2-offset;
			if ( first_index < 0 ) 
				first_index+=p;
			for(int ind=first_index; ind<=2*u;ind+=p ){
				while ( *V[ind]%p==0 ){
					*V[ind]/=p;
					(*exp_vectors[ind])(i)+=1;
				}	
			}
		}

		// process 2 and -1
		
		int offset=0;
		if ( (*V[0]&1)==1)
			offset=1;
		for ( int ind=offset;ind<=2*u;ind+=2){
			while ( (*V[ind]&1 ) !=0 ) {
				*V[ind]=(*V[ind] >> 1);
				(*exp_vectors[ind])(0)+=1;
			}
		}
		for ( int ind=0;ind<=2*u;ind++){
			if ( V[ind]<0 ){
				*V[ind]=-*V[ind];
				(*exp_vectors[ind])(base_size)+=1;
			} else{
				break; // optimization
			}
		}
		// sieveing is finished
	
//		for ( int i=-u;i<=u;i++)
//			std::cout<< V[i+u] << std::endl;

		// initialize matrix for solving
		for ( int i=0;i<=2*u;i++)
			if ( *V[i] == 1 ){
				smooth_numbers.push_back(*U[i]);		
				exponents.push_back(*exp_vectors[i]);
				exponents_binary.push_back(binary(*exp_vectors[i]));
			}		
		// free the memory 
	}	
}
