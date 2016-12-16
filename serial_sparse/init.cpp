#include "parameters.hpp"
#include "init.hpp"
#include "mod_arithm.hpp"

namespace quadratic_sieve{
	int B;	
	int base_size=0;
	int *base;
	int u;
	boost::multiprecision::cpp_int n;
	std::vector<cpp_int> factors;

	/*** Input the number to be factored. The input is stored in 
	 * parameters.hpp variable n. ***/
	bool input_number(){
		//std::string str_n="72328663292947";
		std::string str_n="11869840501";		
		//std::string str_n="3530435723";
		//std::cin >> str_n;
		n=cpp_int(str_n); // the number to be factorized 
		return true;
	}

	/*** Set up parameters, like smoothness bound etc.
	 * The smoothness bound size will depend only on the size of n
	 * At the later stage we can also set up different polynomials for
	 * sieveing process(currently we are only using x^2-n.
	 ***/
	void setup_parameters(){  
	  B=5599;
    u=55500;

//		B=570000;
//		u=1058000;
	}

	/*** For generating all primes, we use Erathostenes sieve. At this stage we also use Tonelli-Shanks
	 * algorithm to check whether the n is quadratic residual modulo respective primes.
	 ***/
	void generate_factor_base(){
		bool *composite=(bool*)malloc(sizeof(bool)*(B+1));  // this can be optimized with bitmasks, but there is no need
																								// since this part is not a bottleneck of algorithm	
		for( int i=0;i<=B;i++)
			composite[i]=false;
		int sqrt_B=sqrt(B);
		for ( int i=2;i<sqrt_B;i++){
			if ( !composite[i] ) {
				int tmp_sum=i*2;  // change multiplication with addition
				while ( tmp_sum <= B){
					composite[tmp_sum]=true;
					tmp_sum+=i;
				}
			}	
		}	
		// prime numbers are generated
	

		// check if n is quadratic residual module candidates for prime numbers	
		// by calculating Legendre symbol
		//
		// We can parallelize this part
		
		bool *good_candidate=(bool*)malloc(sizeof(bool)*(B+1));  // this can be optimized with bitmasks, but there is no need
		for ( int i=0;i<=B;i++)
			good_candidate[i]=false;
		/*** check divisibility ***/
		for ( int i=2;i<=B;i++){
			if ( !composite[i] ){
				while ( n%i==0 ){
					n/=i;
					factors.push_back(i);
					printf("FACTOR %d\n",i);
				}
			}
		}
		std::cout<< n << std::endl;

		for ( int i=2;i<=B;i++){
			if ( !composite[i] && n_is_quadratic_residual(i) ){
				good_candidate[i]=true;	
				base_size++;
			}
		}
		base=(int*)malloc(base_size*sizeof(int));	
		int tmp=0;
		for ( int i=2;i<=B;i++){
			if( !composite[i] && good_candidate[i]) 
				base[tmp++]=i;	
			}
		free(composite);
		free(good_candidate);

	}
}
