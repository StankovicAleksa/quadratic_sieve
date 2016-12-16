#include "mod_arithm.hpp"
#include <stack>

#include "parameters.hpp"
#include <boost/multiprecision/cpp_int.hpp>


using boost::multiprecision::cpp_int;

namespace quadratic_sieve{

	//------------------------------------------------------------------------------------------
	// Fast exponentation procedure, without recursion
	//------------------------------------------------------------------------------------------
	/*** Returns positive number ***/
	int power(int a,int b,int mod){
		// some checks

		long long res=1;  // initial value
		a%=mod; 
		if ( a < 0 ) 
			a+=mod;
		// here we avoid recursion through the use of stack

		long long tmp=a;
		while ( b!=0 ){
			if ( (b&1) != 0 ){
				res*=tmp;
				res%=mod;
			}
			tmp*=tmp;
			tmp%=mod;			
			b=(b>>1);
		}
		return res;
	}


	//------------------------------------------------------------------------------------------
	
	/*** Returns positive number ***/
	int power(int a,cpp_int b,int mod){
		// some checks

		long long res=1;  // initial value
		a%=mod; 
		if ( a < 0 ) 
			a+=mod;
		// here we avoid recursion through the use of stack

		long long tmp=a;
		while ( b!=0 ){
			if ( (b&1) != 0 ){
				res*=tmp;
				res%=mod;
			}
			tmp*=tmp;
			tmp%=mod;			
			b=(b>>1);
		}
		return res;
	}


	//------------------------------------------------------------------------------------------
	
	/*** Returns positive number ***/
	cpp_int power(cpp_int a,int b,cpp_int mod){
		// some checks

		cpp_int res=1;  // initial value
		a=a%mod; 
		if ( a < 0 ) 
			a=a+mod;
		// here we avoid recursion through the use of stack

		cpp_int tmp=a;
		while ( b!=0 ){
			if ( (b&1) != 0 ){
				res*=tmp;
				res%=mod;
			}
			tmp*=tmp;
			tmp%=mod;			
			b=(b>>1);
		}
		return res;
	}


	//------------------------------------------------------------------------------------------
	// We call underlying function that calculates Legendre symbol
	//------------------------------------------------------------------------------------------
	bool n_is_quadratic_residual(int prime){ 
		cpp_int tmp=(n % cpp_int(prime));
		int n_i=tmp.convert_to<int>();
		return is_quadratic_residual(n_i,prime);
	}

	//------------------------------------------------------------------------------------------
	// We just calculate Legendre symbol with fast exponentation
	//------------------------------------------------------------------------------------------
	bool is_quadratic_residual(int n, int prime){		
		int leg_symbol=power(n,(prime-1)/2,prime);
		if ( leg_symbol== 1) 
			return true;
		else if ( leg_symbol == 0) {
			printf("SHOULD NOT HAPPEN");
			return false;
		}
		else
			return false;
	}	
	//------------------------------------------------------------------------------------------
	// Here is the implementation of Tonelli-Shanks algorithm
	//------------------------------------------------------------------------------------------
	std::pair<int,int> get_quadratic_residuals(int p){
		cpp_int tmp=n%p;
		int n_tmp=tmp.convert_to<int>();	
		int sol_1=tonelli_shanks(n_tmp,p);
		int sol_2=p-sol_1;
		return std::make_pair(sol_1,sol_2);
	}

	int tonelli_shanks(int n,int p){
		int Q=p-1;
		int S=0;
		while ((Q&1)==0){
			Q=(Q>>1);
			S++;
		}
		if ( p%4 ==3  ){
			return power(n,(p+1)/4,p);
		}
		/*** find z such that {z \choose p} = 1
		 * We do this by simple search from 2 to p. This can be optimized
		 * but will not significantly impact speed of algorithm 
		 ***/ 
		int z;
		for ( z=2;z<p;z++){
			if (!is_quadratic_residual(z,p) )  /** we will surely find at least one that is not quadratic residual **/
				break;
		}

		/*** Check here for overflows!!! ***/
		long long c=power(z,Q,p);
		long long R=power(n,(Q+1)/2,p);	
		long long t=power(n,Q,p);
		long long M=S;
		while ( t!= 1){
			long long tmp_t=t;
			int i;
			for ( i=1;i<M;i++){
				tmp_t*=tmp_t;
				tmp_t%=p;
				if ( tmp_t==1 ) 
					break;
			}
			cpp_int tmp_exp=1;
			for ( int j=0;j<M-i-1;j++)
				tmp_exp=(tmp_exp << 1);
			long long b=power(c,tmp_exp,p);
			//std::cout << tmp_exp << " " << b << std::endl;
			//if ( b==1 )
				//break;
			R=R*b%p;	
			t=t*b%p;
			t=t*b%p;
			c=b*b%p;
			M=i;
		//	std::cout << b << " " << R << " " << t << " " << c << " " << M << std::endl;  	
		}
		return R;
	}

	/*** Greatest common divisor, iterative ***/
	cpp_int gcd(cpp_int a, cpp_int b){
		cpp_int c;	
		// set the signs
		if ( a< 0 ) a=-a;
		if (b < 0 ) b=-b;
		
		// set a> b
	
		if ( a< b) {
			c=a;
			a=b;
			b=c;
		}

		// find gcd, again with iterative method(avoid recursion)
	
		while ((c=a%b)!=0){
			a=b;b=c;
		}	
		return b;
	}


	//------------------------------------------------------------------------------------------
	// Matrix multiplication
	//------------------------------------------------------------------------------------------


	/*** Multiply matrix with a vector x, save result in y ***/
	void matrix_multiply(vec &y,std::vector<vec> &matrix,vec &x){
		//int width=matrix.size();
		//int height=matrix[0].size();
		std::vector<int>::iterator it_x;
		std::vector<int> &x_data=x.getData();
		y.clear(); y.resize(base_size+1); // make it 0

		for ( it_x=x_data.begin();it_x!=x_data.end();it_x++){
			int j=*it_x;
			y+=matrix[j];
		}
	}	

}
