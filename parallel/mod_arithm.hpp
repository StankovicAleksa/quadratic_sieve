#include "parameters.hpp"
#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

namespace quadratic_sieve {

	/*** finds a^b modulo mod ***/
	int power(int a,int b, int mod);
	int power(int a,cpp_int b, int mod);
	cpp_int power(cpp_int a,int b, cpp_int mod);


	bool is_quadratic_residual(int n, int prime);
	/*** Return true if n is quadratic residual modulo p, false otherwise. ***/
	bool n_is_quadratic_residual(int prime);

	/*** get quadratic residuals ***/
	std::pair<int,int> get_quadratic_residuals(int p);	

	/*** Tonelli-Shanks algorithm ***/
	int tonelli_shanks(int n, int p);

	vec binary(vec1 &v);

	cpp_int gcd(cpp_int a, cpp_int b);
}
