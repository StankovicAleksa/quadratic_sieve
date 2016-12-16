/*** Place for all global variables ***/
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>
using boost::multiprecision::cpp_int;
using namespace boost::numeric::ublas;
typedef vector<int> vec;
namespace quadratic_sieve{  
	extern boost::multiprecision::cpp_int n; // the number to be factorized 
	extern int B;
	extern int *base;  /*** Prime factors base ***/
	extern int base_size;  /*** Size of the prime factors base ***/
	extern int matrix_length;   /*** Length of the matrix of binary vectors ***/
	extern int u;								/*** Parameter important for sieving phase. For more info look at the sc
															scientific report ***/

	extern std::vector<cpp_int> smooth_numbers;
	extern std::vector<vec> exponents;	/*** matrix containing vectors found so far ***/
	extern std::vector<vec> exponents_binary;
	extern std::vector<cpp_int> factors;
	extern int *lin_solution;
	extern int lin_sol_size;
}


