/*** Place for all global variables ***/

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <vector>
#include "sparse_binary_vector.hpp"

#define MY_ROOT 0
using boost::multiprecision::cpp_int;
using namespace boost::numeric::ublas;

typedef quadratic_sieve::SparseBinaryVector vec;
typedef mapped_vector<int> vec1;

namespace quadratic_sieve{  

	/*** GENERAL INFO ***/	
	extern boost::multiprecision::cpp_int n; // the number to be factorized 
	extern boost::multiprecision::cpp_int sqrt_n; // the number to be factorized 

	/*** PRE-SIEVING INFORMATIONS ***/
	extern int B;
	extern int *base;  /*** Prime factors base ***/
	extern int base_size;  /*** Size of the prime factors base ***/
	extern int matrix_length;   /*** Length of the matrix of binary vectors ***/
	extern int u;								/*** Parameter important for sieving phase. For more info look at the sc
															scientific report ***/
	extern std::vector<cpp_int> factors;


	/*** SIEVING INFO ***/
	extern std::vector<cpp_int> smooth_numbers;
	extern std::vector<vec1> exponents;	/*** matrix containing vectors found so far ***/
	extern std::vector<vec> exponents_binary;

	/*** INFORMATION FOR BLOCK-WIEDEMANN ALGORITHM ***/
	extern vec lin_solution;
	extern int lin_sol_size;

	/*** MPI STUFF ***/
	extern int mpi_rank;  // MPI rank
	extern int mpi_size;  // MPI size
}
