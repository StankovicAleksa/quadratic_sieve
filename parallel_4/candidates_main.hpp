#include "parameters.hpp"

namespace quadratic_sieve{

	/*** Block-Wiedemann algorithm-parallel ***/
	void find_choice_vector();
	void parallel_matrix_multiply(vec &y,std::vector<vec> &matrix, vec &x);
	int berlekamp_massey(vec &array, int sz, vec &c);
	bool try_factorization();
}
