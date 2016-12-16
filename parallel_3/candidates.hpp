#include "parameters.hpp"

namespace quadratic_sieve{

	/*** classical Gaussian algorithm ***/
	void find_choice_vector();


	/*** Block-Wiedemann algorithm-parallel ***/
	void p_find_choice_vector();

	void matrix_multiply(vec &y,std::vector<vec> &matrix, vec &x);
	int berlekamp_massey(vec &array, int sz, vec &c);

	bool try_factorization();
}
