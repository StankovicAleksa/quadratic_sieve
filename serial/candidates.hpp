#include "parameters.hpp"

namespace quadratic_sieve{

	/*** classical Gaussian algorithm ***/
	void find_choice_vector();


	/*** Block-Wiedemann algorithm-parallel ***/
	void p_find_choice_vector();

	void matrix_multiply(int *y,std::vector<vec> &matrix, int *x);
	void vectors_multiply(int &c,int *a, int *b, int sz);
	int berlekamp_massey(int *array, int sz, int *c);
	void add_vectors( int * a, int *b, int sz);	

	bool try_factorization();
}
