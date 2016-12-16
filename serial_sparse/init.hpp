#include <iostream>
#include <string>
#include "parameters.hpp"
namespace quadratic_sieve{

	extern int B;  // Smoothness bound

	/*** Input the number n to be factored 
	 *	 The function returns true if the number is sucessfuly retrieved from input source.
	 *	 On failure, it returns false
	 ***/
	bool input_number();

	/*** Setting parameters for execution, like smoothness bound etc. 
	 ***/
	void setup_parameters();

	/*** Generate the base of all primes p_i smaller than the smoothness bound B
	 *   and which satisfy that the n is quadratic residual modulo p_i
	 ***/
	void generate_factor_base();
}
