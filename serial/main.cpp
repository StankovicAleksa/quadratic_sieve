#include <iostream>
#include <string>
#include "init.hpp"
#include "parameters.hpp"
#include "sieving.hpp"
#include "candidates.hpp"

/*** BE CAREFUL WITH VALUES 2 AND -1 ***/

using namespace quadratic_sieve;
boost::multiprecision::cpp_int n; // the number to be factorized 

int main(int args,char *argv[]) {
	input_number(); // input the number
	setup_parameters();  // setup the environment parameters
	generate_factor_base();
//
	//for ( int i=0;i<base_size;i++)
		//printf("%d\n",base[i]);
	
	sieving_phase();	
	find_choice_vector();
	return 0;
}
