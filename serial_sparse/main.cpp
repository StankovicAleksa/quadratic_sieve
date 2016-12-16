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
	clock_t start=clock();
	input_number(); // input the number
	setup_parameters();  // setup the environment parameters
	generate_factor_base();
	printf("Sieving\n");
	sieving_phase();	
	printf("Searching\n");
	find_choice_vector();
	clock_t end=clock();
	printf("TIME is : %lf\n",(float)(end-start)/CLOCKS_PER_SEC);
	return 0;
}
