#include <iostream>
#include <string>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/integer.hpp>
#include "mod_arithm.hpp"
using boost::multiprecision::cpp_int;
using boost::multiprecision::mpz_int;
using namespace boost;
using namespace quadratic_sieve;
int main() {
	cpp_int n("34509701867117");
	int p=46649;
	cpp_int tmp_n=n%p;
	//std::cout<<tmp_n;
	tonelli_shanks(tmp_n.convert_to<int>(),p);
	return 0;
}
