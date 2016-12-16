#include <iostream>
#include <string>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/integer.hpp>
using boost::multiprecision::cpp_int;
using boost::multiprecision::mpz_int;
int main() {
	clock_t start=clock();
	std::string str_n="72328663292947";
	//std::cin >> str_n;
	cpp_int n=cpp_int(str_n); // the number to be factorized 
	//cpp_int i("8504299");
	cpp_int i=2;
	while ( true){
		if ( n%i==0 )
			break;
	}
	clock_t end=clock();
	printf("TIME is : %lf\n",(float)(end-start)/CLOCKS_PER_SEC);
	//std::cout<<tmp_n;
	return 0;
}
