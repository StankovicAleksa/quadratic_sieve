#include <iostream> 
#include <cstdio>
#include <vector>
#include <boost/numeric/ublas/vector_sparse.hpp>
using namespace boost::numeric::ublas;

int main(){
	for (unsigned int i=3;i<4;i++){
		//for (unsigned int i=0;i<exponents.size();i++){
		//std::vector <int > v=exponents_binary[i].getData();
		//for(int j=0;j<v.size();j++) printf("%d ",v[j]);
		//exponents[i].clear();
		int base_size=10;
		mapped_vector<double> v(base_size,base_size); 
		for(unsigned int j=0;j<v.size();++j) v(j)=1; 
		for(unsigned int j=0;j<v.size();++j) {
			std::cout << v(j) << " ";
			int a=v(j);
		 	printf("%d ",a);
		}
		printf("\n");
		//std::cout<<"num:"<<smooth_numbers[i]<<std::endl;
	}	
	return 0;
}
