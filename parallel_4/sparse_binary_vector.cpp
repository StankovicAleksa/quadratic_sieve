#include "sparse_binary_vector.hpp"
#include <cstdio>

namespace quadratic_sieve{
	SparseBinaryVector::SparseBinaryVector(){}
	int SparseBinaryVector::size(){return length;}
	int SparseBinaryVector::non_zero_elements(){return data.size();}

	SparseBinaryVector::SparseBinaryVector(mapped_vector<int> &x){
		for ( mapped_vector<int>::iterator it=x.begin();it!=x.end();it++){
			if ( (*it%2) != 0 )
				data.push_back(it.index());
		}
		length=x.size();
	}		


	void SparseBinaryVector::resize(int n){length=n;}

	void SparseBinaryVector::push_back(int index){data.push_back(index);length++;}
	void SparseBinaryVector::clear(){ data.clear(); }
	std::vector	<int>& SparseBinaryVector::getData(){return data;}

	SparseBinaryVector& SparseBinaryVector::operator +=(SparseBinaryVector & other){
		std::vector < int > newData;	
		std::vector<int> &a=this->getData();
		std::vector<int> &b=other.getData();
		std::vector<int>::iterator it_a=a.begin();
		std::vector<int>::iterator it_b=b.begin();
		while ( it_a != a.end() || it_b!=b.end()){
			if ( it_a==a.end() ){
				newData.push_back(*it_b);
				it_b++;
			}
			else if ( it_b==b.end() ){
				newData.push_back(*it_a);
				it_a++;
			} 
			else { 
				if ( *it_a<*it_b){
					newData.push_back(*it_a);
					it_a++;
				} 
				else if ( *it_b<*it_a){
					newData.push_back(*it_b);
					it_b++;
				} else{
					it_a++;
					it_b++;
				}
			}
		}
		this->data=newData;
		return *this;	
	}
	int SparseBinaryVector::operator *(SparseBinaryVector &other){
	/*** implement multiplication of vectors, c=a*b;  ***/
		int c=0;
		std::vector<int> &a=this->getData();
		std::vector<int> &b=other.getData();	 
		std::vector<int>::iterator it_a=a.begin();
		std::vector<int>::iterator it_b=b.begin();
		while ( it_a != a.end() && it_b!=b.end()){
			if ( *it_a < *it_b)
				it_a++;
			else if ( *it_b < *it_a)
				it_b++;
			else{
				c^=1;
				it_a++;
				it_b++;
			}
		}
		return c;
	}
}
