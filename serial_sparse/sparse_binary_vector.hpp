#ifndef SPARSE_BIN_VECTOR
#define SPARSE_BIN_VECTOR
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <vector>

/*** Sparse vector containing binary operations ***/
using namespace boost::numeric::ublas;


namespace quadratic_sieve{

	class SparseBinaryVector{
		public:
			void resize(int n);
			SparseBinaryVector(mapped_vector<int> &x);
			SparseBinaryVector();
			void push_back(int i);  // sets bit i to 1, we assume that all bits at position >=i
			// are 0
			void clear();
			std::vector<int>& getData();
			int size();
			int non_zero_elements();
			SparseBinaryVector& operator+=(SparseBinaryVector &other);
			int operator*(SparseBinaryVector &other);
		private:
			std::vector<int> data;
			int length;
	};
}
#endif
