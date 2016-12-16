#include "parameters.hpp"


struct b_info{
	int info[4];  /*** info 0 - size of binary chunk of matrix
									*       1 - size of full chunk of matrix
									*       2 - size of array of selected smooth numbers
									*       3 - position in a final matrix***/
	b_info(){
		for ( int i=0;i<4;i++)
			info[i]=0;
	}	
};

/*** data chunk we use to construct matrix ***/
struct b_data{
	int sz;
	int *data;
};


namespace quadratic_sieve{
	
	void sieving_coord();
	void step_matrix_construction_protocol();
	void init_matrix_communication();
	bool all_sieving_recieved();
	void process_data(b_info &info,b_data &data);

}
