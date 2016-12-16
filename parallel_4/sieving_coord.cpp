#include "sieving_coord.hpp"
#include "sieving.hpp"
#include <mpi.h>


namespace quadratic_sieve{
	/*** variables for communication for coordinator ***/
	int *process_phase; // 0 - free, 1 - sieving, 2 - preparing broadcast, 3 - coordinator
	int *sieved_amount; // denotes the cardinality of smooth number in interval given to a node

	/*** Variables for synchronisation ***/
	int current_siev_index=1;   // current location of pointer in sieving process, pointer is smaller than smoothness bound B
	int cur_matrix_width=0;	    // width of the matrix built so far(matrix required for the final stage

	/*** Logical variable ***/
	bool all_sieving_finished=false;

	/*** Request variables ***/
	MPI_Request final_send_bcast;
	MPI_Request* sieving_reqs;  // reports from processes about finishing their sieving
	MPI_Request* availabl_reqs;  // reports from processes about finishing broadcasting

	/*** Configuration variables ***/
	int SIEV_CHUNK_SIZE=100000;  // Denotes chunk size we will send for sieving

	void sieving_coord(){
		init_matrix_communication();		
		/*** Initialization part ***/
		process_phase=new int[mpi_size];  
		sieved_amount=new int[mpi_size];  
		sieving_reqs=new MPI_Request[mpi_size];
		availabl_reqs=new MPI_Request[mpi_size];

		/*** Intialize all processes except the root to be free for work ***/
		for ( int i=0;i<mpi_size;i++){
			process_phase[i]=0;
		}
		process_phase[mpi_rank]=3;


		/*** We do all the tasks of sieveing communication in this loop ***/
		while ( !all_sieving_recieved() ){
			send_sieving_tasks();  // send tasks for sieving of an interval
			process_sieve_reports();  // process reports, giving the information where the data should be put in resulting matrix
			step_matrix_construction_protocol();  // protocol required for construction of the final matrix
		}
		finalize_matrix_communication();
		/*** Cleanup of memory ***/
		delete[]process_phase;
		delete[] sieved_amount;
		delete[] sieving_reqs;
		delete[] availabl_reqs;

	}

	/*** SEND THE TASKS TO ALL THE PROCESSES ***/
	void send_sieving_tasks(){
		/*** IF possible, send sieveing tasks to free processess ***/
		if ( current_siev_index<u ){ 
			for ( int i=0;i<mpi_size;i++){
				if ( process_phase[i]==0 ) {
					int info[2]={current_siev_index,current_siev_index+SIEV_CHUNK_SIZE};
					if ( info[1]>u ) 
						info[1]=u;
					MPI_Request r;
					MPI_Isend(info,2,MPI_INT,i,0,MPI_COMM_WORLD,&r); // send request for sieving
					MPI_Irecv(sieved_amount+i,1,MPI_INT,i,1,MPI_COMM_WORLD,sieving_reqs+i); // ask for result
					current_siev_index+=SIEV_CHUNK_SIZE;
					process_phase[i]=1;
					if ( current_siev_index>u ) /*** we are over the boundary ***/
						break;
				}	
			}
		}
		if ( current_siev_index > u && !all_sieving_finished){
			bool temp_bool=true;
			for ( int i=0;i<mpi_size;i++){
				if ( process_phase[i]==1 ){
					temp_bool=false;
					break;
				}
			}
			if ( temp_bool ){
				all_sieving_finished=true;
				MPI_Ibcast(&cur_matrix_width,1,MPI_INT,MY_ROOT,MPI_COMM_WORLD,&final_send_bcast);
				set_final_matrix_width(cur_matrix_width);
			}
		}
		int flag;MPI_Status status;

		// I think I don't need this
		if ( all_sieving_finished)
			MPI_Test(&final_send_bcast,&flag,&status);
	}

	/*** Coordination of the sieving process, MUST be called by root ***/
	void process_sieve_reports(){
		// local temporary variables
		int flag;
		MPI_Status status;
		int dummy;

		// main loop
		for ( int i=0;i<mpi_size;i++){
			/*** Node is working ***/
			if ( process_phase[i]==1 ){
				MPI_Test(sieving_reqs+i,&flag,&status);
				if( flag==true  ){
				 	//MPI_Request_free(sieving_reqs+i);	
					MPI_Request req_temp;
					MPI_Isend(&cur_matrix_width,1,MPI_INT,i,2,MPI_COMM_WORLD,&req_temp);	
				  MPI_Irecv(&dummy,1,MPI_INT,i,3,MPI_COMM_WORLD,availabl_reqs+i);
					cur_matrix_width+=sieved_amount[i]; // add the number off sieved 
					process_phase[i]=0;  // mark the process as broadcasting the information
				}
			} else if ( process_phase[i]==2 ){
				MPI_Test(availabl_reqs+i,&flag,&status);
				if (flag == true ){
					//MPI_Request_free(availabl_reqs+i);
					process_phase[i]=0;  // mark the process as free again! 
				}
			}
		}
	}
}
