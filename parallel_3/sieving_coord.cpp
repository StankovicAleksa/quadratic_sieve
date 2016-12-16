#include "sieving_coord.hpp"
#include "sieving.hpp"
#include <mpi.h>

/*** THINK AT THIS PLACE ABOUT USING STATIC BLOCKS! ***/

namespace quadratic_sieve{


	/*** variables for communication from coordinator ***/
	int *process_phase; // 0 - free, 1 - sieving, 2 - preparing broadcast, 3 - coordinator
	int *sieved_amount; // denotes the cardinality of smooth number in interval given to a node
	
	int current_siev_index=0;
	int cur_matrix_width=0;
	bool all_sieving_finished=false;

	MPI_Request final_send_bcast;
	MPI_Request* sieving_reqs;  // reports from processes about finishing their sieving
	MPI_Request* availabl_reqs;  // reports from processes about finishing broadcasting
	
	int SIEV_CHUNK_SIZE=100000;

	void sieving_coord(){
		init_matrix_communication();		
		/*** INITIALIZATION PART ***/
		process_phase=new int[mpi_size];  
		sieved_amount=new int[mpi_size];  
		sieving_reqs=new MPI_Request[mpi_size];
		availabl_reqs=new MPI_Request[mpi_size];

		for ( int i=0;i<mpi_size;i++){
			process_phase[i]=0;
		}
		process_phase[mpi_rank]=3;


		/*** We do all the tasks of sieveing communication in this loop ***/
		while ( !all_sieving_recieved() ){
			send_sieving_tasks();
			process_sieve_reports();
			step_matrix_construction_protocol();
		}
	
		/*** Here we need barrier, to synchronize all processes.
			* Barrier is not really required, but it is nice to have it, so that we know we go into 
			* the next phase together! 
			* ***/


	}

	/*** SEND THE TASKS TO ALL THE PROCESSES ***/
	void send_sieving_tasks(){
		/*** IF possible, send sieveing tasks to free processess ***/
		if ( current_siev_index<u ){ 
			for ( int i=0;i<mpi_size;i++){
				if ( process_phase[i]==0 ) {
					int info[2]={current_siev_index,current_siev_index+SIEV_CHUNK_SIZE};
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
				printf("FINISH!!");
				MPI_Ibcast(&cur_matrix_width,1,MPI_INT,MY_ROOT,MPI_COMM_WORLD,&final_send_bcast);
				set_final_matrix_width(cur_matrix_width);
			}
		}
		int flag;MPI_Status status;
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
					MPI_Request req_temp;
					MPI_Isend(&cur_matrix_width,1,MPI_INT,i,2,MPI_COMM_WORLD,&req_temp);	
				  MPI_Irecv(&dummy,1,MPI_INT,i,3,MPI_COMM_WORLD,availabl_reqs+i);
					cur_matrix_width+=sieved_amount[i]; // add the number off sieved 
					process_phase[i]=0;  // mark the process as broadcasting the information
				}
			} else if ( process_phase[i]==2 ){
				MPI_Test(availabl_reqs+i,&flag,&status);
				if (flag == true )
					process_phase[i]=0;  // mark the process as free again! 
			}
		}
	}
}
