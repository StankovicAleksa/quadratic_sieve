#include <cstdio>
#include <vector>
#include <mpi.h>
using namespace std;
int a;
MPI_Request req[3];
int flag;MPI_Status status;
int mpi_size,mpi_rank;
int main(int argc, char*argv[]){
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	a=mpi_rank;
	if ( mpi_rank == 0 || mpi_rank== 1){
		
		for ( int i=0;i<3;i++)
			MPI_Ibcast(&a,1,MPI_INT,i,MPI_COMM_WORLD,req+i);
	}
	else{
		for ( int i=0;i<10000000;i++);
		for ( int i=0;i<3;i++)
			MPI_Ibcast(&a,1,MPI_INT,i,MPI_COMM_WORLD,req+i);
	}
	MPI_Wait(req+1,&status);
	printf("%d\n",a);
	MPI_Finalize();	

	return 0;
}
