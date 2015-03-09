/* Simple send-receive example */
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h>
#include "util.h"

int main(int argc, char *argv[])
{
  int rank,size,i=0;
  long N;
  timestamp_type time1,time2;
  N=atol(argv[1]);

  char hostname[1024];
  gethostname(hostname,1024);


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int message=0;
  MPI_Status status;

  get_timestamp(&time1);

  if(rank==0){
    MPI_Send(&message, 1, MPI_INT, rank+1, i, MPI_COMM_WORLD);
  }
  while(i<N){
   
   
    MPI_Recv(&message, 1, MPI_INT, (rank-1+size)%size, i, MPI_COMM_WORLD, &status);
    message=message+rank;
    //printf("Rank %d hosted on %s,the %d try, message received is %d\n",rank,hostname,i+1,message);

    //printf("Message from %d\n",(rank-1+size)%size);
    //printf("Message to %d\n",(rank+1)%size);
    if(rank==0){
      if(i<N-1){
	MPI_Send(&message, 1, MPI_INT, (rank+1)%size, i+1, MPI_COMM_WORLD);
      }}
      else{
      MPI_Send(&message, 1, MPI_INT, (rank+1)%size, i, MPI_COMM_WORLD);
      }
    i++;
  }
  get_timestamp(&time2);

  double elapsed=timestamp_diff_in_seconds(time1,time2);
  double singleCommunicationTime=elapsed/size/N;
  MPI_Finalize();

  if(rank==0){
    printf("Last message is %d\n",message);
    printf("Total elapsed time is %f seconds\n",elapsed);
    printf("Average time for single communication: %.12f seconds\n",singleCommunicationTime);
}

  return 0;
}
