# include <stdio.h>
# include <string.h>
# include <math.h>
# include "util.h"
# include <mpi.h>
# include <time.h>
# include <stdlib.h>

int main(int argc, char ** argv)
{
 
  
  

  int mpisize,mpirank,N,i,iter,MaxIter,my_n;
  double h,aij,aii,my_residue,residue,residue0,uTop,uBot,uTemp;
  double *uOld,*uNew;
  //timestamp_type time1,time2;
  double time1,time2;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  if(argc!=3){
    fprintf(stderr,"Function needs vector size as input argument!\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  N= atol(argv[1]);
  h=1./(N+1);
  aii=2/h/h;
  aij=-1/h/h;
  
  
  MaxIter= atol(argv[2]);
    
 
  if(N % mpisize !=0){
    fprintf(stderr,"Vector size not divisible!\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  // printf("rank %d/%d reporting for duty\n",mpirank,mpisize);

  my_n = N/mpisize;
  residue0=sqrt(N);
  residue=residue0;

  uOld=(double *)malloc(my_n*sizeof(double));
  uNew=(double *)malloc(my_n*sizeof(double));
  uTop=0.;
  uBot=0.;
  MPI_Status status;


  for(i=0;i<my_n;i++){
    uOld[i]=0.;
    uNew[i]=0.;
  }

  //get_timestamp(&time1);
  iter=0;
  time1=MPI_Wtime();

  //clock_t t;
  //t=clock();


  do{
    if(mpirank==0){
      uNew[0]=(1-uOld[1]*aij)/aii;
      uNew[my_n-1]=(1-uOld[my_n-2]*aij-uBot*aij)/aii;
    }
    else if(mpirank==mpisize-1){
      uNew[0]=(1-uTop*aij-uOld[1]*aij)/aii;
      uNew[my_n-1]=(1-uOld[my_n-2]*aij)/aii;
    }
    else{
      uNew[0]=(1-uTop*aij-uOld[1]*aij)/aii;
      uNew[my_n-1]=(1-uOld[my_n-2]*aij-uBot*aij)/aii;
    }

    for(i=1;i<my_n-1;i++){
      uNew[i]=(1-uOld[i-1]*aij-uOld[i+1]*aij)/aii;
    }
    for(i=0;i<my_n;i++){
      uOld[i]=uNew[i];
    }

  
  
  if(mpirank==0){
    uTemp=uNew[my_n-1];
    MPI_Send(&uTemp,1,MPI_DOUBLE,mpirank+1,99,MPI_COMM_WORLD);
    MPI_Recv(&uBot,1,MPI_DOUBLE, mpirank+1,99,MPI_COMM_WORLD,&status);

  }else if(mpirank==mpisize-1){
    uTemp=uNew[0];
    MPI_Send(&uTemp,1,MPI_DOUBLE,mpirank-1,99,MPI_COMM_WORLD);
    MPI_Recv(&uTop,1,MPI_DOUBLE, mpirank-1,99,MPI_COMM_WORLD,&status);
  }else{
    uTemp=uNew[0];    
    MPI_Send(&uTemp,1,MPI_DOUBLE,mpirank-1,99,MPI_COMM_WORLD);
    MPI_Recv(&uTop,1,MPI_DOUBLE, mpirank-1,99,MPI_COMM_WORLD,&status);
    uTemp=uNew[my_n-1];
    MPI_Send(&uTemp,1,MPI_DOUBLE,mpirank+1,99,MPI_COMM_WORLD);
    MPI_Recv(&uBot,1,MPI_DOUBLE, mpirank+1,99,MPI_COMM_WORLD,&status);
   
   
}
  
  /*
  if(mpirank==0){
    my_residue=(uNew[0]*aii+uNew[1]*aij-1)*(uNew[0]*aii+uNew[1]*aij-1);
    for(i=1;i<my_n-1;i++){
      my_residue=my_residue+(uNew[i-1]*aij+uNew[i]*aii+uNew[i+1]*aij-1)*(uNew[i-1]*aij+uNew[i]*aii+uNew[i+1]*aij-1);
} 
    my_residue=my_residue+(uNew[my_n-2]*aij+uNew[my_n-1]*aii+uBot*aij-1)*(uNew[my_n-2]*aij+uNew[my_n-1]*aii+uBot*aij-1);
    MPI_Allreduce(&my_residue, &residue,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    residue=sqrt(residue);

  }else if(mpirank==mpisize-1){
    my_residue=(uNew[my_n-2]*aij+uNew[my_n-1]*aii-1)*(uNew[my_n-2]*aij+uNew[my_n-1]*aii-1);
    for(i=1;i<my_n-1;i++){
      my_residue=my_residue+(uNew[i-1]*aij+uNew[i]*aii+uNew[i+1]*aij-1)*(uNew[i-1]*aij+uNew[i]*aii+uNew[i+1]*aij-1);
    }
    my_residue=my_residue+(uTop*aij+uNew[0]*aii+uNew[1]*aij-1)*(uTop*aij+uNew[0]*aii+uNew[1]*aij-1);
    MPI_Allreduce(&my_residue, &residue,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    residue=sqrt(residue);
  }else{

    my_residue=(uTop*aij+uNew[0]*aii+uNew[1]*aij-1)*(uTop*aij+uNew[0]*aii+uNew[1]*aij-1);
    for(i=1;i<my_n-1;i++){
      my_residue=my_residue+(uNew[i-1]*aij+uNew[i]*aii+uNew[i+1]*aij-1)*(uNew[i-1]*aij+uNew[i]*aii+uNew[i+1]*aij-1);
    }
    my_residue=my_residue+(uNew[my_n-2]*aij+uNew[my_n-1]*aii+uBot*aij-1)*(uNew[my_n-2]*aij+uNew[my_n-1]*aii+uBot*aij-1);
    MPI_Allreduce(&my_residue, &residue,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    residue=sqrt(residue);
   
}
  */
  iter++;



  /*
  if(mpirank==0){
    printf("The %ld th try. Residue is: %f\n",iter,residue);
}
  */

  }while(residue>residue0*0.000001 && iter<MaxIter);

  //get_timestamp(&time2);
  //t=clock()-t;

  time2=MPI_Wtime();


  //double elapsed = timestamp_diff_in_seconds(time1,time2);
  //double elapsed=(double) t/CLOCKS_PER_SEC;
  double elapsed=time2-time1;
  if(mpirank==0){
    printf("Time elapsed is %f seconds.\n",elapsed);
  }


  free(uNew);
  free(uOld);
 

  MPI_Finalize();
  return 0;
}


