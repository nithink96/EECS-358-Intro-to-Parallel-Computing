/*  This was modified from Lab1. So I have included only the main() and gauss() functions.

 1) Processor 0 sends rows to other processors(One to all scatter). Rows are distributed in terms of interleaved scheduling, which means processor 1 will receive 
  Row[1], Row[1+numprocs], Row[1+2*numprocs]... After that, each processor has the rows that it will use in the future.
  
  2) Through each loop, the processor that hold the 'norm' row, which means (myid == norm%numprocs), 
  will be responsible for broadcasting Row[norm] to all other processors. This processor will also compute the rows in its range.
  
  3) Once the processor recieves rthe norm row, they will begin computation. 
  4) A barrier here is kept until all processors complete their work, then next iteration will continue.*/

int main(int argc, char **argv) {
  
  double start_time = 0.0, end_time;
  int len;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  //ID = argv[argc-1];
  //argc--;
  MPI_INIT(&argc, &argv):
  MPI_COMM_SIZE(MPI_COMM_WORLD, &procs);
  MPI_COMM_RANK(MPI_COMM_WORLD, &myid);
  MPI_GET_PROCESSOR_NAME(processor_name, &len);
  /* Process program parameters */
  if(myid==0)
  {
  parameters(argc, argv);

  /* Initialize A and B */
  initialize_inputs();

  /* Print input matrices */
  print_inputs();
  start_time = MPI_WTIME();
}

  /* Gaussian Elimination */
  gauss();

  if(myid == 0)
  {
      print_X();
      end_time = MPI_WTIME();
      printf("Tome = %f\n", end_time - start_time);
  }
  MPI_FINALIZE():
  return 0;

}

/* ------------------ Above Was Provided --------------------- */

/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, procs, A[][], B[], and X[],
 * defined in the beginning of this code.  X[] is initialized to zeros.
 */

void gauss() {
  int i, norm, row, col;  /* Normalization row, and zeroing
                           * element row and col */
  float multiplier;
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

  for(norm = 0; norm < N-1; norm++) {
	 if(norm == 0){
      if(myid == 0){
        for(i = 1; i < procs; i++) {
              for(row=i;row<N;row+=procs){
            MPI_Send(&A[row],N,MPI_FLOAT,i,0,MPI_COMM_WORLD);
            MPI_Send(&B[row],1,MPI_FLOAT,i,1,MPI_COMM_WORLD);
          }
        }
      }
     else{
        for(row=myid;row<N;row+=procs){
          MPI_Status Areceiver,Breceiver;
          MPI_Recv(&A[row],N,MPI_FLOAT,0,0,MPI_COMM_WORLD,&Areceiver);      
          MPI_Recv(&B[row],1,MPI_FLOAT,0,1,MPI_COMM_WORLD,&Breceiver);
        }
      }
    }

    //calculating
    if(myid==norm%procs){
      for(i=0;i<procs;i++){
      if(i!=myid){
      MPI_Send(&A[norm],N,MPI_FLOAT,i,2,MPI_COMM_WORLD);
      MPI_Send(&B[norm],1,MPI_FLOAT,i,3,MPI_COMM_WORLD);
      }
      }
      for(row=norm+procs;row<N;row+=procs){
          multiplier=A[row][norm]/A[norm][norm];
          for(col=norm;col<N;col++){
            A[row][col]-=A[norm][col]*multiplier;
          }
          B[row]-=B[norm]*multiplier;
        }
    }
    else{
        MPI_Status ReceiveAStatus,ReceiveBStatus;
        MPI_Recv(&A[norm],N,MPI_FLOAT,norm%procs,2,MPI_COMM_WORLD,&ReceiveAStatus);
        MPI_Recv(&B[norm],1,MPI_FLOAT,norm%procs,3,MPI_COMM_WORLD,&ReceiveBStatus);
        for(row=myid;row<N;row+=procs){
          if(row>norm){
          multiplier=A[row][norm]/A[norm][norm];
          for(col=norm;col<N;col++){
            A[row][col]-=A[norm][col]*multiplier;
          }
          B[row]-=B[norm]*multiplier;
        }
      }
    }

      MPI_Barrier(MPI_COMM_WORLD);
    }




    
  /* (Diagonal elements are not normalized to 1.  This is treated in back
   * substitution.)
   */

  /* Back substitution */
  for (row = N - 1; row >= 0; row--) {
    X[row] = B[row];
    for (col = N-1; col > row; col--) {
      X[row] -= A[row][col] * X[col];
    }
    X[row] /= A[row][row];
  }
}//
