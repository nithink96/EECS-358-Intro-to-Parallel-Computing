void gauss() {
  int norm, row, col;  /* Normalization row, and zeroing
                        * element row and col */
  float multiplier;

  //printf("Computing Serially.\n");

  /* Gaussian elimination */
  for (norm = 0; norm < N - 1; norm++) {
    /* setting row, multiplier and col as local variables  to be changed during the last 
      iteration of the loop */
    #pragma omp parallel private(row, multiplier, col) num_threads(procs)
    // using dynamic scheduling
        #pragma omp for schedule(dynamic)
        for (row = norm + 1; row < N; row++) {
      multiplier = A[row][norm] / A[norm][norm];
      for (col = norm; col < N; col++) {
        A[row][col] -= A[norm][col] * multiplier;
      }
      B[row] -= B[norm] * multiplier;
    }
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
}
