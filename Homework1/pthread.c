int RowCount;
pthread_mutex_t RowCount_lock;

void *SubGauss(void * n) {
  int row, col, norm, tmpCount, chunk;
  float multiplier;

  norm = *((int *)n);
  while(RowCount < N) {
        pthread_mutex_lock(&RowCount_lock);// ensuring sequential computation to avoid incorrect values
          tmpCount = RowCount;
          chunk = (N - tmpCount + 1)/(2*procs) + 1;//calculating chunk size
          RowCount += chunk;
        pthread_mutex_unlock(&RowCount_lock);

        for (row = tmpCount + 1; row < tmpCount + chunk + 1 && row < N; row++) {
      multiplier = A[row][norm] / A[norm][norm];
      for (col = norm; col < N; col++) {
        A[row][col] -= A[norm][col] * multiplier;
      }
      B[row] -= B[norm] * multiplier;
    }
  }
}

void gauss() {
  int i, norm, row, col;  /* Normalization row, and zeroing
                           * element row and col */
  pthread_t threads[procs];

  pthread_mutex_init(&RowCount_lock,NULL);

  for(norm = 0; norm < N-1; norm++) {
        RowCount = norm;
    for(i = 0; i < procs; i++) {
      pthread_create(&threads[i], NULL, &SubGauss, (void *)(&norm));
    }
    for(i = 0; i < procs; i++) {
      pthread_join(threads[i], NULL);//suspends execution of the current thread and waits for other thread to complete
    }
  }
  pthread_mutex_destroy(&RowCount_lock);
  /* (Diagonal elements are not normalized to 1.  This is treated in back
   * substitution.)
   */

  /* Back substitution */
  for (row = N - 1; row >= 0; row--) {
    X[row] = B[row];
    for (col = N-1; col > row; col--) {
      X[row] -= A[row][col] * X[col];

