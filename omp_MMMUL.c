#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#define N 1024

int main(){
  int i, j, k, NumT, MyID;
  double t, x[N][N], y[N][N], z[N][N];

  t=omp_get_wtime();
  #pragma omp parallel for
  for( i = 0; i < N; i++ ){
    for( j = 0; j < N; j++ ){
      for( k = 0; k < N; k++ ){
        x[i][j] = y[i] [k] + z[k][j];
      }
   }
 }
 i = lrand48() % (N); j = lrand48() % (N);
 t=omp_get_wtime() - t;
 printf( "t = %f sec. x[%d][%d]=%f¥n",t, i, j, x[i][j] );
}