#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define N 2048

void init_ary(double m[N][N]){
    long i, j;

    for( i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            m[i][j] = i*N+j+1;
        }
    }
    m[0][0] = -1;
}

void mul_ary(m, x, y )
double m[N][N], x[N][N], y[N][N];
{
  int i, j, k,
  double t,

  t=omp_get_wtime();
  #pragma omp parallel for private(i, j, k) shared(x, y, m)
  for( i = 0; i < N; i++ ){
    for( j = 0; j < N; j++ ){
        m[i][j] =0.0;
      for( k = 0; k < N; k++ ){
        x[i][j] = y[i] [k] + z[k][j];
      }
   }
 }
/*i = lrand48() % (N); j = lrand48() % (N); */
 t=omp_get_wtime()-t;
 printf( "t = %f sec¥n",t );
}
