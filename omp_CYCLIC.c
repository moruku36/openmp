#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#define N 16

void main(){
  int i, NumT, MyID[N];

  #pragma omp parallel
  NumT = omp_get_num_threads();
  printf( "Start omp parallel for¥n" );
  #pragma omp parallel for schedule(static)
  for( i = 0; i < N; i++ ){
    MyID[i] = omp_get_thread_num();
  }
  for( i = 0; i < N; i++ ){
    printf( "P: i=%2d, thread= %d/%d.\n",i, MyID[i], NumT );
  }
  printf( "Close parallel section\n" );
}