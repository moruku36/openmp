#ifdef _OPENMP
#include <omp.h>
#endif #include <stdio.h>
#define N 16

void main(){
  int i, NumT, MyID;
  
  NumT = omp_get_num_threads();
  MyID = omp_get_thread_num();
  printf("S: I am %d/%d.¥n", MyID, NumT );
  #pragma omp parallel for
  for( i = 0; i < N; i++ ){
    NumT = omp_get_num_threads();
    MyID = omp_get_thread_num();
    printf("P: i=%d, thread= %d/%d.¥n",i, MyID,NumT );
}
  NumT = omp_get_num_threads();
  MyID = omp_get_thread_num();
  printf("S: I am %d/%d.¥n", MyID, NumT );
}