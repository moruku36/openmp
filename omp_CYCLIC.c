#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define ITER (512*1024*1024)

main(argc, argv )
    int   argc;
    char  *argv[];
{
        int numthreads, MyID;
        int i, j, sum=0;
        double x, y;

        #pragma omp parallel
        numthreads = omp_get_num_threads();

  for( i = 0; i < ITER; i++ ){
    x = drand48();
    y = drand48();
    if(x*x+y*y < 1) sum++;
  }

  printf( "PI=%f¥n", 4.0*sum/ITER );
}