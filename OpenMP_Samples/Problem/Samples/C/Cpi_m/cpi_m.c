/*
c**********************************************************************
c   pi.c - compute pi by integrating f(x) = 4/(1 + x**2)     
c     
c   Each node: 
c    1) receives the number of rectangles used in the approximation.
c    2) calculates the areas of it's rectangles.
c    3) Synchronizes for a global summation.
c   Node 0 prints the result.
c
c  Variables:
c
c    pi  the calculated result
c    n   number of points of integration.  
c    x           midpoint of each rectangle's interval
c    f           function to integrate
c    sum,pi      area of rectangles
c    tmp         temporary scratch space for global summation
c    i           do loop index
c****************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


int main(int argc, char* argv[]) {

     double PI25DT = 3.141592653589793238462643;
     double pi, h, sum, x;
     double t1, t2;
     long   n, i;
       

        printf("Enter the number of intervals \n");
        scanf("%d",&n);
        printf("n=%d \n",n);
     
/*                                 check for quit signal */

     if ( n <= 0 ) {
        exit(1);
     }


     t1 = omp_get_wtime();


/*
                                  calculate the interval size */
      h = 1.0 / (double)n;

      sum  = 0.0;


 
#pragma omp parallel for private(x) reduction(+:sum)
      for (i = 1; i<=n; i++) {
	 x = h * ((double)i - 0.5);
         sum = sum + 4.0 / (1.0 + x*x);
      }

      pi = h * sum;

/*                                 prints the answer. */
      printf("  pi is approximately: %18.16lf Error is: %18.16lf \n", 
                 pi, fabs(pi-PI25DT));

      t2 = omp_get_wtime();  

      printf(" time = %lf [sec.] \n",t2-t1);



}

