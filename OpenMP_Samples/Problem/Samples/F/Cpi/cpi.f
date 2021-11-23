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

       program main
       use omp_lib

       double precision  PI25DT 
       parameter (PI25DT= 3.141592653589793238462643)
       double precision  pi, h, sum, x
       integer n, i

         print *, "Enter the number of intervals"
         read (*,*) n
         print *, "n=",n

c      === check for quit signal 
       if ( n .le. 0 )  then
         stop
       endif

c      === calculate the interval size 
       h = 1.0 / dble(n)

       sum  = 0.0;

!$omp parallel do private(x) reduction(+:sum)
       do i=1, n
         x = h * (dble(i) - 0.5)
         sum = sum + 4.0 / (1.0 + x*x)
       enddo
!$omp end parallel do

       pi = h * sum

c      ===  prints the answer
         print *, "  pi is approximately:", pi, 
     &         "Error is:", dabs(pi-PI25DT)

       stop
       end




      



