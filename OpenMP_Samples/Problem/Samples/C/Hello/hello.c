#include <stdio.h>
#include "omp.h"

int  main(int argc, char* argv[]) {

     int    i;
      
#pragma omp parallel for
     for (i=0; i<10; i++) {
       printf("Hello parallel world!  i:%d \n", i);
     }

}


