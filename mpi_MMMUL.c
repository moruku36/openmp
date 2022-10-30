#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 2048


void init_ary(double m[SIZE][SIZE]){
    int i, j;

    for( i = 0; i < SIZE; i++ ){
        for( j = 0; j < SIZE; j++ ){
            m[i][j] = i*SIZE+j+1;
        }
    }
    m[0][0] = -1;
}



void mul_ary( m, x, y )         /*  M = X * Y  */
        double m[SIZE][SIZE], x[SIZE][SIZE], y[SIZE][SIZE];
{
    int i, j, k;

    for( i = 0; i < SIZE; i++ ){
        for( j = 0; j < SIZE; j++ ){
            m[i][j] = 0.0;
            for( k = 0; k < SIZE; k++ ){
                m[i][j] += x[i][k] * y[k][j];
            }
        }
    }
}



int main( int argc, char *argv[] ){
    int    MyID = 0, NumP;
    int    i=5,j=6;
    double (*a)[SIZE], (*x)[SIZE], (*y)[SIZE];
    double proctime, sum=0.0;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &NumP );
    MPI_Comm_rank( MPI_COMM_WORLD, &MyID );

    if( MyID == 0 ){
        a = malloc( (long)sizeof(double) * SIZE * SIZE );
        x = malloc( (long)sizeof(double) * SIZE * SIZE );
        y = malloc( (long)sizeof(double) * SIZE * SIZE );

        if( (long)a+(long)x+(long)y == 0 )
            printf( "malloc error\n" );

        init_ary( x );
        init_ary( y );

        proctime = MPI_Wtime();
        mul_ary( a, x, y );
        proctime = MPI_Wtime() - proctime;

        for( i = 0; i < SIZE; i++ )
            for( j = 0; j < SIZE; j++ )
                sum += a[i][j];

        printf( "sum = %f\n", sum );
        printf( "x[%d][%d](=%f) * y[%d][%d](=%f) = a[%d][%d](=%f)\n",
                i,j,x[i][j], i,j,y[i][j], i,j,a[i][j] );
        printf( "Time = %f Sec.\n", proctime );
    }
    MPI_Finalize();
    return(0);
}
