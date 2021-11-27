
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define ITER (512*1024*1024)

double parallel_rng(unsigned long long* seed)
{
    unsigned long long random = *seed;

    random ^= random >> 12;
    random ^= random << 25;
    random ^= random >> 27;

    *seed = random;

    return ((double) random * 0x2545F4914F6CDD1DULL) / ((double) 0xFFFFFFFFFFFFFFFFULL);
}

int main(int argc, char* argv[])
{
    int i;
    int sum = 0;
    int partial_sum;
    unsigned long long seed;

#pragma omp parallel private(partial_sum, seed)
    {
        partial_sum = 0;

#pragma omp critical
        {
            seed = (unsigned long long) rand();
            seed <<= 16;
            seed += (unsigned long long) rand();
            seed <<= 16;
            seed += (unsigned long long) rand();
            seed <<= 16;
            seed += (unsigned long long) rand();
        }

#pragma omp for
        for(i = 0; i < ITER; ++i)
        {
            double x = parallel_rng(&seed);
            double y = parallel_rng(&seed);

            if(x * x + y * y < 1.0)
            {
                ++partial_sum;
            }
        }

#pragma omp critical
        {
            sum += partial_sum;
        }
    }

    printf( "PI=%f¥n", 4.0*sum/ITER );
}
