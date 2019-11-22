#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define NTHREADS 12
#define tolerancia 0.0000000000001

int main(){
omp_set_num_threads(NTHREADS);
	
int N, i;
double *random, a, b;
FILE *randArray;
N = 6000000;
random  = (double*)malloc(N*sizeof(double));
srand(1234);
randArray = fopen("randomarray.txt","w");
for(i=0;i<N;i++){
	a = rand()*(1.0/RAND_MAX);
	b = rand()*(1.0/RAND_MAX);
	random[i] = pow(-2.0*log(a + tolerancia) + 3.0*tolerancia,0.5)*cos(2.0*M_PI*b);
	fprintf(randArray, "%.16f\n", random[i]);
}
fclose(randArray);
}
