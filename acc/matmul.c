
#include <stdio.h>
#include <stdlib.h>

#ifdef MATMUL_ACC
#include <openacc.h>
#endif

#ifdef MATMUL_OMP
#include <omp.h>
#endif


#ifdef MATMUL_SEQ
void matmul(float *a, float *b, float *x, int len){

}
#endif

#ifdef MATMUL_ACC
void matmul(float *a, float *b, float *x, int len){

}
#endif

#ifdef MATMUL_OMP
void matmul(float *a, float *b, float *x, int len){

}
#endif

void info(){
#ifdef MATMUL_SEQ
  printf("Running matmul sequential\n");
#endif

#ifdef MATMUL_SSE
  printf("Running matmul using vectorization (sse)\n");
#endif

#ifdef MATMUL_ACC
  printf("Running matmul using OpenACC\n");
#endif

#ifdef MATMUL_OMP
  #pragma omp parallel
  printf("Running matmul using OpenMP (Thread %d of %d threads)\n", 
          omp_get_thread_num(),
          omp_get_num_threads());
#endif
}

void do_compute(int n){
  float *a = NULL, *b = NULL, *x = NULL;

  matmul(a, b, x, n);
}

int main(int argc, char **argv){

  int n = 100;

  if(argc > 1){
    n = atoi(argv[1]);
  }

#ifdef MATMUL_OMP
  int t = 1;

  if(argc > 2){
    t = atoi(argv[2]);
  }
  omp_set_num_threads(t);
#endif

  info();

  printf("Multiply %dx%d matrices.\n", n, n);

  do_compute(n);

  return 0;

}

