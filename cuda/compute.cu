#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
extern "C" {
#include "compute.h"
}

#define CUDA_SAFE_CALL(ret) { gpuAssert((ret), __FILE__, __LINE__); }

// http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

/* Does the reduction step and return if the convergence has setteled */
static int fill_report(const struct parameters *p, struct results *r,
                        size_t h, size_t w, 
                        double * a,
                        double * b,
                        double iter,
                        struct timeval *before)
{
    /* compute min/max/avg */
    double tmin = DBL_MAX, tmax = DBL_MIN;
    double sum = 0.0;
    double maxdiff = 0.0;
    struct timeval after;

    /* We have said that the final reduction does not need to be included. */
    gettimeofday(&after, NULL);
 
    for (size_t i = 1; i < h - 1; ++i)
        for (size_t j = 1; j < w - 1; ++j) 
        {
            double v = a[(i) * w + j];
            double v_old = b[(i) * w + j];
            double diff = fabs(v - v_old);
            sum += v;
            if (tmin > v) tmin = v;
            if (tmax < v) tmax = v;
            if (diff > maxdiff) maxdiff = diff;
        }

    r->niter = iter;
    r->maxdiff = maxdiff;
    r->tmin = tmin;
    r->tmax = tmax;
    r->tavg = sum / (p->N * p->M);

    r->time = (double)(after.tv_sec - before->tv_sec) + 
        (double)(after.tv_usec - before->tv_usec) / 1e6;

    return (maxdiff >= p->threshold) ? 0 : 1;
}

static const double c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
static const double c_cdiag = 0.25 / (M_SQRT2 + 1.0);

void __global__ do_calc(double *c, double *src, double *dst, size_t w, size_t h)
{
    unsigned x = blockIdx.x * blockDim.x + threadIdx.x + 1;
    unsigned y = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if (x >= w-1)
        return;
    if (y >= h-1)
        return;

    double w_local = c[(y) * w + x];
    double restw = 1.0 - w_local;

    dst[(y) * w + x] = w_local * src[(y) * w + x] + 

	    (src[(y+1) * w + x  ] + src[(y-1) * w + x  ] + 
	     src[(y  ) * w + x+1] + src[(y  ) * w + x-1]) * (restw * c_cdir) +

	    (src[(y-1) * w + x-1] + src[(y-1) * w + x+1] + 
	     src[(y+1) * w + x-1] + src[(y+1) * w + x+1]) * (restw * c_cdiag);
}

void __global__ do_smear(double *dst, size_t w, size_t h) {
    unsigned y = blockIdx.x + 1;
    if (y >= h-1)
        return;

    dst[(y) * w + w-1] = dst[(y) * w + 1];
    dst[(y) * w + 0] = dst[(y) * w + w-2];
}

void do_compute(const struct parameters* p, struct results *r)
{
    size_t i, j;

    /* alias input parameters */
    const double *tinit = (const double *)p->tinit;
    const double *cinit = (const double *)p->conductivity;

    /* allocate grid data */
    const size_t h = p->N + 2;
    const size_t w = p->M + 2;
    double *g1 = (double *) malloc(h * w * sizeof(double));
    double *g2 = (double *) malloc(h * w * sizeof(double));

    /* allocate halo for conductivities */
    double *c = (double *) malloc(h * w * sizeof(double));

    /* set initial temperatures and conductivities */
    for (i = 1; i < h - 1; ++i)
        for (j = 1; j < w - 1; ++j) 
        {
            g1[(i) * w + j] = tinit[(i-1) * p->M + j-1];
            c[(i) * w + j] = cinit[(i-1) * p->M + j-1];
        }
    /* smear outermost columns to border */
    for (j = 1; j < w-1; ++j) {
        g1[(0) * w + j] = g2[(0) * w + j] = g1[(1) * w + j];
        g1[(h-1) * w + j] = g2[(h-1) * w + j] = g1[(h-2) * w + j];
    }
    /* smear outermost rows to borders */
    for (i = 0; i < h; ++i)
    {
        g1[(i) * w + w-1] = g2[(i) * w + w-1] = g1[(i) * w + 1];
        g1[(i) * w + 0] = g2[(i) * w + 0] = g1[(i) * w + w-2];
    }

    /* compute */
    size_t iter;
    double *src = g2;
    double *dst = g1;

    double *d_c, *d_dst, *d_src;
    CUDA_SAFE_CALL( cudaMalloc(&d_src, sizeof(double) * w * h ) );
    CUDA_SAFE_CALL( cudaMalloc(&d_dst, sizeof(double) * w * h ) );
    CUDA_SAFE_CALL( cudaMalloc(&d_c, sizeof(double) * w * h ) );
    CUDA_SAFE_CALL( cudaMemcpy(d_c, c, sizeof(double) * w * h, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy(d_dst, dst, sizeof(double) * w * h, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy(d_src, src, sizeof(double) * w * h, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaDeviceSynchronize() );
	
    struct timeval before;
    gettimeofday(&before, NULL);

    for (iter = 1; iter <= p->maxiter; ++iter)
    {
        /* swap source and destination */
        { double *tmp = src; src = dst; dst = tmp; }
        { double *tmp = d_src; d_src = d_dst; d_dst = tmp; }

        const unsigned gridSize = 16;
        dim3 blockSize((h-2+gridSize-1)/gridSize, (w-2+gridSize-1)/gridSize);
        dim3 threadSize(gridSize, gridSize);
        do_calc<<<blockSize, threadSize>>>(d_c, d_src, d_dst, w, h);
        CUDA_SAFE_CALL( cudaDeviceSynchronize() );
        do_smear<<<h-2, 1>>>(d_dst, w, h);
        CUDA_SAFE_CALL( cudaDeviceSynchronize() );

        /* conditional reporting */
        if (iter % p->period == 0) {
            CUDA_SAFE_CALL( cudaMemcpy(dst, d_dst, sizeof(double) * w * h, cudaMemcpyDeviceToHost) );
            CUDA_SAFE_CALL( cudaMemcpy(src, d_src, sizeof(double) * w * h, cudaMemcpyDeviceToHost) );
            CUDA_SAFE_CALL( cudaDeviceSynchronize() );
            if(fill_report(p, r, h, w, dst, src, iter, &before)) {iter++; break;}
            if(p->printreports) report_results(p, r);
        }
    }

    /* report at end in all cases */
    iter--;
    CUDA_SAFE_CALL( cudaMemcpy(dst, d_dst, sizeof(double) * w * h, cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaMemcpy(src, d_src, sizeof(double) * w * h, cudaMemcpyDeviceToHost) );
    CUDA_SAFE_CALL( cudaDeviceSynchronize() );

    fill_report(p, r, h, w, dst, src, iter, &before);

    free(c);
    free(g2);
    free(g1);
}
