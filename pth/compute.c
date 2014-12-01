#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h> // printf
#include <string.h> // memcpy
#include <float.h> // DBL
#include <errno.h>
#include <assert.h>
#include "fail.h"
#include "compute.h"
#include "pthread.h"

#define TRUE 1

/* 
 * double dirnmul = sqrt(2)/(sqrt(2) + 1) / 4.0;
 * double diagnmul = 1.0/(sqrt(2) + 1) / 4.0;
 *
 * (going massively overboard on the precision just to be sure)
 */
const double dirnmul = 0.14644660940672626914249576657311990857124328613281;
const double diagnmul = 0.10355339059327377249086765687025035731494426727295;

inline void do_calc(const double * restrict cnd, double * restrict dst, size_t x, const double * restrict row, size_t width) {
	double ourcnd = cnd[x]; // the point's conductivity
	double leftover = 1.0 - ourcnd;
	double directcnd = leftover * dirnmul; // direct neighbours weighted conductivity
	double diagcnd = leftover * diagnmul; // diagonal neighbours weighted conductivity

	/* 
	*  computation of the resulting temprature based on own conductivity 
	*  plus average weighted of direct neighbours plus average weighted of
	*  diagonal neighbours. If this comment made you tired here take a break:
	*  http://imgs.xkcd.com/comics/ballmer_peak.png
	*/
	double result = row[x] * ourcnd;

	double directsum = row[x-width] + row[x+width] + row[x-1] + row[x+1];
	result += directcnd * directsum;
	double diagsum = row[x-1-width] + row[x+1-width] + row[x-1+width] + row[x+1+width];
	result += diagcnd * diagsum;
	/*
	*  reflect the point temprature result on the t_next state;
	*  this allows independent computations between the points
	*/
	dst[x] = result; 
}

pthread_barrier_t barrier;
const struct parameters *g_p;
size_t g_width;
volatile int g_do_reduction;

struct thread_params {
	size_t srow;
	size_t erow;
	size_t iter;
	// state
	double *t_prev, *t_next;

	// reduction output
	double tmin, tmax, maxdiff, tavg;
};

void thread_work(const struct thread_params *tparams)
{
	const struct parameters *p = g_p;
	const double * restrict cond = p->conductivity;
	double *t_prev = tparams->t_prev;
	double *t_next = tparams->t_next;
	const size_t width = g_width;

	// N rows, M columns
	// iterate over non-constant rows
	for (size_t y = tparams->srow; y < tparams->erow; ++y) {

		double * restrict dst = &t_next[width*y];
		const double * restrict cnd_dst = &cond[(y-1) * p->M] - 1;
		//double * volatile dst = &t_next[width*y];
		//const double * volatile cnd_dst = &cond[(y-1) * p->M] - 1;

		const double * restrict row = &t_prev[width*y];

		// traverse the non-smeared columns in the row
		for (size_t x = 1; x < p->M+1; ++x) {
			do_calc(cnd_dst, dst, x, row, width);
		}

		// smear into first/last columns
		dst[0] = dst[p->M];
		dst[p->M+1] = dst[1];
	}
}

void *thread_reduction_work(struct thread_params *tparams)
{
	const struct parameters *p = g_p;
	double *t_prev = tparams->t_prev;
	double *t_next = tparams->t_next;
	const size_t width = g_width;

	double tmin = DBL_MAX;
	double tmax = DBL_MIN;
	double maxdiff = DBL_MIN;
	double tavg = 0.0;

	for (size_t y = tparams->srow; y < tparams->erow; ++y) {
		const double * restrict row = &t_prev[width * y];
		const double * restrict dst = &t_next[width * y];

		for (size_t x = 1; x < p->M+1; ++x) {
			double result = dst[x];

			// update minimum temprature if needed
			if (result < tmin)
				tmin = result;
			// update maximum temprature if needed
			if (result > tmax)
				tmax = result;
			// update the sum (not yet average)
			tavg += result;
			// update the maximum temprature difference if needed
			double diff = fabs(row[x] - result);
			if (diff > maxdiff)
				maxdiff = diff;
		}
	}

	tparams->tmin = tmin;
	tparams->tmax = tmax;
	tparams->maxdiff = maxdiff;
	tparams->tavg = tavg;
}

void *thread_main(struct thread_params *tparams)
{
	pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
	const struct parameters *p = g_p;
	while (1)
	{
		thread_work(tparams);
		pthread_barrier_wait(&barrier); // wait for computation
		tparams->iter++;
		if (tparams->iter % p->period == 0 || tparams->iter >= p->maxiter - 1) { //do reduction
			thread_reduction_work(tparams);
			//pthread_barrier_wait(&barrier);
		}

		// swap our local copy
		double *tmp = tparams->t_prev;
		tparams->t_prev = tparams->t_next;
		tparams->t_next = tmp;
	}
}

/*
 * Function:  do_compute 
 * --------------------
 * computes the heat dissipation of a cylindrical 
 * surface presented as a N x M matrix:   
 *
 *  p: The struct containing the parameters of the input
 *  
 *  r: The struct to store the interesting results
 *       
 */
void do_compute(const struct parameters* p, struct results *r) {
	// require the input is 'nice'
	assert((p->N % 2 == 0) && (p->M % 2 == 0));
	//omp_set_num_threads(p->nthreads);
	// start timing
	struct timeval tv_start, tv_curr, tv_diff;
	int ret = gettimeofday(&tv_start, NULL);
	if (ret == -1) die(strerror(errno));
	
	g_p = p;

	//const size_t unaligned_width = p->M + 2; // including smearing columns at start/end
	//const size_t width = unaligned_width + 2 - (unaligned_width % 2); // aligned to 16 bytes (two doubles)
	const size_t width = p->M + 2;
	g_width = width;

	// allocate two temperature matrices, with extra rows on both sides and at top/bottom (i.e. (N+2)x(M+2))
	double *t_prev, *t_next;
	t_prev = malloc(sizeof(double) * width * (p->N + 2));
	if (t_prev == 0) die("Out of memory");
	t_next = malloc(sizeof(double) * width * (p->N + 2));
	if (t_next == 0) die("Out of memory");

	// copy the initial state into t_prev
	double *memtgt = t_prev;
	memcpy(memtgt + 1, p->tinit, sizeof(double) * p->M); // dup first row
	for (size_t y = 0; y < p->N; y++) {
		memtgt = t_prev + (width * (y+1));
		memcpy(memtgt + 1, p->tinit + (p->M * y), sizeof(double) * p->M); // dup normal row
	}
	memtgt = t_prev + (width * (p->N+1));
	memcpy(memtgt + 1, p->tinit + (p->M * (p->N - 1)), sizeof(double) * p->M); // dup last row

	// smear *all* rows (including initial ones)
	for (size_t y = 0; y < p->N+2; y++) {
		t_prev[width * y] = t_prev[(width * y) + p->M];
		t_prev[(width * y) + p->M + 1] = t_prev[(width * y) + 1];
	}

	// copy constant rows into the second copy
	memcpy(t_next, t_prev, sizeof(double) * width);
	memcpy(t_next + (width * (p->N+1)), t_prev + (width * (p->N+1)), sizeof(double) * width);

	pthread_barrier_init(&barrier, NULL, p->nthreads);

	// iteration count
	size_t iter = 0;
	pthread_t threads[p->nthreads];
	struct thread_params params[p->nthreads];
	size_t itersize = (p->N+1)/p->nthreads;

	for (int i=0; i<p->nthreads; i++)
	{
		params[i].iter = 0;
		params[i].srow = 1 + itersize*i;
		params[i].erow = 1 + itersize*(i+1);
		if (i == p->nthreads-1)
			params[i].erow = p->N+1;
		printf("thread %d from %d to %d\n", i, params[i].srow, params[i].erow);
		params[i].t_prev = t_prev;
		params[i].t_next = t_next;
		if(i!=0) pthread_create(&threads[i], NULL, &thread_main, &params[i]);
	}
		
	/* master thread is thread 0 */
	struct thread_params *tparams = &params[0];

	while (TRUE) { 
		int done = 0;
		g_do_reduction = 0;
		if (iter++ >= p->maxiter - 1) done = 1; // do the final iteration and finish
		g_do_reduction = done; // if we are done do a final reduction
		if (iter % p->period == 0) g_do_reduction = 1; // do a midrun reduction 

		thread_work(tparams);

		if(g_do_reduction) thread_reduction_work(tparams);

		// wait for this round of computation to be done
		pthread_barrier_wait(&barrier);

		if (g_do_reduction) {
			// reset results
			r->niter = iter;

			r->tmin = DBL_MAX;
			r->tmax = DBL_MIN;
			r->maxdiff = DBL_MIN;
			r->tavg = 0.0;


			for (int i=0; i<p->nthreads; i++) {
				if (params[i].tmin < r->tmin)
					r->tmin = params[i].tmin;
				if (params[i].tmax > r->tmax)
					r->tmax = params[i].tmax;
				if (params[i].maxdiff > r->maxdiff)
					r->maxdiff = params[i].maxdiff;
				r->tavg += params[i].tavg;
			}
			r->tavg = r->tavg / (p->N * p->M);
		}

		if (g_do_reduction) {
			// get current time difference for the result report
			ret = gettimeofday(&tv_curr, NULL);
			if (ret == -1) die(strerror(errno));
			timersub(&tv_curr, &tv_start, &tv_diff);
			r->time = tv_diff.tv_sec + (tv_diff.tv_usec / 1000000.0);

			// compute average temprature for the result report 
			if (r->maxdiff < p->threshold) done = 1;

			// report state for debugging reasons and create current state frame
			if(p->printreports) {
/*				begin_picture(iter, p->M, p->N, p->io_tmin, p->io_tmax);
				double *tptr = t_prev;
				for (size_t y = 0; y < p->N; ++y) {
					for (size_t x = 0; x < p->M; ++x) {
						draw_point(x, y, *tptr++);
					}
				}
				end_picture();*/
				report_results(p, r);
			}
		}

		if (done) // done; reached maxiter or maxdiff is smaller than threshold
			break;


		double *tmp = tparams->t_prev;
		tparams->t_prev = tparams->t_next;
		tparams->t_next = tmp;
	}

	for (int i=1; i<p->nthreads; i++)
	{
		pthread_cancel(threads[i]);
		pthread_join(threads[i], NULL);
	}

	printf("Execution time: %f\n", r->time);

	// report stuff
	FILE *fp;
	fp = fopen("results.log", "a+");
	fprintf(fp, "%f %.6e %zu %zu %zu\n", r->time, 
		(double)p->N * (double)p->M * (double)(r->niter * 12 + (double)r->niter / p->period) / r->time
		,p->N, p->M, p->nthreads);

	free(t_prev);
	free(t_next);
}

