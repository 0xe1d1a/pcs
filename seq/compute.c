#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> // memcpy
#include <float.h> // DBL
#include <errno.h>
#include "fail.h"
#include "compute.h"

#define TRUE 1

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
	double dirnmul = sqrt(2)/(sqrt(2) + 1);
	double diagnmul = 1.0/(sqrt(2) + 1);

	const double *cond = p->conductivity;

	// allocate two temperature matrices, copy the initial state into t_prev
	double *t_prev, *t_next;
	t_prev = malloc(sizeof(double) * p->M * p->N);
	if (t_prev == 0) die("Out of memory");
	t_next = malloc(sizeof(double) * p->M * p->N);
	if (t_next == 0) die("Out of memory");


	memcpy(t_prev, p->tinit, sizeof(double) * p->M * p->N); //@Alyssa i'd like to discuss memory efficiency as posed to the forum

	// iteration count
	size_t iter = 0;
	
	// start timing
	struct timeval tv_start, tv_curr, tv_diff;
	int ret = gettimeofday(&tv_start, NULL);
	if (ret == -1) die(strerror(errno));
	
	r->maxdiff = DBL_MAX; // something comfortably over the threshold

	while (TRUE) { 
		int done = 0;
		int do_reduction = 0;
		if (iter++ >= p->maxiter - 1) done = 1; // do the final iteration and finish
		do_reduction = done; // if we are done do a final reduction
		if (iter % p->period == 0) do_reduction = 1; // do a midrun reduction 

		if (do_reduction) {
			r->niter = iter;
			r->tmin = DBL_MAX;
			r->tmax = DBL_MIN;
			r->maxdiff = DBL_MIN;
			r->tavg = 0.0; // we don't care about the fp inaccuracy, right?
		}


		// N rows, M columns
		double *dst = t_next;
		// traverse rows
		for (size_t y = 0; y < p->N; ++y) {
			const double *cnd = &cond[p->M * y]; 

			const double *row = &t_prev[p->M * y];
			// boundary conditions (from the initial state)
			const double *rowup = &t_prev[p->M * (y - 1)];
			if (y == 0)
				rowup = &p->tinit[0];
			const double *rowdown = &t_prev[p->M * (y + 1)];
			if (y == p->N-1)
				rowdown = &p->tinit[p->M * (p->N-1)];
			// traverse the row (columns)
			for (size_t x = 0; x < p->M; ++x, dst++) {
				// wraparound idiom
				size_t prev = (x == 0 ? p->M-1 : x-1);
				size_t next = (x == p->M-1 ? 0 : x+1);

				double ourcnd = cnd[x]; // the point's conductivity
				double directcnd = (1.0 - ourcnd) * dirnmul; // direct neighbours weighted conductivity
				double diagcnd = (1.0 - ourcnd) * diagnmul; // diagonal neighbours weighted conductivity

				/* 
				*  computation of the resulting temprature based on own conductivity 
				*  plus average weighted of direct neighbours plus average weighted of
				*  diagonal neighbours. If this comment made you tired here take a break:
				*  http://imgs.xkcd.com/comics/ballmer_peak.png
				*/
				double result = row[x] * ourcnd;

				double directsum = rowup[x] + rowdown[x] + row[prev] + row[next];
				result += directcnd * directsum/4.0;
				double diagsum = rowup[prev] + rowup[next] + rowdown[prev] + rowdown[next];
				result += diagcnd * diagsum/4.0;
				/*
				*  reflect the point temprature result on the t_next state;
				*  this allows independent computations between the points
				*/
				*dst = result; 

				if (do_reduction) {
					// update minimum temprature if needed
					if (result < r->tmin)
						r->tmin = result;
					// update maximum temprature if needed
					if (result > r->tmax)
						r->tmax = result;
					// update the sum (not yet average)
					r->tavg += result;
					// update the maximum temprature difference if needed
					double diff = abs(row[x] - *dst);
					if (diff > r->maxdiff)
						r->maxdiff = diff;
				}
			}
		}
		if (do_reduction) {
			// get current time difference for the result report
			ret = gettimeofday(&tv_curr, NULL);
			if (ret == -1) die(strerror(errno));
			timersub(&tv_curr, &tv_start, &tv_diff);
			r->time = tv_diff.tv_sec + (tv_diff.tv_usec / 1000000.0);
			// compute average temprature for the result report 
			r->tavg /= (p->N * p->M);
			if (r->maxdiff < p->threshold) done = 1;

			// report state for debugging reasons and create current state frame
			if(p->printreports) {
				begin_picture(iter, p->M, p->N, p->io_tmin, p->io_tmax);
				double *tptr = t_prev;
				for (size_t y = 0; y < p->N; ++y) {
					for (size_t x = 0; x < p->M; ++x) {
						draw_point(x, y, *tptr++);
					}
				}
				end_picture();
				report_results(p, r);
			}
		}

		if (!done) {
			// swap the buffers 
			double *t_temp;
			t_temp = t_prev;
			t_prev = t_next;
			t_next = t_temp;
		}
		else // done; reached maxiter or maxdiff is smaller than threshold
			break;
	}

	free(t_prev);
	free(t_next);
}
