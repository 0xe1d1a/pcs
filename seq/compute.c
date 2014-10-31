#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> // memcpy
#include <float.h> // DBL
#include "compute.h"

void do_compute(const struct parameters* p, struct results *r) {
	double dirnmul = sqrt(2)/(sqrt(2) + 1);
	double diagnmul = 1.0/(sqrt(2) + 1);

	const double *cond = p->conductivity;

	// allocate two temperature matrices, copy the initial state into t_prev
	double *t_prev, *t_next;
	t_prev = malloc(sizeof(double) * p->M * p->N);
	t_next = malloc(sizeof(double) * p->M * p->N);
	memcpy(t_prev, p->tinit, sizeof(double) * p->M * p->N);

	// iteration count
	size_t iter = 0;
	
	
	
	r->maxdiff = DBL_MAX; // something comfortably over the threshold
	//while (r->maxdiff >= p->threshold && iter++ < p->maxiter) {
	while (1) {
		int done = 0;
		int do_reduction = 0;
		if (iter++ >= p->maxiter) done = 1;
		do_reduction = done;
		if (iter % p->period == 0) do_reduction = 1;

		if (do_reduction) {
			// FIXME: update time
			r->niter = iter;
			r->tmin = DBL_MAX;
			r->tmax = DBL_MIN;
			r->maxdiff = DBL_MIN;
			r->tavg = 0.0; // we don't care about the fp inaccuracy, right?
		}


		// N rows, M columns
		double *dst = t_next;
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

			for (size_t x = 0; x < p->M; ++x, dst++) {
				// wraparound
				size_t prev = (x == 0 ? p->M-1 : x-1);
				size_t next = (x == p->M-1 ? 0 : x+1);

				double ourcnd = cnd[x];
				double directcnd = (1.0 - ourcnd) * dirnmul;
				double diagcnd = (1.0 - ourcnd) * diagnmul;

				double result = row[x] * ourcnd;

				double directsum = rowup[x] + rowdown[x] + row[prev] + row[next];
				result += directcnd * directsum/4.0;
				double diagsum = rowup[prev] + rowup[next] + rowdown[prev] + rowdown[next];
				result += diagcnd * diagsum/4.0;

				*dst = result;

				// FIXME: don't do this if we don't need it
				if (result < r->tmin)
					r->tmin = result;
				if (result > r->tmax)
					r->tmax = result;
				r->tavg += result;

				double diff = abs(row[x] - *dst);
				if (diff > r->maxdiff)
					r->maxdiff = diff;
			}
		}
		if (do_reduction) {
			r->tavg /= (p->N * p->M);
			if (r->maxdiff < p->threshold) done = 1;
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

		if (done) break;
		// swap the buffers
		double *t_temp;
		t_temp = t_prev;
		t_prev = t_next;
		t_next = t_temp;
	}

	report_results(p, r);

	free(t_prev);
	free(t_next);
}
