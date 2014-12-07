#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "compute.h"

#ifdef GEN_PICTURES
static void do_draw(const struct parameters *p,
		size_t key, size_t h, size_t w,
		double (*restrict g))
{
	begin_picture(key, w-2, h-2, p->io_tmin, p->io_tmax);
	size_t i, j;
	for (i = 1; i < h-1; ++i)
		for (j = 1; j < w-1; ++j)
			draw_point(j-1, i-1, g[(i) * w + j]);
	end_picture();
}
#endif

void do_compute(const struct parameters* p, struct results *r)
{
	size_t i, j;

	/* alias input parameters */
	const double (*restrict tinit) = (const double *)p->tinit;
	const double (*restrict cinit) = (const double *)p->conductivity;

	/* allocate grid data */
	const size_t h = p->N + 2;
	const size_t w = p->M + 2;
	double (*restrict g1) = (double (*restrict)) malloc(h * w * sizeof(double));
	double (*restrict g2) = (double (*restrict)) malloc(h * w * sizeof(double));

	/* allocate halo for conductivities */
	double (*restrict c) = (double (*restrict)) malloc(h * w * sizeof(double));

	struct timeval before;

	static const double c_cdir = 0.25 * M_SQRT2 / (M_SQRT2 + 1.0);
	static const double c_cdiag = 0.25 / (M_SQRT2 + 1.0);

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
	double (*restrict src) = g2;
	double (*restrict dst) = g1;

	gettimeofday(&before, NULL);

	// need to copy in both dst (which becomes src) and src (because it has the constant smearing)
	#pragma acc enter data copyin(c[0:w*h], dst[0:w*h], src[0:w*h])
	for (iter = 1; iter <= p->maxiter; ++iter)
	{
#ifdef GEN_PICTURES
		const unsigned drawrate = 10;
		#pragma acc update host(src[0:w*h])
		if (iter % drawrate == 0)
			do_draw(p, iter/drawrate, h, w, src);
#endif

		/* swap source and destination */
		{ void *tmp = src; src = dst; dst = tmp; }

		/* compute */
		#pragma acc parallel loop gang num_gangs(h-2) num_workers(4) vector_length(24) independent present(c[0:w*h], dst[0:w*h], src[0:w*h])
		for (i = 1; i < h - 1; ++i)
		{
			#pragma acc loop worker independent
			for (j = 1; j < w - 1; ++j)
			{
				double w_local = c[(i) * w + j];
				double restw = 1.0 - w_local;

				dst[(i) * w + j] = w_local * src[(i) * w + j] + 

					(src[(i+1) * w + j  ] + src[(i-1) * w + j  ] + 
					 src[(i  ) * w + j+1] + src[(i  ) * w + j-1]) * (restw * c_cdir) +

					(src[(i-1) * w + j-1] + src[(i-1) * w + j+1] + 
					 src[(i+1) * w + j-1] + src[(i+1) * w + j+1]) * (restw * c_cdiag);

			}

			/* copy left and right column to opposite border */
			dst[(i) * w + w-1] = dst[(i) * w + 1];
			dst[(i) * w + 0] = dst[(i) * w + w-2];
		}

		/* conditional reporting */
		if (iter % p->period == 0 || iter == p->maxiter) {
			r->niter = iter;

			struct timeval after;
			gettimeofday(&after, NULL);

			double tmin = DBL_MAX;
			double tmax = DBL_MIN;
			double maxdiff = DBL_MIN;
			double tavg = 0.0;

			// iterate over non-constant rows
			#pragma acc parallel loop gang num_gangs(h-2) num_workers(4) vector_length(24) reduction(min: tmin) reduction(max: tmax) reduction(max: maxdiff) reduction(+: tavg) present(dst[0:w*h], src[0:w*h])
			for (size_t y = 1; y < w - 1; ++y) {
				#pragma acc loop worker independent
				for (size_t x = 1; x < h - 1; ++x) {
					double result = dst[w*y + x];

					// update minimum temprature if needed
					if (result < tmin)
						tmin = result;
					// update maximum temprature if needed
					if (result > tmax)
						tmax = result;
					// update the sum (not yet average)
					tavg += result;
					// update the maximum temprature difference if needed
					double diff = fabs(src[w*y + x] - result);
					if (diff > maxdiff)
						maxdiff = diff;
				}
			}

			r->tmin = tmin;
			r->tmax = tmax;
			r->maxdiff = maxdiff;
			r->tavg = tavg / (p->N * p->M);

			r->time = (double)(after.tv_sec - before.tv_sec) + 
				(double)(after.tv_usec - before.tv_usec) / 1e6;

			int done = (maxdiff < p->threshold || iter == p->maxiter);

			if (p->printreports)
				report_results(p, r);

			if (done)
				break;
		}
	}
	// We don't need this, but a 'real' application would (and it's outside our timings, so harmless).
	#pragma acc exit data copyout(dst[0:w*h], src[0:w*h])

	free(c);
	free(g2);
	free(g1);
}
