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

/* Does the reduction step and return if the convergence has setteled */
static int fill_report(const struct parameters *p, struct results *r,
                        size_t h, size_t w, 
                        double (*restrict a),
                        double (*restrict b),
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
    	//#pragma acc parallel loop gang independent present(c[0:w*h], dst[0:w*h], src[0:w*h])
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
        if (iter % p->period == 0) {
            #pragma acc update host(dst[0:w*h], src[0:w*h])
            if(fill_report(p, r, h, w, dst, src, iter, &before)) {iter++; break;}
            if(p->printreports) report_results(p, r);
        }
    }
    #pragma acc exit data copyout(dst[0:w*h], src[0:w*h])

    /* report at end in all cases */
    iter--;
    fill_report(p, r, h, w, dst, src, iter, &before);

    free(c);
    free(g2);
    free(g1);
}
