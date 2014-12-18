use common;
use Time;
use BlockDist;
use AdvancedIters;

record heatReductionResults
{
    type eltType;
    var tmin, tmax, tsum: eltType;
}

/* based on reference: modules/internal/ChapelReduce.chpl
   also helpful: http://faculty.knox.edu/dbunde/teaching/chapel/#Custom%20Reductions (thanks Franz) */
class heatReduction : ReduceScanOp
{
    type eltType;
    var tmin: eltType = max(eltType);
    var tmax: eltType = min(eltType);
    var tsum: chpl__sumType(eltType);

    proc accumulate(val: eltType)
    {
        if (val < tmin) then tmin = val;
        if (val > tmax) then tmax = val;
        tsum += val;
    }

    proc combine(other: heatReduction)
    {
        if (other.tmin < tmin) then tmin = other.tmin;
        if (other.tmax > tmax) then tmax = other.tmax;
        tsum += other.tsum;
    }

    proc generate()
    {
        return new heatReductionResults(eltType, tmin, tmax, tsum);
    }
}

proc do_compute(p : params) 
{
    /* Results structure */
    var r : results;
    /* Timer */
    var t : Timer;

    var iteration = 0;                    	// iteration counter
    var delta: real;                  		// measure of convergence 
    var epsilon = p.threshold;			// set threshold
    const dir = 0.25 * sqrt(2) / (sqrt(2) + 1.0); 
    const diag = 0.25 / (sqrt(2) + 1.0);
    const down = (-1,0), up = (1,0), right = (0,1), left = (0,-1);
    const upright = (1,1), downright = (-1,1), upleft = (1,-1), downleft = (-1,-1);

    const ProblemSpace = {1..p.N, 1..p.M};       // domain for grid points
    const BigDomain = {0..p.N+1, 0..p.M+1};   	// domain including boundary points
    var src, dst: [BigDomain] real;  		// source, destination temprature arrays: 
    var c : [BigDomain] real;	     		// conductivity array

    // copy initial values
    dst[ProblemSpace] = p.tinit;
    c[ProblemSpace] = p.tcond;

    for i in 1..p.M {				//smear upper-lower bounds
	dst[0,i] = dst[1,i];
        src[0,i] = dst[1,i];

	dst[p.N+1,i] = dst[p.N,i];
	src[p.N+1,i] = dst[p.N,i];
    }
    for i in 0..p.N+1 { 				//smear wrap-around bounds
	dst[i,0] = dst[i,p.M];
	src[i,0] = dst[i,p.M];

	dst[i,p.M+1] = dst[i,1];
	src[i,p.M+1] = dst[i,1];
    }
  
    t.start(); 					//start timer
    for iteration in 1..p.maxiter {
	dst <=> src;				//swap buffers
        forall i in 1..p.N {
          for j in 1..p.M {
            var ij = (i,j);

            var cond = c(ij);
	    var weight = 1.0 - cond;
	    dst(ij) = cond * src(ij) + 
	        (src(ij+up) + src(ij+down) + src(ij+left) + src(ij+right)) * (weight * dir) +
	        (src(ij+upleft) + src(ij+upright) + src(ij+downleft) + src(ij+downright)) * (weight * diag);
          }
	}
	forall i in 1..p.N {
            dst[i,0] = dst[i,p.M];
            dst[i,p.M+1] = dst[i,1];
	}

	if (iteration % p.period == 0 || iteration == p.maxiter) {

            /* *** first variant: for loop *** */

/*	    var tmin = max(real);
	    var tmax = min(real);
	    var maxdiff = min(real);
	    var tavg = 0.0;
	    for i in 1..p.N {
            	for j in 1..p.M {  
		    var res = dst[i,j];
		    var old = src[i,j];
		    if (res < tmin) then tmin = res;
		    if (res > tmax) then tmax = res;
		    tavg += res;
		    var diff = fabs(old - res);
		    if (diff > maxdiff) then maxdiff = diff;
		}
	    }
	    r.tmin = tmin;
	    r.tmax = tmax;
	    r.tavg = tavg / (p.N * p.M);
	    r.maxdiff = maxdiff;*/

            /* *** second variant: Chapel reduction */

/*            r.maxdiff = max reduce [ij in ProblemSpace] abs(dst[ij] - src[ij]);
	    r.tmin = min reduce dst[ProblemSpace];
	    r.tmax = max reduce dst[ProblemSpace];
	    r.tavg = (+ reduce dst[ProblemSpace]) / (p.N * p.M);*/
            
            /* *** third variant: custom Chapel reduction */

            r.maxdiff = max reduce [ij in ProblemSpace] abs(dst[ij] - src[ij]);
            var reduction = heatReduction reduce dst[ProblemSpace];
            r.tmin = reduction.tmin;
            r.tmax = reduction.tmax;
            r.tavg = reduction.tsum / (p.N * p.M);

            /* shared code :) */

	    r.niter = iteration;
	    r.time = t.elapsed();
	    var done = (r.maxdiff < p.threshold || iteration == p.maxiter);
	    if (0) then report_results(p, r);
	}
    }

    /* Stop timer, set final timing */
    t.stop();
    r.time = t.elapsed();

    return r;
}

