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
    const TwoBigDomains = {0..p.N+1, 0..p.M+1, 0..1} dmapped Block(boundingBox={0..p.N+1, 0..p.M+1, 0..0}); // domain including both src+dest
    var data: [TwoBigDomains] real;  		// source, destination temprature array
    var src = 0, dst = 1;			// indexes into array
    var c : [BigDomain] real;	     		// conductivity array

    var initsrc: [BigDomain] => data[..,..,src];
    var initdst: [BigDomain] => data[..,..,dst];

    // copy initial values
    initdst[ProblemSpace] = p.tinit;
    c[ProblemSpace] = p.tcond;

    for i in 1..p.M {				//smear upper-lower bounds
	initdst[0,i] = initdst[1,i];
        initsrc[0,i] = initdst[1,i];

	initdst[p.N+1,i] = initdst[p.N,i];
	initsrc[p.N+1,i] = initdst[p.N,i];
    }
    for i in 0..p.N+1 { 				//smear wrap-around bounds
	initdst[i,0] = initdst[i,p.M];
	initsrc[i,0] = initdst[i,p.M];

	initdst[i,p.M+1] = initdst[i,1];
	initsrc[i,p.M+1] = initdst[i,1];
    }
  
    t.start(); 					//start timer
    for iteration in 1..p.maxiter {
	dst <=> src;				//swap buffers
        var cursrc: [BigDomain] => data[..,..,src];
        var curdst: [BigDomain] => data[..,..,dst];
        forall i in 1..p.N {
          for j in 1..p.M {
            var ij = (i,j);

            var cond = c(ij);
	    var weight = 1.0 - cond;
	    curdst(ij) = cond * cursrc(ij) + 
	        (cursrc(ij+up) + cursrc(ij+down) + cursrc(ij+left) + cursrc(ij+right)) * (weight * dir) +
	        (cursrc(ij+upleft) + cursrc(ij+upright) + cursrc(ij+downleft) + cursrc(ij+downright)) * (weight * diag);
          }
	}
	forall i in 1..p.N {
            curdst[i,0] = curdst[i,p.M];
            curdst[i,p.M+1] = curdst[i,1];
	}

	if (iteration % p.period == 0 || iteration == p.maxiter) {

            /* *** first variant: for loop *** */

/*	    var tmin = max(real);
	    var tmax = min(real);
	    var maxdiff = min(real);
	    var tavg = 0.0;
	    for i in 1..p.N {
            	for j in 1..p.M {  
		    var res = curdst[i,j];
		    var old = cursrc[i,j];
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

/*            r.maxdiff = max reduce [ij in ProblemSpace] abs(curdst[ij] - cursrc[ij]);
	    r.tmin = min reduce curdst[ProblemSpace];
	    r.tmax = max reduce curdst[ProblemSpace];
	    r.tavg = (+ reduce curdst[ProblemSpace]) / (p.N * p.M);*/
            
            /* *** third variant: custom Chapel reduction */

            r.maxdiff = max reduce [ij in ProblemSpace] abs(curdst[ij] - cursrc[ij]);
            var reduction = heatReduction reduce curdst[ProblemSpace];
            r.tmin = reduction.tmin;
            r.tmax = reduction.tmax;
            r.tavg = reduction.tsum / (p.N * p.M);

            /* shared code :) */

	    r.niter = iteration;
	    r.time = t.elapsed();
	    var done = (r.maxdiff < p.threshold || iteration == p.maxiter);
	    if (p.printreports) then report_results(p, r);
	}
    }

    /* Stop timer, set final timing */
    t.stop();
    r.time = t.elapsed();

    return r;
}

