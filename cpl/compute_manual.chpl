use common;
use Time;
use BlockDist;
use AdvancedIters;

config const overlapSize = 4;

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

    var localeStarts: [0..numLocales-1] int;
    var localeWidths: [0..numLocales-1] int;
    for i in 0..numLocales-1 {
        localeWidths[i] = p.N/numLocales;
        localeStarts[i] = 1+localeWidths[0]*i;
    }
    localeWidths[numLocales-1] = p.N - localeWidths[0]*(numLocales-1);

    const TwoBigDomains = {0..p.N+1, 0..p.M+1, 0..1} dmapped Block(boundingBox={1..p.N, 1..p.M, 0..0}); // domain including both src+dest
    const BigDomain = TwoBigDomains[0..p.N+1, 0..p.M+1, 0];		   	// domain including boundary points
    const ProblemSpace = BigDomain[1..p.N, 1..p.M];
    var data: [TwoBigDomains] real;  		// source, destination temprature array
    var src = 0, dst = 1;			// indexes into array
    var c : [BigDomain] real;	     		// conductivity array

    var initsrc: [BigDomain] => data[..,..,src];
    var initdst: [BigDomain] => data[..,..,dst];

    writeln(TwoBigDomains.dist);
    writeln(ProblemSpace.dist);

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

        coforall loca in 0..numLocales-1 do
            on Locales(loca) {
                var ourdata: [{1..localeWidths[loca]+overlapSize+overlapSize, 0..p.M+1, 0..1}] real;
                if loca == 0 then {
                    // leftmost side
                    ourdata[overlapSize+1..localeWidths[loca]+overlapSize+overlapSize,..,..] = data[0..localeWidths[loca]+overlapSize,..,..];
                    // overlap from right side
                    ourdata[1..overlapSize,..,..] = data[p.M-overlapSize..p.M,..,..];
                } else if loca == numLocales - 1 then {
                    // rightmost side
                    ourdata[1..localeWidths[loca]+overlapSize,..,..] = data[p.M-localeWidths[loca]-overlapSize..p.M,..,..];
                    // overlap from left side
                    ourdata[localeWidths[loca]+overlapSize+1..localeWidths[loca]+overlapSize+overlapSize,..,..] = data[1..overlapSize,..,..];
                } else {
                    // copy as a single piece
                    ourdata[1..localeWidths[loca]+overlapSize+overlapSize,..,..] = data[localeStarts[loca]-overlapSize..localeStarts[loca]+localeWidths[loca]+overlapSize,..,..];
                }

                var osrc = src, odst = dst;
                for ouriter in 1..overlapSize {
                    var oursrc: [{1..localeWidths[loca]+overlapSize+overlapSize, 0..p.M+1}] => ourdata[..,..,osrc];
                    var ourdst: [{1..localeWidths[loca]+overlapSize+overlapSize, 0..p.M+1}] => ourdata[..,..,odst];

                    // in first round we can't recalculate at the edge, second round not one node inward, etc
                    forall ij in {1+ouriter..localeWidths[loca]+overlapSize+overlapSize-ouriter, 0..p.M+1} {
                        var cond = c(ij);
      	                var weight = 1.0 - cond;
                  	    ourdst(ij) = cond * oursrc(ij) + 
                 	        (oursrc(ij+up) + oursrc(ij+down) + oursrc(ij+left) + oursrc(ij+right)) * (weight * dir) +
                	        (oursrc(ij+upleft) + oursrc(ij+upright) + oursrc(ij+downleft) + oursrc(ij+downright)) * (weight * diag);
                    }

                    osrc <=> odst;
        	}
                osrc <=> odst;

                data[localeStarts[loca]..localeStarts[loca]+localeWidths[loca],..,src] = ourdata[overlapSize..overlapSize+localeWidths[loca],..,osrc];
                data[localeStarts[loca]..localeStarts[loca]+localeWidths[loca],..,dst] = ourdata[overlapSize..overlapSize+localeWidths[loca],..,odst];
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

