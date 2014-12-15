use common;
use Time;



proc do_compute(p : params) 
{
    /* Results structure */
    var r : results;
    /* Timer */
    var t : Timer;
    var iteration = 0,                    	// iteration counter
    delta: real;                      		// measure of convergence 
    var epsilon = p.threshold;			// set threshold
    const dir = 0.25 * sqrt(2) / (sqrt(2) + 1.0); 
    const diag = 0.25 / (sqrt(2) + 1.0);
    const down = (-1,0), up = (1,0), right = (0,1), left = (0,-1);
    const upright = (1,1), downright = (-1,1), upleft = (1,-1), downleft = (-1,-1);

    const ProblemSpace = {1..p.N, 1..p.M},    	// domain for grid points
          BigDomain = {0..p.N+1, 0..p.M+1};   	// domain including boundary points
    var src, dst: [BigDomain] real;  		// source, destination temprature arrays: 
    var c : [BigDomain] real;	     		// conductivity array
    //copy init values
    for i in 1..p.N {
	for j in 1..p.M {
	    dst[i,j] = p.tinit[i,j];
	    c[i,j] = p.tcond[i,j];
	}
    }
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
  
   writeln(c);
    t.start(); 					//start timer
    for iteration in 1..p.maxiter {
	dst <=> src;				//swap buffers
	writeln("Iteration ----------- ", iteration);
	writeln("dst");
	writeln(dst);
	writeln();
	writeln();
	writeln("src");
	writeln(src);
	writeln();
	for i in 1..p.N {
	    for j in 1..p.M {			//start compute
		    var ij = (i,j);
		    var cond = c(ij);
		    var weight = 1.0 - cond;
		    dst(ij) = weight * src(ij) + 
		        (src(ij+up) + src(ij+down) + src(ij+left) + src(ij+right)) * (weight * dir) +
		        (src(ij+upleft) + src(ij+upright) + src(ij+downleft) + src(ij+downright)) * (weight * diag);
	    }
            dst[i,0] = dst[i,p.M];
            dst[i,p.M+1] = dst[i,1];
	}
	if (iteration % p.period == 0 || iteration == p.maxiter) {
	    r.niter = iteration;
	    var tmin = 1.0/0.0;
	    var tmax = -tmin;
	    var maxdiff = -tmin;
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
            r.maxdiff = maxdiff;
	    r.tavg = tavg / (p.N * p.M);
	    r.time = t.elapsed();

	    var done = (maxdiff < p.threshold || iteration == p.maxiter);
	    if (0) then report_results(p, r);
	    //if done then break;
	}
    }

    /* Stop timer, set final timing */
    t.stop();
    r.time = t.elapsed();

    return r;
}

