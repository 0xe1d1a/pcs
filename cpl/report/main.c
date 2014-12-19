for iteration in 1..p.maxiter {
	dst <=> src;	//swap buffers
	forall i in 1..p.N {			
		for j in 1..p.M {
		    /* computation here */
		}	
	}

	forall i in 1..p.N {
		dst[i,0] = dst[i,p.M];
		dst[i,p.M+1] = dst[i,1];
	}