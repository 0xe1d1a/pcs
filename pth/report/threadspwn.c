	for (int i=0; i<p->nthreads; i++)
	{
		params[i].srow = 1 + itersize*i; //thread starting row
		params[i].erow = 1 + itersize*(i+1); //thread ending row 
		if (i == p->nthreads-1)	//special case for bounded row
			params[i].erow = p->N+1;
		params[i].t_prev = t_prev; //own copy of previous snapshot
		params[i].t_next = t_next; //own copy of next snapshot
		pthread_create(&threads[i], NULL, &thread_main, &params[i]);
	}

	while (TRUE) { 
		pthread_barrier_wait(&barrier);
		/* master thread operates here to control reductions */
	}