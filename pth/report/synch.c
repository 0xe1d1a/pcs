void *thread_main(struct thread_params *tparams)
{
	//so we can kill them whenever we like!
	pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
	while (1)
	{
		pthread_barrier_wait(&barrier);
		thread_work(tparams); //computations
		pthread_barrier_wait(&barrier); // wait for computation
		if (g_do_reduction) {
			thread_reduction_work(tparams); //reduction calculations
			pthread_barrier_wait(&barrier); //wait for reductions
		}

		// swap our local copy
		double *tmp = tparams->t_prev;
		tparams->t_prev = tparams->t_next;
		tparams->t_next = tmp;
	}
}