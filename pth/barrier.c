#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h> // printf
#include <string.h> // memcpy
#include <float.h> // DBL
#include <errno.h>
#include <assert.h>
#include <pthread.h>
#include <getopt.h>
#include <semaphore.h>

struct semaphore_barrier {
	sem_t in_sem;
	sem_t out_sem;
	//sem_t mutex_sem;
	int count;
	int counter;
};

int semaphore_barrier_init(struct semaphore_barrier *barrier, const pthread_barrierattr_t *attr, unsigned int count) {
	sem_init(&barrier->in_sem, 0, 0);
	sem_init(&barrier->out_sem, 0, 0);
	//sem_init(&barrier->mutex_sem, 0, 1);
	barrier->counter = count;
	barrier->count = 0;

	return 0;
}

int semaphore_barrier_wait(struct semaphore_barrier *barrier) {
	int newcount;

	//sem_wait(&barrier->mutex_sem);
	newcount = __sync_fetch_and_add(&barrier->count, 1);
	if (newcount+1 == barrier->counter) {
		for (unsigned i = 0; i < barrier->counter; ++i)
			sem_post(&barrier->in_sem);
	}
	//sem_post(&barrier->mutex_sem);
	sem_wait(&barrier->in_sem);

	//sem_wait(&barrier->mutex_sem);
	newcount = __sync_fetch_and_add(&barrier->count, -1);
	if (newcount-1 == 0) {
		for (unsigned i = 0; i < barrier->counter; ++i)
			sem_post(&barrier->out_sem);
	}
	//sem_post(&barrier->mutex_sem);
	sem_wait(&barrier->out_sem);

	return 0;
}

struct mutexcv_barrier {
	pthread_mutex_t mutex;
	pthread_cond_t cond;
	unsigned count;
	unsigned counter;
};

int mutexcv_barrier_init(struct mutexcv_barrier *barrier, const pthread_barrierattr_t *attr, unsigned int count) {
	pthread_mutex_init(&barrier->mutex, NULL);
	pthread_cond_init(&barrier->cond, NULL);
	barrier->counter = count;
	barrier->count = 0;

	return 0;
}

int mutexcv_barrier_wait(struct mutexcv_barrier *barrier) {
	pthread_mutex_lock(&barrier->mutex);
	barrier->count++;
	if (barrier->count == barrier->counter) {
		barrier->count = 0;
		pthread_cond_broadcast(&barrier->cond);
	} else {
		pthread_cond_wait(&barrier->cond, &barrier->mutex);
	}
	pthread_mutex_unlock(&barrier->mutex);

	return 0;
}

struct barrier_funcs {
	// Return values are always zero, and attr is ignored.
	int (*barrier_init)(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count);
	int (*barrier_wait)(pthread_barrier_t *barrier);
};

struct barrier_data {
	const struct barrier_funcs funcs;
	union {
		pthread_barrier_t barrier;
		struct semaphore_barrier sem_barrier;
		struct mutexcv_barrier mut_barrier;
	} barrier;
	const unsigned count;
};

void *worker_func(struct barrier_data *data) {
	for (unsigned i = 0; i < data->count; ++i) {
		data->funcs.barrier_wait((void *)&data->barrier);
	}

	return NULL;
}

void usage(const char *pname) {
	printf("%s [-s|-m] (semaphore/mutex; defaults to posix barriers) [-p nthreads] [-c iterationcount]\n", pname);
}

int main(int argc, char **argv) {
	unsigned nthreads = 1;
	unsigned count = 1e8;
	struct barrier_funcs funcs;
	funcs.barrier_init = &pthread_barrier_init;
	funcs.barrier_wait = &pthread_barrier_wait;

	char ch;
	while ((ch = getopt(argc, argv, "p:c:sm")) != -1) {
		switch (ch) {
			case 'p': nthreads = strtol(optarg, 0, 10); break;
			case 'c': count = strtol(optarg, 0, 10); break;
			case 's':
				funcs.barrier_init = (void *)&semaphore_barrier_init;
				funcs.barrier_wait = (void *)&semaphore_barrier_wait;
				break;
			case 'm':
				funcs.barrier_init = (void *)&mutexcv_barrier_init;
				funcs.barrier_wait = (void *)&mutexcv_barrier_wait;
				break;
			case 'h': default: usage(argv[0]); return 1;
		}
	}

	struct barrier_data data = { funcs, {0}, count };

	// start timing (just accept the overhead)
	struct timeval tv_start, tv_curr, tv_diff;
	gettimeofday(&tv_start, NULL);

	funcs.barrier_init((void *)&data.barrier, NULL, nthreads);

	pthread_t threads[nthreads];
	for (unsigned int i = 0; i < nthreads; ++i) {
		pthread_create(&threads[i], NULL, (void *)&worker_func, &data);
	}

	for (unsigned int i = 0; i < nthreads; ++i) {
		pthread_join(threads[i], NULL);
	}

	gettimeofday(&tv_curr, NULL);
	timersub(&tv_curr, &tv_start, &tv_diff);
	double time_taken = tv_diff.tv_sec + (tv_diff.tv_usec / 1000000.0);

	printf("Execution time: %f\n", time_taken);
}
