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

unsigned int total_iters = 100;
unsigned int time_in_section = 1e7;
unsigned int time_doing_work = 1e7;
unsigned int write_proportion = 1;
unsigned int mem_size = 8;
char *mem_ptr = 0, *mem_ptr2 = 0;

__attribute__((transaction_safe)) int wait_in_barrier() {
	volatile int t = 0;
	for (unsigned i = 0; i < time_in_section; ++i) {
		t++;
	}
	return t;
}

pthread_mutex_t mymutex = PTHREAD_MUTEX_INITIALIZER;
int mutexlock_access(int write) {
	pthread_mutex_lock(&mymutex);
	wait_in_barrier();
	pthread_mutex_unlock(&mymutex);
	return 0;
}

sem_t mysem;
int semaphore_access(int write) {
	sem_wait(&mysem);
	wait_in_barrier();
	sem_post(&mysem);
	return 0;
}

pthread_rwlock_t myrwlock = PTHREAD_RWLOCK_INITIALIZER;
int rwlock_access(int write) {
	if (write)
		pthread_rwlock_wrlock(&myrwlock);
	else
		pthread_rwlock_rdlock(&myrwlock);
	wait_in_barrier();
	pthread_rwlock_unlock(&myrwlock);
	return 0;
}

int countdown(int from) {
	if (from == 0) return 0;
	return countdown(from-1);
}

int stm_access(int write) {
	int ret = 0;
	__transaction_atomic {
		if (write)
			memcpy(mem_ptr2, mem_ptr, mem_size);
		else if (*mem_ptr2) // force read
			ret = countdown(500);
		wait_in_barrier();
	}
	return ret;
}

int (*do_access)(int write) = &mutexlock_access;

void do_work() {
	volatile int t = 0;
	for (unsigned i = 0; i < time_doing_work; ++i) {
		t++;
	}
}

// this is the main function for all the threads
void *worker_func() {
	unsigned int seed = time(NULL);
	for (unsigned i = 0; i < total_iters; ++i) {
		do_work();
		// write_proportion = number of reads done for every write
		// (e.g., if 1, then we do 1 read for every write, approximately)
		// (and 0 = always write)
		unsigned check = rand_r(&seed);
		int write = ((check % (write_proportion + 1)) == 0);
		do_access(write);
	}

	return NULL;
}

void usage(const char *pname) {
	printf("%s [-p nthreads] [-c iterationcount] [-t timeforwork] [-b timeinsection] [-r writeproportion] [-s memsize] [-S|-R|-T] (semaphores/rwlock/STM)\n", pname);
}

int main(int argc, char **argv) {
	unsigned nthreads = 1;
	
	char ch;
	while ((ch = getopt(argc, argv, "p:c:t:b:r:s:SRT")) != -1) {
		switch (ch) {
			case 'p': nthreads = strtol(optarg, 0, 10); break;
			case 'c': total_iters = strtol(optarg, 0, 10); break;
			case 't': time_doing_work = strtol(optarg, 0, 10); break;
			case 'b': time_in_section = strtol(optarg, 0, 10); break;
			case 'r': write_proportion = strtol(optarg, 0, 10); break;
			case 's': mem_size = strtol(optarg, 0, 10); break;
			case 'S': do_access = &semaphore_access; break;
			case 'R': do_access = &rwlock_access; break;
			case 'T': do_access = &stm_access; break;
			case 'h': default: usage(argv[0]); return 1;
		}
	}

	mem_ptr = malloc(mem_size);
	mem_ptr2 = malloc(mem_size);
	sem_init(&mysem, 0, 1);

	printf("%d threads, %d iters, %d time_doing_work, %d time_in_section, %d write proportion, %d buffer size\n",
		nthreads, total_iters, time_doing_work, time_in_section, write_proportion, mem_size);

	// start timing (just accept the overhead)
	struct timeval tv_start, tv_curr, tv_diff;
	gettimeofday(&tv_start, NULL);

	pthread_t threads[nthreads];
	for (unsigned int i = 0; i < nthreads; ++i) {
		pthread_create(&threads[i], NULL, (void *)&worker_func, NULL);
	}

	// wait for all the threads to finish
	for (unsigned int i = 0; i < nthreads; ++i) {
		pthread_join(threads[i], NULL);
	}

	gettimeofday(&tv_curr, NULL);
	timersub(&tv_curr, &tv_start, &tv_diff);
	double time_taken = tv_diff.tv_sec + (tv_diff.tv_usec / 1000000.0);

	// tidy up (pointlessly)
	sem_destroy(&mysem);
	free(mem_ptr);
	free(mem_ptr2);

	printf("Execution time: %f\n", time_taken);
}
