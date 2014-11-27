#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <time.h>
#define TRUE 1

typedef int symbol;

typedef struct {
	symbol *buffer;
	int nextin;
	int nextout;
	size_t size;
	sem_t full, empty;
	int buf_id;
	//pthread_mutex_t mutex;
} bounded_buffer;

int thread_id = 0;

void init_bbuffer(bounded_buffer *bbuffer, size_t length) {
	bbuffer->buffer = (int *) malloc(length * sizeof(int));
	bbuffer->size = length;
	bbuffer->nextin = 0;
	bbuffer->nextout = 0;
	sem_init(&bbuffer->full, 0, 0);
	sem_init(&bbuffer->empty, 0, bbuffer->size);
	bbuffer->buf_id = thread_id;
}


void push(symbol s, bounded_buffer *send_bbuffer) {
	send_bbuffer->buffer[send_bbuffer->nextin] = s;
	send_bbuffer->nextin = (send_bbuffer->nextin+1) % send_bbuffer->size;
}

symbol fetch(bounded_buffer *rec_bbuffer) {
	symbol s  = rec_bbuffer->buffer[rec_bbuffer->nextout];
	rec_bbuffer->nextout = (rec_bbuffer->nextout + 1) % rec_bbuffer->size;
	return s;
}

void *printer(void *p_buffer) {
	bounded_buffer *rec_bbuffer = (bounded_buffer*) p_buffer;
	int my_id = thread_id++;
	symbol s;
	while (TRUE) {
		sem_wait(&rec_bbuffer->full);
		s = fetch(rec_bbuffer);
		sem_post(&rec_bbuffer->empty);
		if (s!= -1 && s!=-2) printf("%d\n", s);
		if (s == -2) return;
	}
}

void *comparator(void *p_buffer) {
	bounded_buffer *rec_bbuffer = (bounded_buffer*) p_buffer;
	symbol myVal = 501;
	symbol s;
	symbol tosend;	
	pthread_t next;
	int my_id = thread_id++;
	bounded_buffer send_bbuffer;
	int have_successor = !TRUE;	
	while (TRUE) {
		sem_wait(&rec_bbuffer->full);
		//pthread_mutex_lock(&rec_bbuffer->mutex);
		s = fetch(rec_bbuffer);
		//pthread_mutex_unlock(&rec_bbuffer->mutex);
		sem_post(&rec_bbuffer->empty);
		if (myVal == 501) {
			myVal = s;
			//printf("Thread: %u - My val is %d nothing to send\n", my_id, myVal);
			continue;
		}
		if (s == -1 ) {
			if(have_successor == !TRUE) {
				have_successor = TRUE;
				init_bbuffer(&send_bbuffer, rec_bbuffer->size);
				pthread_create(&next, NULL, printer, (void*) &send_bbuffer);
			}
			sem_wait(&send_bbuffer.empty);
			push(s, &send_bbuffer);
			sem_post(&send_bbuffer.full);
			sem_wait(&send_bbuffer.empty);
			push(myVal, &send_bbuffer);
			sem_post(&send_bbuffer.full);
			while(TRUE) {
				sem_wait(&rec_bbuffer->full);
				s = fetch(rec_bbuffer);
				sem_post(&rec_bbuffer->empty);
				sem_wait(&send_bbuffer.empty);
				push(s, &send_bbuffer);
				sem_post(&send_bbuffer.full);
				if(s == -2) {
					pthread_join(next, NULL);
					return NULL;
				}
			}
		}
		else {
			if (myVal > s) {
				tosend = s;
			}
			else {
				tosend = myVal;
				myVal = s;
			}
			//printf("Thread: %u - My val is %d, to send is %d\n",my_id, myVal, tosend);
			if(have_successor == !TRUE) {
				have_successor = TRUE;
				init_bbuffer(&send_bbuffer, rec_bbuffer->size);
				pthread_create(&next, NULL, comparator, (void*) &send_bbuffer);
			}
			sem_wait(&send_bbuffer.empty);
			//pthread_mutex_lock(&send_bbuffer.mutex);
			push(tosend, &send_bbuffer);
			//pthread_mutex_unlock(&send_bbuffer.mutex);
			sem_post(&send_bbuffer.full);
		}
	}
}

void pipesort(size_t length) {
	bounded_buffer send_bbuffer;
	init_bbuffer(&send_bbuffer, length);

	pthread_t next;
	pthread_create(&next, NULL, comparator, (void*) &send_bbuffer);

	printf("Generating %d numbers\n", length);
	srand(time(NULL));
	int i;
	for (i=0; i<length; i++) {
		symbol s = rand() % 500;
		printf("Number: %d\n", s);
		sem_wait(&send_bbuffer.empty);
		//pthread_mutex_lock(&send_bbuffer.mutex);
		push(s, &send_bbuffer);
		//pthread_mutex_unlock(&send_bbuffer.mutex);
		sem_post(&send_bbuffer.full);
	}
	//sleep(1);
	printf("Sending END symbol 1\n");
	sem_wait(&send_bbuffer.empty);
	push(-1, &send_bbuffer);
	sem_post(&send_bbuffer.full);
	printf("Sending END symbol 2\n");
	sem_wait(&send_bbuffer.empty);
	push(-2, &send_bbuffer);
	sem_post(&send_bbuffer.full);

	pthread_join(next, NULL);
}

int main(int argc, char **argv) {
	if (argc != 2) {
		printf("Invalid params\n");
		return 1;
	}
	pipesort(atoi(argv[1]));
	return 0;
}

