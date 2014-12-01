#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <time.h>
#define TRUE 1

typedef int symbol;

// wrap-around/FIFO data structure
typedef struct {
	symbol *buffer;
	int nextin;
	int nextout;
	size_t size;
	sem_t full, empty;
	int buf_id;
} bounded_buffer;

// global thread ID (for debugging), safe because we add the threads only after making the buffers
int thread_id = 0;

// function to create a new bounded buffer (for between threads)
void init_bbuffer(bounded_buffer *bbuffer, size_t length) {
	bbuffer->buffer = (int *) malloc(length * sizeof(int));
	bbuffer->size = length;
	bbuffer->nextin = 0;
	bbuffer->nextout = 0;
	sem_init(&bbuffer->full, 0, 0);
	sem_init(&bbuffer->empty, 0, bbuffer->size);
	bbuffer->buf_id = thread_id;
}

// push a symbol into the bounded buffer
void push(symbol s, bounded_buffer *send_bbuffer) {
	send_bbuffer->buffer[send_bbuffer->nextin] = s;
	send_bbuffer->nextin = (send_bbuffer->nextin+1) % send_bbuffer->size;
}

// fetch a symbol from the bounded buffer
symbol fetch(bounded_buffer *rec_bbuffer) {
	symbol s  = rec_bbuffer->buffer[rec_bbuffer->nextout];
	rec_bbuffer->nextout = (rec_bbuffer->nextout + 1) % rec_bbuffer->size;
	return s;
}

// printer thread (at the end)
void *printer(void *p_buffer) {
	bounded_buffer *rec_bbuffer = (bounded_buffer*) p_buffer;
	thread_id++;
	symbol s;
	while (TRUE) {
		// get next symbol
		sem_wait(&rec_bbuffer->full);
		s = fetch(rec_bbuffer);
		sem_post(&rec_bbuffer->empty);
		// not an END symbol, print the value
		if (s!= -1 && s!=-2) printf("%d\n", s);
		// second END symbol; this causes all the pthread_join calls to succeed in order
		if (s == -2) return NULL;
	}
}

// comparator (normal) thread
void *comparator(void *p_buffer) {
	bounded_buffer *rec_bbuffer = (bounded_buffer*) p_buffer;
	symbol myVal = -3; // special 'no value seen yet' value
	symbol s;
	symbol tosend;
	pthread_t next;
	thread_id++;
	bounded_buffer send_bbuffer;
	int have_successor = !TRUE;	
	while (TRUE) {
		// get the next value from our incoming buffer..
		sem_wait(&rec_bbuffer->full);
		s = fetch(rec_bbuffer);
		sem_post(&rec_bbuffer->empty);
		if (myVal == -3) {
			myVal = s;
			continue;
		}
		if (s == -1 ) { //ending state (first END)
			if (have_successor == !TRUE) {
				// (just in case we're the last thread)
				have_successor = TRUE;
				init_bbuffer(&send_bbuffer, rec_bbuffer->size);
				pthread_create(&next, NULL, printer, (void*) &send_bbuffer);
			}
			// push the end symbol plus our symbol
			sem_wait(&send_bbuffer.empty);
			push(s, &send_bbuffer);
			sem_post(&send_bbuffer.full);
			sem_wait(&send_bbuffer.empty);
			push(myVal, &send_bbuffer);
			sem_post(&send_bbuffer.full);
			while(TRUE) {
				// keep pushing until we're done
				sem_wait(&rec_bbuffer->full);
				s = fetch(rec_bbuffer);
				sem_post(&rec_bbuffer->empty);
				sem_wait(&send_bbuffer.empty);
				push(s, &send_bbuffer);
				sem_post(&send_bbuffer.full);
				if(s == -2) { //kill state (second END)
					// join the next thread in sequence and return
					pthread_join(next, NULL);
					return NULL;
				}
			}
		}
		else {
			//general state; pick largest value
			if (myVal > s) {
				tosend = s;
			}
			else {
				tosend = myVal;
				myVal = s;
			}
			if (have_successor == !TRUE) {
				// we don't have a successor yet, make one!
				have_successor = TRUE;
				init_bbuffer(&send_bbuffer, rec_bbuffer->size);
				pthread_create(&next, NULL, comparator, (void*) &send_bbuffer);
			}
			// push other symbol to successor
			sem_wait(&send_bbuffer.empty);
			push(tosend, &send_bbuffer);
			sem_post(&send_bbuffer.full);
		}
	}
}

void pipesort(size_t length, size_t bufsize) {
	// make a buffer for communication with the thread we're about to make
	bounded_buffer send_bbuffer;
	init_bbuffer(&send_bbuffer, bufsize);

	pthread_t next;
	// make a single comparator thread to begin
	pthread_create(&next, NULL, comparator, (void*) &send_bbuffer);

	printf("Generating %d numbers (bufsize %d)\n", (int)length, (int)bufsize);
	// seed with some arbitrary value to get reproducible results (ok, it shouldn't matter)
	srand(5678);
	int i;
	for (i=0; i<length; i++) {
		// generate a symbol
		symbol s = rand() % 5000;
#ifdef PIPEDEBUG
		printf("Number: %d\n", s);
#endif
		// push into the first buffer
		sem_wait(&send_bbuffer.empty);
		push(s, &send_bbuffer);
		sem_post(&send_bbuffer.full);
	}

	// send two END symbols into the buffer
	printf("Sending END symbol 1\n");
	sem_wait(&send_bbuffer.empty);
	push(-1, &send_bbuffer);
	sem_post(&send_bbuffer.full);
	printf("Sending END symbol 2\n");
	sem_wait(&send_bbuffer.empty);
	push(-2, &send_bbuffer);
	sem_post(&send_bbuffer.full);

	// wait for the first comparator thread to die before finishing
	pthread_join(next, NULL);
}

int main(int argc, char **argv) {
	if (argc != 3) {
		printf("Invalid params\n");
		return 1;
	}
	pipesort(atoi(argv[1]), atoi(argv[2]));
	return 0;
}

