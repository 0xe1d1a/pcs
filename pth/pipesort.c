#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <time.h>
#define TRUE 1

typedef int symbol;
symbol *buffer;
size_t buffer_size;
int counter;

/* the semaphores */
sem_t full, empty;

/* push a symbol through the pipeline */
int push(symbol s) {
   /* When the buffer is not full add the item
      and increment the counter*/
   if(counter < buffer_size) {
      buffer[counter] = s;
      counter++;
      return TRUE;
   }
   /* Error the buffer is full */
   else {
      return !TRUE;
   }
}

int fetch(symbol *s) {
   /* When the buffer is not empty remove the item
      and decrement the counter */
   if(counter > 0) {
      *s = buffer[(counter-1)];
      counter--;
      return TRUE;
   }
   else { /* Error buffer empty */
      return !TRUE;
   }
}

void init() {
        counter = 0;
        /* Create simple seed from time */
        srand(time(NULL));
        /* Create the full semaphore and initialize to 0 */
        sem_init(&full, 0, 0);
        /* Create the empty semaphore and initialize to BUFFER_SIZE */
        sem_init(&empty, 0, buffer_size);
}

void *comparator() {
	while(TRUE) {
		/* aquire the full lock */
      		sem_wait(&full);
		symbol *s;
		int res = fetch(s);
		if(res) printf("Fetched symbol: %d\n", *s);
		else {
			printf("Buffer is empty\n");
			//break;
		}
		sem_post(&empty);
	}
} 



/* producer thread */ 
void generator() {
	init();
	pthread_t mthread;
	pthread_create(&mthread, NULL, comparator, NULL); //fire the first thread
	while(TRUE) {
		symbol s = rand();
		/* acquire the empty lock */
      		sem_wait(&empty);
		int res = push(s);
		if (res) printf("Generated symbol: %d\n", s);
		else { 
			printf("Buffer is full\n");
			//break;
		}
		sem_post(&full);
	}	
}


int main(int argc, char **argv) {
	if (argc != 2) {
		printf("Invalid params\n");
		return 1;
	}
	printf("Length is %s\n", argv[1]);	
	buffer_size = atoi(argv[1]);
	symbol tmp[buffer_size];
	buffer = tmp;
	
	generator();
}
