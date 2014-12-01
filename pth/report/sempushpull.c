/* idiom to send data to the next thread */
sem_wait(&send_bbuffer.empty);
push(tosend, &send_bbuffer);
sem_post(&send_bbuffer.full);

/* idiom to receive data from the previous thread */
sem_wait(&rec_bbuffer->full);
s = fetch(rec_bbuffer);
sem_post(&rec_bbuffer->empty);