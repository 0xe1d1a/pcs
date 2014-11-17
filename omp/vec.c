
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>

int max_depth = 2;

/* Ordering of the vector */
typedef enum Ordering {ASCENDING, DESCENDING, RANDOM} Order;

int debug = 0;

// shared non-in-place merge function (merge entries from src into dest)
void TopDownMerge(int *src, int begin, int middle, int end, int *dest) {
	int i0 = begin;
	int i1 = middle;
	for (int j = begin; j < end; ++j) {
		if (i0 < middle && (i1 >= end || src[i0] <= src[i1])) {
			dest[j] = src[i0++];
		} else {
			dest[j] = src[i1++];
		}
	}
}

// sequential variant (deeper levels)
void TopDownSplitMergeSeq(int *src, int begin, int end, int *dest) {
	if (end - begin < 2)
		return;

	int middle = (end + begin)/2;
	TopDownSplitMergeSeq(src, begin, middle, dest);
	TopDownSplitMergeSeq(src, middle, end, dest);
	TopDownMerge(src, begin, middle, end, dest);
	memcpy(src + begin, dest + begin, sizeof(int) * (end - begin));
}

// parallelized variant (highest levels)
void TopDownSplitMerge(int *src, int begin, int end, int *dest, int depth) {
	if (end - begin < 2)
		return;

	// Once we've reached the desired depth, move to the sequential version.
	if (depth > max_depth) {
		TopDownSplitMergeSeq(src, begin, end, dest);
		return;
	}

	int middle = (end + begin)/2;
#pragma omp parallel sections num_threads(2)
{
#pragma omp section
	TopDownSplitMerge(src, begin, middle, dest, depth + 1);
#pragma omp section
	TopDownSplitMerge(src, middle, end, dest, depth + 1);
}
	// Merge the two sorted array into a sorted section of the temporary buffer,
	// then copy that temporary buffer back to our main results array.
	TopDownMerge(src, begin, middle, end, dest);
	memcpy(src + begin, dest + begin, sizeof(int) * (end - begin));
}

/* Sort vector v of l elements using mergesort */
void msort(int *v, long l) {
	int *tmpbuf = malloc(l * sizeof(int));
	TopDownSplitMerge(v, 0, l, tmpbuf, 0);
	free(tmpbuf);
}

void print_v(int *v, long l) {
  printf("\n");
  for(long i = 0; i < l; i++) {
    if(i != 0 && (i % 10 == 0)) {
      printf("\n");
    }
    printf("%d ", v[i]);
  }
  printf("\n");
}

int main(int argc, char **argv) {

  int c;
  int seed = 42;
  long length = 50;
  long maxlength = 1e6;
  long minlength = 1e5;
  Order order = ASCENDING;
  int **vector;
  int *vlengths;
  unsigned threadcount = -1;

  omp_set_nested(1);

  /* Read command-line options. */
  while((c = getopt(argc, argv, "n:adrgl:s:t:")) != -1) {
    switch(c) {
      case 'a':
        order = ASCENDING;
        break;
      case 'd':
        order = DESCENDING;
        break;
      case 'r':
        order = RANDOM;
        break;
      case 'l':
        length = atol(optarg);
        break;
      case 'n':
        max_depth = atol(optarg);
        break;
      case 't':
        threadcount = atol(optarg);
        // doesn't do the right thing for nesting, we used OMP_THREAD_LIMIT instead
//        omp_set_num_threads(threadcount);
        break;
      case 'g':
        debug = 1;
        break;
      case 's':
        seed = atoi(optarg);
      case '?':
        if(optopt == 'l' || optopt == 's') {
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        }
        else if(isprint(optopt)) {
          fprintf(stderr, "Unknown option '-%c'.\n", optopt);
        }
        else {
          fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
        }
        return -1;
      default:
        return -1;
      }
  }

  /* Seed such that we can always reproduce the same random vector */
  srand(seed);

  /* Allocate outer vector (and associated lengths). */
  vector = (int**)malloc(length*sizeof(int*));
  vlengths = (int *)malloc(length*sizeof(int));
  if(vector == NULL) {
    fprintf(stderr, "Malloc failed...\n");
    return -1;
  }

  for (long j = 0; j < length; j++) {
    /* Pick a random length. */
    long vlength = rand() % (maxlength - minlength) + minlength;
    vlengths[j] = vlength;
    /* Allocate inner vector. */
    vector[j] = (int*)malloc(vlength*sizeof(int));

    /* Fill vector. */
    switch(order){
      case ASCENDING:
        for(long i = 0; i < vlength; i++) {
          vector[j][i] = (int)i;
        }
        break;
      case DESCENDING:
        for(long i = 0; i < vlength; i++) {
          vector[j][i] = (int)(vlength - i);
        } 
        break;
      case RANDOM:
        for(long i = 0; i < vlength; i++) {
          vector[j][i] = rand();
        }
        break;
    }
  }

  struct timeval tv_start, tv_curr, tv_diff;
  gettimeofday(&tv_start, NULL);

  /* Sort */
  #pragma omp parallel for schedule(dynamic, 1)
// num_threads(2)
  for (long j = 0; j < length; j++) {
    if (debug)
      fprintf(stderr, "sorting #%d (size %d)\n", j, vlengths[j]);
    msort(vector[j], vlengths[j]);
  }

  gettimeofday(&tv_curr, NULL);
  timersub(&tv_curr, &tv_start, &tv_diff);
  float time = tv_diff.tv_sec + (tv_diff.tv_usec / 1000000.0);
  printf("%f %d all-tasks\n", time, (int)threadcount);

  return 0;
}

