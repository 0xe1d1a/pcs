
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <omp.h>

/* Ordering of the vector */
typedef enum Ordering {ASCENDING, DESCENDING, RANDOM} Order;

int debug = 0;

void TopDownMerge(int *src, int begin, int middle, int end, int *dest) {
	int i0 = begin;
	int i1 = end;
	for (int j = begin; j < end; ++j) {
		if (i0 < middle && (i1 >= end || src[i0] <= src[i1])) {
			dest[j] = src[i0++];
			i0++;
		} else {
			dest[j] = src[i1++];
		}
	}
}

void TopDownSplitMerge(int *src, int begin, int end, int *dest) {
	if (end - begin < 2)
		return;

	int middle = (end + begin)/2;
	TopDownSplitMerge(src, begin, middle, dest);
	TopDownSplitMerge(src, middle, end, dest);
	TopDownMerge(src, begin, middle, end, dest);
	memcpy(src + begin, dest + begin, sizeof(int) * (end - begin));
}

/* Sort vector v of l elements using mergesort */
void msort(int *v, long l){
	int *tmpbuf = malloc(l * sizeof(int));
	TopDownSplitMerge(v, 0, l, tmpbuf);
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
  long length = 1e4;
  Order order = ASCENDING;
  int *vector;

  /* Read command-line options. */
  while((c = getopt(argc, argv, "adrgl:s:")) != -1) {
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

  /* Allocate vector. */
  vector = (int*)malloc(length*sizeof(int));
  if(vector == NULL) {
    fprintf(stderr, "Malloc failed...\n");
    return -1;
  }

  /* Fill vector. */
  switch(order){
    case ASCENDING:
      for(long i = 0; i < length; i++) {
        vector[i] = (int)i;
      }
      break;
    case DESCENDING:
      for(long i = 0; i < length; i++) {
        vector[i] = (int)(length - i);
      } 
      break;
    case RANDOM:
      for(long i = 0; i < length; i++) {
        vector[i] = rand();
      }
      break;
  }

  if(debug) {
    print_v(vector, length);
  }

  /* Sort */
  msort(vector, length);

  if(debug) {
    print_v(vector, length);
  }

  return 0;
}

