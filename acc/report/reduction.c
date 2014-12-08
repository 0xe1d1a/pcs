// iterate over non-constant rows
#pragma acc parallel loop gang independent reduction(min: tmin)
 reduction(max: tmax) reduction(max: maxdiff) 
 reduction(+: tavg) present(dst[0:w*h], src[0:w*h])
for (size_t y = 1; y < h - 1; ++y) {
	#pragma acc loop vectors independent
	for (size_t x = 1; x < w - 1; ++x) {
		/* reduction computation here */
	}
}
