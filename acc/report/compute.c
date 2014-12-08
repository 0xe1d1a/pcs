#pragma acc parallel loop gang independent
 present(c[0:w*h], dst[0:w*h], src[0:w*h])
for (size_t y = 1; y < h - 1; ++y)
{
	#pragma acc loop worker independent
	for (size_t x = 1; x < w - 1; ++x)
	{
		/* computation */
	}

	/* smearing of borders */
}
