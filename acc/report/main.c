#pragma acc enter data copyin(c[0:w*h], dst[0:w*h], src[0:w*h])
for (size_t iter = 1; iter <= p->maxiter; ++iter)
{
	/* per-iteration code */
}
