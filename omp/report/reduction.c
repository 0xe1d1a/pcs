#pragma omp parallel for reduction(min: tmin) reduction(max: tmax) \
                         reduction(max: maxdiff) reduction(+: tavg)
for (size_t y = 1; y < p->N+1; ++y)
{
    /*
        reduction computation here
    */
}