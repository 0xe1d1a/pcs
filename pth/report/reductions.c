for (int i=0; i<p->nthreads; i++) {
	if (params[i].tmin < r->tmin)
		r->tmin = params[i].tmin;
	if (params[i].tmax > r->tmax)
		r->tmax = params[i].tmax;
	if (params[i].maxdiff > r->maxdiff)
		r->maxdiff = params[i].maxdiff;
	r->tavg += params[i].tavg;
}
r->tavg = r->tavg / (p->N * p->M);
