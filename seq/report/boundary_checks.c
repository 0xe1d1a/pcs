const double *row = &t_prev[p->M * y];
// boundary conditions (from the initial state)
const double *rowup = &t_prev[p->M * (y - 1)];
if (y == 0)
	rowup = &p->tinit[0];
const double *rowdown = &t_prev[p->M * (y + 1)];
if (y == p->N-1)
	rowdown = &p->tinit[p->M * (p->N-1)];
