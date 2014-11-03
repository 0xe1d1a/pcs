// first column
do_calc(cnd_dst++, dst++, 0, p->M-1, 1, row, rowup, rowdown);
// traverse the remaining columns in the row
for (size_t x = 1; x < p->M-1; ++x) {
	do_calc(cnd_dst++, dst++, x, x-1, x+1, row, rowup, rowdown);
}
// last column
do_calc(cnd_dst++, dst++, p->M-1, p->M-2, 0, row, rowup, rowdown);
