/* actually constants in the final code */
double dirnmul = sqrt(2)/(sqrt(2) + 1) / 4.0;
double diagnmul = 1.0/(sqrt(2) + 1) / 4.0;

double ourcnd = *cnd; // the point's conductivity
double leftover = 1.0 - ourcnd;
double directcnd = leftover * dirnmul; // direct neighbours
double diagcnd = leftover * diagnmul; // diagonal neighbours

double result = row[x] * ourcnd;
double directsum = rowup[x] + rowdown[x] + row[prev] + row[next];
result += directcnd * directsum;
double diagsum = rowup[prev] + rowup[next] + rowdown[prev] + rowdown[next];
result += diagcnd * diagsum;
