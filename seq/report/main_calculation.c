double result = row[x] * ourcnd;
double directsum = rowup[x] + rowdown[x] + row[prev] + row[next];
result += directcnd * directsum/4.0;
double diagsum = rowup[prev] + rowup[next] + rowdown[prev] + rowdown[next];
result += diagcnd * diagsum/4.0;