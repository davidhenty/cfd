#ifndef USEVLA
void writedatafiles( int m, int n, double **psi, int scale);
#else
void writedatafiles( int m, int n, double psi[m][n], int scale);
#endif

void writeplotfile(int m, int n, int scale);

void hue2rgb(double hue, int *r, int *g, int *b);

double colfunc(double x);

double gettime(void);
