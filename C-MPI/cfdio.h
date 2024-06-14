#include <mpi.h>

#ifndef USEVLA
void writedatafiles( int m, int n, double **psi, int scale,  MPI_Comm comm);
#else
void writedatafiles( int m, int n, double psi[m+2][n+2], int scale,  MPI_Comm comm);
#endif

void writeplotfile(int m, int n, int scale);

void hue2rgb(double hue, int *r, int *g, int *b);

double colfunc(double x);
