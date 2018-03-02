#include <mpi.h>

void writedatafiles(double **psi, int m, int n, int scale, MPI_Comm comm);

void writeplotfile(int m, int n, int scale);

void hue2rgb(double hue, int *r, int *g, int *b);

double colfunc(double x);
