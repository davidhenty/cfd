#ifndef USEVLA
void boundarypsi(int m, int n, double **psi, int b, int h, int w,
		 MPI_Comm comm);

void boundaryzet(int m, int n, double **zet, double **psi, MPI_Comm comm);

void haloswap(int m, int  n, double **x, MPI_Comm comm);
#else
void boundarypsi(int m, int n, double psi[m+2][n+2], int b, int h, int w,
		 MPI_Comm comm);

void boundaryzet(int m, int n, double zet[m+2][n+2], double psi[m+2][n+2], MPI_Comm comm);

void haloswap(int m, int  n, double x[m+2][n+2], MPI_Comm comm);
#endif
