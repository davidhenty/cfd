void boundarypsi(int m, int n, double **psi, int b, int h, int w,
		 MPI_Comm comm);

void boundaryzet(int m, int n, double **zet, double **psi, MPI_Comm comm);

void haloswap(int m, int  n, double **x, MPI_Comm comm);
