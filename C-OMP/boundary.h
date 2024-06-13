#ifndef USEVLA
void boundarypsi(int m, int n, double **psi, int b, int h, int w);
void boundaryzet(int m, int n, double **zet, double **psi);
#else
void boundarypsi(int m, int n, double psi[m+2][n+2], int b, int h, int w);
void boundaryzet(int m, int n, double psi[m+2][n+2], double zet[m+2][n+2]);
#endif
