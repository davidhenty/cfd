#ifndef USEVLA
void jacobistep(int m, int n, double **psinew, double **psi);
#else
void jacobistep(int m, int n, double psinew[m+2][n+2], double psi[m+2][n+2]);
#endif

#ifndef USEVLA
void jacobistepvort(int m, int n,
                    double **zetnew, double **psinew,
		    double **zet,    double **psi,
		    double re);
#else
void jacobistepvort(int m, int n,
                    double zetnew[m+2][n+2], double psinew[m+2][n+2],
		    double zet[m+2][n+2],    double psi[m+2][n+2],
		    double re);
#endif

#ifndef USEVLA
double deltasq(int m, int n, double **newarr, double **oldarr);
#else
double deltasq(int m, int n, double newarr[m+2][n+2], double oldarr[m+2][n+2]);
#endif
