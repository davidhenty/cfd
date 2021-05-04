void jacobistep(double **psinew, double **psi, int m, int n,
                int jstart, int jstop);

void jacobistepvort(double **zetnew, double **psinew,
		    double **zet,    double** psi,
		    int m, int n, double re);

double deltasq(double **newarr, double **oldarr, int m, int n);
