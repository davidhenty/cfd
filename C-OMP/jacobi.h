void jacobistep(int m, int n, double **psinew, double **psi);

void jacobistepvort(int m, int n,
                    double **zetnew, double **psinew,
		    double **zet,    double** psi,
		    double re);

double deltasq(int m, int n, double **newarr, double **oldarr);
