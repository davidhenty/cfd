#include <stdio.h>

#include "jacobi.h"

void jacobistep(int m, int n, double **psinew, double **psi)
{
  int i, j;

#pragma omp parallel for default(none), private(i,j), shared(psi,psinew,m,n)
  for(i=1;i<=m;i++)
    {
      for(j=1;j<=n;j++)
	{
	  psinew[i][j]=0.25*(psi[i-1][j]+psi[i+1][j]+psi[i][j-1]+psi[i][j+1]);
        }
    }
}

void jacobistepvort(int m, int n,
                    double **zetnew, double **psinew,
		    double **zet,    double **psi,
		    double re)
{
  int i, j;

#pragma omp parallel for default(none), \
  private(i,j), shared(psi,psinew,zet,m,n)
  for(i=1;i<=m;i++)
    {
      for(j=1;j<=n;j++)
	{
	  psinew[i][j]=0.25*(  psi[i-1][j]+psi[i+1][j]+psi[i][j-1]+psi[i][j+1]
			     - zet[i][j] );
	}
    }

#pragma omp parallel for default(none), \
  private(i,j), shared(zet,zetnew,psi,m,n,re)
  for(i=1;i<=m;i++)
    {
      for(j=1;j<=n;j++)
	{
	  zetnew[i][j]=0.25*(zet[i-1][j]+zet[i+1][j]+zet[i][j-1]+zet[i][j+1])
	    - re/16.0*(
		       (  psi[i][j+1]-psi[i][j-1])*(zet[i+1][j]-zet[i-1][j])
		       - (psi[i+1][j]-psi[i-1][j])*(zet[i][j+1]-zet[i][j-1])
		       );
	}
    }
}

double deltasq(int m, int n, double **newarr, double **oldarr)
{
  int i, j;

  double dsq=0.0;
  double tmp;

#pragma omp parallel for default(none), \
  private(i,j,tmp), shared(newarr,oldarr,m,n), reduction(+:dsq)
  for(i=1;i<=m;i++)
    {
      for(j=1;j<=n;j++)
	{
	  tmp = newarr[i][j]-oldarr[i][j];
	  dsq += tmp*tmp;
        }
    }

  return dsq;
}
