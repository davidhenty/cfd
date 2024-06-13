#include "boundary.h"
#include <stdio.h>

//grid is parallelised in the x direction

#ifndef USEVLA
void boundarypsi(int m, int n, double **psi, int b, int h, int w)
#else
void boundarypsi(int m, int n, double psi[m+2][n+2], int b, int h, int w)  
#endif
{

  int i,j;

  //BCs on bottom edge

  for (i=b+1;i<=b+w-1;i++)
    {
      psi[i][0] = (double)(i-b);
    }

  for (i=b+w;i<=m;i++)
    {
      psi[i][0] = (double)(w);
    }

  //BCS on RHS

  for (j=1; j <= h; j++)
    {
      psi[m+1][j] = (double) w;
    }

  for (j=h+1;j<=h+w-1; j++)
    {
      psi[m+1][j]=(double)(w-j+h);
    }
}

#ifndef USEVLA
void boundaryzet(int m, int n, double **zet, double **psi)
#else
void boundaryzet(int m, int n, double zet[m+2][n+2], double psi[m+2][n+2])
#endif  
{
  int i,j;

  //set top/bottom BCs:

  for (i=1;i<m+1;i++)
    {
      zet[i][0]   = 2.0*(psi[i][1]-psi[i][0]);
      zet[i][n+1] = 2.0*(psi[i][n]-psi[i][n+1]);
    }

  //set left BCs:

  for (j=1;j<n+1;j++)
    {
      zet[0][j] = 2.0*(psi[1][j]-psi[0][j]);
    }

  //set right BCs

  for (j=1;j<n+1;j++)
    {
      zet[m+1][j] = 2.0*(psi[m][j]-psi[m+1][j]);
    }
}
