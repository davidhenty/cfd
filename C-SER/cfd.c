#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "arraymalloc.h"
#include "boundary.h"
#include "jacobi.h"
#include "cfdio.h"

int main(int argc, char **argv)
{
  int printfreq=1000; //output frequency
  double error, bnorm;
  double tolerance=0.0; //tolerance for convergence. <=0 means do not check

  //main arrays
  double **psi, **zet;
  //temporary versions of main arrays
  double **psitmp, **zettmp;

  //command line arguments
  int scalefactor, numiter;

  double re; // Reynold's number - must be less than 3.7

  //simulation sizes
  int bbase=10;
  int hbase=15;
  int wbase=5;
  int mbase=32;
  int nbase=32;

  int irrotational = 1, checkerr = 0;

  int m,n,b,h,w;
  int iter;
  int i,j;

  double tstart, tstop, ttot, titer;

  //do we stop because of tolerance?
  if (tolerance > 0)
    {
      checkerr = 1;
    }

  //check command line parameters and parse them

  if (argc <3|| argc >4)
    {
      printf("Usage: cfd <scale> <numiter> [reynolds]\n");
      return 0;
    }

  scalefactor=atoi(argv[1]);
  numiter=atoi(argv[2]);

  if (argc == 4)
    {
      re=atof(argv[3]);
      irrotational=0;
    }
  else
    {
      re=-1.0;
    }

  if(!checkerr)
    {
      printf("Scale Factor = %i, iterations = %i\n",scalefactor, numiter);
    }
  else
    {
      printf("Scale Factor = %i, iterations = %i, tolerance= %g\n",scalefactor,numiter,tolerance);
    }

  if (irrotational)
    {
      printf("Irrotational flow\n");
    }
  else
    {
      printf("Reynolds number = %f\n",re);
    }

  //Calculate b, h & w and m & n
  b = bbase*scalefactor;
  h = hbase*scalefactor;
  w = wbase*scalefactor;
  m = mbase*scalefactor;
  n = nbase*scalefactor;

  re = re / (double)scalefactor;

  printf("Running CFD on %d x %d grid in serial\n",m,n);

  //allocate arrays

  psi    = (double **) arraymalloc2d(m+2,n+2,sizeof(double));
  psitmp = (double **) arraymalloc2d(m+2,n+2,sizeof(double));

  //zero the psi array
  for (i=0;i<m+2;i++)
    {
      for(j=0;j<n+2;j++)
	{
	  psi[i][j]=0.0;
	}
    }

  if (!irrotational)
    {
      //allocate arrays

      zet =   (double **) arraymalloc2d(m+2,n+2,sizeof(double));
      zettmp =(double **) arraymalloc2d(m+2,n+2,sizeof(double));

      //zero the zeta array

      for (i=0;i<m+2;i++)
	{
	  for(j=0;j<n+2;j++)
	    {
	      zet[i][j]=0.0;
	    }
	}
    }
  
  //set the psi boundary conditions

  boundarypsi(m,n,psi,b,h,w);

  //compute normalisation factor for error

  bnorm=0.0;

  for (i=0;i<m+2;i++)
    {
      for (j=0;j<n+2;j++)
	{
	  bnorm += psi[i][j]*psi[i][j];
	}
    }

  if (!irrotational)
    {
      //update zeta BCs that depend on psi
      boundaryzet(m,n,zet,psi);

      //update normalisation

      for (i=0;i<m+2;i++)
	{
	  for (j=0;j<n+2;j++)
	    {
	      bnorm += zet[i][j]*zet[i][j];
	    }
	}
    }

  bnorm=sqrt(bnorm);

  //begin iterative Jacobi loop

  printf("\nStarting main loop...\n\n");
  
  tstart=gettime();

  for(iter=1;iter<=numiter;iter++)
    {
      //calculate psi for next iteration

      if (irrotational)
	{
	  jacobistep(m,n,psitmp,psi);
	}
      else
	{
	  jacobistepvort(m,n,zettmp,psitmp,zet,psi,re);
	}

      //calculate current error if required

      if (checkerr || iter == numiter)
	{
	  error = deltasq(m,n,psitmp,psi);

	  if(!irrotational)
	    {
	      error += deltasq(m,n,zettmp,zet);
	    }

	  error=sqrt(error);
	  error=error/bnorm;
	}

      //copy back

      for(i=1;i<=m;i++)
	{
	  for(j=1;j<=n;j++)
	    {
	      psi[i][j]=psitmp[i][j];
	    }
	}

      if (!irrotational)
	{
	  for(i=1;i<=m;i++)
	    {
	      for(j=1;j<=n;j++)
		{
		  zet[i][j]=zettmp[i][j];
		}
	    }
	}

      if (!irrotational)
	{
	  //update zeta BCs that depend on psi
	  boundaryzet(m,n,zet,psi);
	}

      //quit early if we have reached required tolerance

      if (checkerr)
	{
	  if (error < tolerance)
	    {
	      printf("Converged on iteration %d\n",iter);
	      break;
	    }
	}

      //print loop information

      if(iter%printfreq == 0)
	{
	  if (!checkerr)
	    {
	      printf("Completed iteration %d\n",iter);
	    }
	  else
	    {
	      printf("Completed iteration %d, error = %g\n",iter,error);
	    }
	}
    }

  if (iter > numiter) iter=numiter;

  tstop=gettime();

  ttot=tstop-tstart;
  titer=ttot/(double)iter;


  //print out some stats

  printf("\n... finished\n");
  printf("After %d iterations, the error is %g\n",iter,error);
  printf("Time for %d iterations was %g seconds\n",iter,ttot);
  printf("Each iteration took %g seconds\n",titer);

  //output results

  writedatafiles(m,n,psi,scalefactor);

  writeplotfile(m,n,scalefactor);

  //free un-needed arrays
  free(psi);
  free(psitmp);

  if (!irrotational)
    {
      free(zet);
      free(zettmp);
    }

  printf("... finished\n");

  return 0;
}
