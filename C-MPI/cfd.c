#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "arraymalloc.h"
#include "boundary.h"
#include "jacobi.h"
#include "cfdio.h"

int main(int argc, char **argv)
{
  int printfreq=1000; //output frequency
  double localerror, error, localbnorm, bnorm;
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

  int m,n,lm,b,h,w;
  int iter;
  int i,j;

  double tstart, tstop, ttot, titer;

  //parallelisation parameters
  int rank, size;
  MPI_Comm comm;


  //do we stop because of tolerance?
  if (tolerance > 0)
    {
      checkerr = 1;
    }

  comm=MPI_COMM_WORLD;

  MPI_Init(&argc,&argv);

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  //check command line parameters and parse them

  if (argc <3|| argc >4)
    {
      if (rank == 0) printf("Usage: cfd <scale> <numiter> [reynolds]\n");
      MPI_Finalize();
      return 0;
    }

  if (rank == 0)
    {
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
    }


  //broadcast runtime params to other processors
  MPI_Bcast(&scalefactor,1,MPI_INT,0,comm);
  MPI_Bcast(&numiter,1,MPI_INT,0,comm);
  MPI_Bcast(&re,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(&irrotational,1,MPI_INT,0,comm);

  //Calculate b, h & w and m & n
  b = bbase*scalefactor;
  h = hbase*scalefactor;
  w = wbase*scalefactor;
  m = mbase*scalefactor;
  n = nbase*scalefactor;

  re = re / (double)scalefactor;

  //calculate local size
  lm = m/size;

  //consistency check
  if (size*lm != m)
    {
      if (rank == 0)
        {
          printf("ERROR: m=%d does not divide onto %d processes\n", m, size);
        }
      MPI_Finalize();
      return -1;
    }

  if (rank == 0)
    {
      printf("Running CFD on %d x %d grid using %d process(es) \n",m,n,size);
    }

  //allocate arrays

  psi    = (double **) arraymalloc2d(lm+2,n+2,sizeof(double));
  psitmp = (double **) arraymalloc2d(lm+2,n+2,sizeof(double));

  //zero the psi array
  for (i=0;i<lm+2;i++)
    {
      for(j=0;j<n+2;j++)
        {
          psi[i][j]=0.;
        }
    }

  if (!irrotational)
    {
      //allocate arrays

      zet =   (double **) arraymalloc2d(lm+2,n+2,sizeof(double));
      zettmp =(double **) arraymalloc2d(lm+2,n+2,sizeof(double));

      //zero the zeta array

      for (i=0;i<lm+2;i++)
        {
          for(j=0;j<n+2;j++)
            {
              zet[i][j]=0.;
            }
        }
    }
  
  //set the psi boundary conditions

  boundarypsi(psi,lm,n,b,h,w,comm);

  //compute normalisation factor for error

  localbnorm=0.;

  for (i=0;i<lm+2;i++)
    {
      for (j=0;j<n+2;j++)
        {
          localbnorm += psi[i][j]*psi[i][j];
        }
    }

  //boundary swap of psi

  haloswap(psi,lm,n,comm);

  if (!irrotational)
    {
      //update zeta BCs that depend on psi
      boundaryzet(zet,psi,lm,n,comm);

      //update normalisation

      for (i=0;i<lm+2;i++)
        {
          for (j=0;j<n+2;j++)
            {
              localbnorm += zet[i][j]*zet[i][j];
            }
        }

      //boundary swap zeta
      haloswap(zet,lm,n,comm);
    }

  //get global bnorm
  MPI_Allreduce(&localbnorm,&bnorm,1,MPI_DOUBLE,MPI_SUM,comm);

  bnorm=sqrt(bnorm);

  //begin iterative Jacobi loop

  if (rank == 0)
    {
      printf("\nStarting main loop...\n\n");
    }

  //barrier for accurate timing - not needed for correctness

  MPI_Barrier(comm);

  tstart=MPI_Wtime();

  for(iter=1;iter<=numiter;iter++)
    {
      if (irrotational)
        {
          jacobistep(psitmp,psi,lm,n);
        }
      else
        {
          jacobistepvort(zettmp,psitmp,zet,psi,lm,n,re);
        }

      //calculate current error if required

      if (checkerr || iter == numiter)
        {
          localerror = deltasq(psitmp,psi,lm,n);

          if(!irrotational)
            {
              localerror += deltasq(zettmp,zet,lm,n);
            }

          MPI_Allreduce(&localerror,&error,1,MPI_DOUBLE,MPI_SUM,comm);
          error=sqrt(error);
          error=error/bnorm;
        }

      //copy back

      for(i=1;i<=lm;i++)
        {
          for(j=1;j<=n;j++)
            {
              psi[i][j]=psitmp[i][j];
            }
        }

      if (!irrotational)
        {
          for(i=1;i<=lm;i++)
            {
              for(j=1;j<=n;j++)
                {
                  zet[i][j]=zettmp[i][j];
                }
            }
        }

      //do a boundary swap

      haloswap(psi,lm,n,comm);

      if (!irrotational)
        {
          haloswap(zet,lm,n,comm);
          boundaryzet(zet,psi,lm,n,comm);
        }

      //quit early if we have reached required tolerance

      if (checkerr)
        {
          if (error < tolerance)
            {
              if (rank == 0)
                {
                  printf("Converged on iteration %d\n",iter);
                }
              break;
            }
        }

      //print loop information

      if(iter%printfreq == 0)
        {
          if (rank==0)
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
    }

  if (iter > numiter) iter=numiter;

  MPI_Barrier(comm);
  tstop=MPI_Wtime();

  ttot=tstop-tstart;
  titer=ttot/(double)iter;


  //print out some stats
  if (rank == 0)
    {
      printf("\n... finished\n");
      printf("After %d iterations, the error is %g\n",iter,error);
      printf("Time for %d iterations was %g seconds\n",iter,ttot);
      printf("Each iteration took %g seconds\n",titer);
    }

  //output results

  writedatafiles(psi,lm,n, scalefactor,comm);

  if (rank == 0) writeplotfile(m,n,scalefactor);

  //free un-needed arrays
  free(psi);
  free(psitmp);

  if (!irrotational)
    {
      free(zet);
      free(zettmp);
    }

  MPI_Finalize();

  if (rank == 0)
    {
      printf("... finished\n");
    }

  return 0;
}
