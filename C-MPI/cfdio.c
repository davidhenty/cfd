#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cfdio.h"
#include "arraymalloc.h"

#ifndef USEVLA
void writedatafiles(int m, int n, double **psi, int scale, MPI_Comm comm)
#else
void writedatafiles(int m, int n, double psi[m+2][n+2], int scale, MPI_Comm comm)
#endif
{
  typedef double vecvel[2];
  typedef int    vecrgb[3];

  vecvel **vel;
  vecrgb **rgb;

  double modvsq, hue;
  int size, rank, irank, i,j, ix, iy;
  int nvel, nrgb;

  int tag=1;

  MPI_Status status;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  if (rank==0) printf("\n\nWriting data files ...\n");

  vel = (vecvel **) arraymalloc2d(m,n,sizeof(vecvel));
  rgb = (vecrgb **) arraymalloc2d(m,n,sizeof(vecrgb));

  //calculate velocities and hues

  double v1, v2;

  for (i=0;i<m;i++)
    {
      for (j=0;j<n;j++)
	{
	  vel[i][j][0] =  (psi[i+1][j+2]-psi[i+1][j])/2.0;
	  vel[i][j][1] = -(psi[i+2][j+1]-psi[i][j+1])/2.0;

	  v1 = vel[i][j][0];
	  v2=  vel[i][j][1];

	  modvsq = v1*v1 + v2*v2;

	  hue = pow(modvsq,0.4);

	  hue2rgb(hue,&(rgb[i][j][0]),&(rgb[i][j][1]),&(rgb[i][j][2]));
	}
    }

  //receive ddata

  if (rank == 0)
    {
      FILE *cfile, *vfile;

      cfile=fopen("colourmap.dat","w");
      vfile=fopen("velocity.dat","w");

      for (irank=0;irank<size;irank++)
	{
	  if (irank != 0)
	    {
	      MPI_Recv(&rgb[0][0][0],3*m*n,MPI_INT,irank,tag,comm,&status);
	      MPI_Recv(&vel[0][0][0],2*m*n,MPI_DOUBLE,irank,tag,comm,&status);
	    }

	  for (i=0;i<m;i++)
	      {
		ix = irank*m+i+1;

                for (j=0;j<n;j++)
		  {
		    iy = j+1;

                    fprintf(cfile,"%i %i %i %i %i\n", ix, iy,
			    rgb[i][j][0], rgb[i][j][1],rgb[i][j][2]);

                    if ((ix-1)%scale == (scale-1)/2 &&
			(iy-1)%scale == (scale-1)/2    )
		      {

                        fprintf(vfile,"%i %i %f %f\n",
				ix,iy,vel[i][j][0],vel[i][j][1]);
		      }
		  }
	      }
	  }

        fclose(vfile);
        fclose(cfile);
      }
    else
      {
	MPI_Ssend(&rgb[0][0][0],3*m*n,MPI_INT,0,tag,comm);
	MPI_Ssend(&vel[0][0][0],2*m*n,MPI_DOUBLE,0,tag,comm);
      }

  free(rgb);
  free(vel);

  if (rank ==0) printf("... done!\n");
}

void writeplotfile(int m, int n, int scale)
{
  FILE *gnuplot;

  gnuplot = fopen("cfd.plt","w");

  fprintf(gnuplot,"set size square\n");
  fprintf(gnuplot,"set key off\n");
  fprintf(gnuplot,"unset xtics\n");
  fprintf(gnuplot,"unset ytics\n");

  fprintf(gnuplot,"set xrange [%i:%i]\n",1-scale,m+scale);
  fprintf(gnuplot,"set yrange [%i:%i]\n",1-scale,n+scale);

  fprintf(gnuplot,"plot \"colourmap.dat\" w rgbimage, \"velocity.dat\" u 1:2:(%d*0.75*$3/sqrt($3**2+$4**2)):(%d*0.75*$4/sqrt($3**2+$4**2)) with vectors  lc rgb \"#7F7F7F\"",scale,scale);

  fclose(gnuplot);

  printf("\nWritten gnuplot script 'cfd.plt'\n");
}


void hue2rgb(double hue, int *r, int *g, int *b)
{
  int rgbmax = 255;

  *r = (int)(rgbmax*colfunc(hue-1.0));
  *g = (int)(rgbmax*colfunc(hue-0.5));
  *b = (int)(rgbmax*colfunc(hue    ));
}


double colfunc(double x)
{
  double absx;

  double x1=0.2;
  double x2=0.5;

  absx=fabs(x);

  if (absx > x2)
    {
      return 0.0;
    }
  else if (absx < x1)
    {
      return 1.0;
    }
  else
    {
      return 1.0-pow((absx-x1)/(x2-x1),2);
    }
}
