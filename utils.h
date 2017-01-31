#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#define NR_END 1
#define I3D(Nx,Ny,Nz,i,j,k,n) ((i)+(Nx)*(j)+(Nx)*(Ny)*(k)+(Nx)*(Ny)*(Nz)*(n))

double *dvector (long N)
{
	double *v;
	v = (double *)malloc(sizeof(double)*N);
	return v;
}

void zero_matrix(double *a, int Nx, int Ny, int Nz, int Ng)
{
	int i; int j; int k; int n;
      for (n=0; n < Ng; n++) {
         for (k = 0; k < Nz; k++) {
 	    for (j = 0; j < Ny; j++) {
	       for (i = 0; i < Nx; i++) {
			a[I3D(Nx, Ny, Nz, i, j, k, n)] = 0.0;
	       }
	    }
         }
      }
}

void print_mat(FILE *fptr,
	double *a, int Nx, int Ny, int Nz, int Ng)
{
	int i; int j; int k; int n;
       for (n=0; n < Ng; n++) {
          for (k = 0; k < Nz; k++) {
	     for (j = 0; j < Ny; j++) {
		 for (i = 0; i < Nx; i++) {
			fprintf(fptr, "%.16f ", a[I3D(Nx, Ny, Nz, i, j, k, n)]);
		 }
	     }
          }
	    fprintf(fptr, "\n");
      }
            fclose(fptr);
}

void initialize(double *a, int Nx, int Ny, int Nz, int Ng)
{
	int i; int j; int k; int n;
      for (n=0; n < Ng; n++) {
         for (k = 1; k < Nz-1; k++) {
	     for (j = 1; j < Ny-1; j++) {
		 for (i = 1; i < Nx-1; i++) {
			a[I3D(Nx, Ny, Nz, i, j, k, n)] = ((double) rand() / (RAND_MAX));
		}
	     }
         }
      }  
}
int iDivUp(int a, int b){ 
	return (((a) + (b) - 1)/(b));
}
