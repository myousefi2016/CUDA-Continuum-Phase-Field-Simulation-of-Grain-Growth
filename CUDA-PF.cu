#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <cuda.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "utils.h"
#include "dev_matrix.h"

#define I3D(Nx,Ny,Nz,i,j,k,n) ((i)+(Nx)*(j)+(Nx)*(Ny)*(k)+(Nx)*(Ny)*(Nz)*(n))

#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#define BLOCK_SIZE_Z 16

__global__ void pf3d_gpu(int Nx, int Ny, int Nz, int Ng, double L, double alpha, double beta, double gamma, double kappa, double ht, double hx, double hy, double hz, double *in, double *out, double *pf) 
{
	int i, j, k, n, m, P, Pm, W, E, S, N, U, B, temp;
	double lxx, lyy, lzz, sum;
        sum=0.0;
	// find i and j indices of this thread
	i = blockIdx.x*(BLOCK_SIZE_X) + threadIdx.x;
	j = blockIdx.y*(BLOCK_SIZE_Y) + threadIdx.y;
        k = blockIdx.z*(BLOCK_SIZE_Z) + threadIdx.z;

	// find indices into linear memory  
        for (n=0; n < Ng; n++) {
        P = I3D(Nx, Ny, Nz, i, j, k, n);
	W = I3D(Nx, Ny, Nz, i-1, j, k, n); E = I3D(Nx, Ny, Nz, i+1, j, k, n);
	S = I3D(Nx, Ny, Nz, i, j-1, k, n); N = I3D(Nx, Ny, Nz, i, j+1, k, n);
        B = I3D(Nx, Ny, Nz, i, j, k-1, n); U = I3D(Nx, Ny, Nz, i, j, k+1, n);
	// check that thread is within domain (not on boundary or outside domain)
	if (i > 0 && i < Nx-1 && j > 0 && j < Ny-1 && k>0 && k<Nz-1) {
		lxx = (in[W] - 2.0*in[P] + in[E])/pow(hx,2);
		lyy = (in[N] - 2.0*in[P] + in[S])/pow(hy,2);
                lzz = (in[U] - 2.0*in[P] + in[B])/pow(hz,2);

              for (m=0; m < Ng && m!=n; m++) {
                Pm = I3D(Nx, Ny, Nz, i, j, k, m);
                sum=sum+pow(in[Pm],2);
              }

        out[P] = in[P]+ht*(alpha*L*in[P]-beta*L*pow(in[P],3)-2*gamma*L*in[P]*sum+kappa*L*(lxx+lyy+lzz));
	}
        if (i==0) {
        out[P] = in[E];
        }
        if (j==0) {
        out[P] = in[N];
        }
        if (k==0) {
        out[P] = in[U];
        }
        if (i==Nx-1) {
        temp = I3D(Nx, Ny, Nz, 0, j, k, n);
        out[P] = in[temp];
        }
        if (j==Ny-1) {
        temp = I3D(Nx, Ny, Nz, i, 0, k, n);
        out[P] = in[temp];
        }
        if (k==Nz-1) {
        temp = I3D(Nx, Ny, Nz, i, j, 0, n);
        out[P] = in[temp];
        }
        pf[P] = pf[P] + pow(out[P],2);
        }
}

int main()
{
	int Nx, Ny, Nz, Nt, Ng;
        double *uh_old, *uh_new, *pfh, *tmp_h;
        int iter;
        double L, kappa, alpha, beta, gamma;
        double hx; double hy; double hz; double ht;
	dim3 numBlocks, threadsPerBlock;

	Nx = 32;
	Ny = 32;
        Nz = 32;
        Nt = 10000;
        Ng = 5;
        hx = 2.0;
        hy = 2.0;
        hz = 2.0;
        ht = 0.25;
        L = 1.0;
        kappa = 2.0;
        alpha = 1.0;
        beta = 1.0;
        gamma = 1.0;
	uh_old = dvector(Nx*Ny*Nz*Ng); uh_new = dvector(Nx*Ny*Nz*Ng); pfh = dvector(Nx*Ny*Nz);

	zero_matrix(uh_old, Nx, Ny, Nz, Ng);
	zero_matrix(uh_new, Nx, Ny, Nz, Ng);
        zero_matrix(pfh, Nx, Ny, Nz, 1);

	// initial
	initialize(uh_old, Nx, Ny, Nz, Ng);

        dev_matrix<double> ud_old(Nx, Ny, Nz, Ng); ud_old.set(uh_old, Nx, Ny, Nz, Ng);
	dev_matrix<double> ud_new(Nx, Ny, Nz, Ng); ud_new.set(uh_new, Nx, Ny, Nz, Ng);
        dev_matrix<double> tmp_d(Nx, Ny, Nz, Ng);
        dev_matrix<double> pfd(Nx, Ny, Nz, 1);
        
	numBlocks = dim3(iDivUp(Nx,BLOCK_SIZE_X), iDivUp(Ny,BLOCK_SIZE_Y), iDivUp(Nz,BLOCK_SIZE_Z));
	threadsPerBlock = dim3(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);

        vtkSmartPointer<vtkImageData> imageData =
          vtkSmartPointer<vtkImageData>::New();

        imageData->SetDimensions(Nx, Ny, Nz);

        vtkSmartPointer<vtkDoubleArray> phase =
          vtkSmartPointer<vtkDoubleArray>::New();

        phase->SetNumberOfComponents(1);
        phase->SetNumberOfTuples(Nx * Ny * Nz);

	for (iter = 0; iter < Nt; iter++) {
		pf3d_gpu<<<numBlocks, threadsPerBlock>>>(Nx, Ny, Nz, Ng, L, alpha, beta, gamma, kappa, ht, hx, hy, hz, ud_old.getData(), ud_new.getData(), pfd.getData());
        tmp_d = ud_new;
        ud_new = ud_old;
        ud_old = tmp_d;
        tmp_h = uh_new;
        uh_new = uh_old;
        uh_old = tmp_h;
        char myfile[16];
        sprintf(myfile, "myfile_%d.vti", iter);
        pfd.get(pfh, Nx, Ny, Nz, 1);
        for (i=0; i < Nx; i++) {
        for (j=0; j < Ny; j++) {
        for (k=0; k < Nz; k++) {
        P = I3D(Nx, Ny, Nz, i, j, k, 0);
        phase->SetValue(P, pfh[P]);
        }
        }
        }
        imageData->GetPointData()->AddArray(phase);
        phase->SetName("Phase Field");
        vtkSmartPointer<vtkXMLImageDataWriter> writer =
          vtkSmartPointer<vtkXMLImageDataWriter>::New();

        writer->SetFileName(myfile);
    #if VTK_MAJOR_VERSION <= 5
        writer->SetInputConnection(imageData->GetProducerPort());
    #else
        writer->SetInputData(imageData);
    #endif
        writer->Write();
	} 
	cudaThreadSynchronize();
}


