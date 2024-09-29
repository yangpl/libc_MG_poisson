/* Demo for conjugate gradient method, BiCGStab and preconditioned BiCGStab
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"

void gmg_init(int nx, int ny, int nz, double dx, double dy, double dz);
void gmg_apply(int n, double *b, double *x);//solve Ax=b by multigrid
void gmg_close();

  
int main(int argc, char *argv[])
{
  int nx, ny, nz, n, verb;
  double dx, dy, dz, x2, y2, z2;
  int i, j, k, kk;
  double *x, *xt, *b;
  
  initargs(argc,argv);

  if(!getparint("verb", &verb)) verb = 1;/* verbosity */
  if(!getparint("nx", &nx)) nx = 128;/* dimension in x */
  if(!getparint("ny", &ny)) ny = 128;/* dimension in y */
  if(!getparint("nz", &nz)) nz = 128;/* dimension in z */
  if(!getpardouble("dx", &dx)) dx = 1./nx;
  if(!getpardouble("dy", &dy)) dy = 1./ny;
  if(!getpardouble("dz", &dz)) dz = 1./nz;
  printf("problem size [nx, ny, nz]=[%d, %d, %d]\n", nx, ny, nz);
  printf("grid spacing [dx, dy, dz]=[%g, %g, %g]\n", dx, dy, dz);
  
  n = (nx+1)*(ny+1)*(nz+1);

  b = alloc1double(n);
  x = alloc1double(n);
  xt = alloc1double(n);

  for(k=0; k<=nz; k++){
    z2 = k*dz;
    z2 = z2*z2;
    for(j=0; j<=ny; j++){
      y2 = j*dy;
      y2 = y2*y2;
      for(i=0; i<=nx; i++){
	x2 = i*dx;
	x2 = x2*x2;
	kk = i + (nx+1)*(j + (ny+1)*k);
	xt[kk] = (x2-x2*x2)*(y2-y2*y2)*(z2-z2*z2);
      }
    }
  }
  memset(x, 0, n*sizeof(double));//x=0 as initialization
  memset(b, 0, n*sizeof(double));//x=0 as initialization
  b[n/2] = 1.;
  
  gmg_init(nx, ny, nz, dx, dy, dz);//preconditioner initialization
  gmg_apply(n, b, x);//direct solve Ax=b by multigrid
  gmg_close();

  
  free1double(b);
  free1double(x);
  free1double(xt);

  return 0;
}
