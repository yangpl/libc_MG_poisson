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
  int nx, ny, nz, niter, n, verb;
  double dx, dy, dz, tol;
  int i, j, k;
  double *x, *b;
  
  initargs(argc,argv);

  if(!getparint("verb", &verb)) verb = 1;/* verbosity */
  if(!getparint("nx", &nx)) nx = 128;/* dimension in x */
  if(!getparint("ny", &ny)) ny = 128;/* dimension in y */
  if(!getparint("nz", &nz)) nz = 128;/* dimension in z */
  if(!getparint("niter", &niter)) niter = 10;
  if(!getpardouble("tol", &tol)) tol = 1e-6;
  
  dx = 1./nx;
  dy = 1./ny;
  dz = 1./nz;
  n = (nx+1)*(ny+1)*(nz+1);

  b = alloc1double(n);
  x = alloc1double(n);
  memset(x, 0, n*sizeof(double));//x=0 as initialization
  memset(b, 0, n*sizeof(double));//x=0 as initialization
  i = nx/2;
  j = ny/2;
  k = nz/2;
  b[i + (nx+1)*(j + (ny+1)*k)] = 1;

  gmg_init(nx, ny, nz, dx, dy, dz);//preconditioner initialization
  gmg_apply(n, b, x);//direct solve Ax=b by multigrid
  gmg_close();

  free1double(b);
  free1double(x);


  return 0;
}
