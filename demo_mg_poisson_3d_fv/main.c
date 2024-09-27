/* Demo for conjugate gradient method, BiCGStab and preconditioned BiCGStab
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "sparse.h"
#include "solver.h"

typedef void (*op_t)(int, double*, double*); //define type for linear operator

void gmg_init(int nx, int ny, int nz, double dx, double dy, double dz);
void gmg_apply(int n, double *b, double *x);//solve Ax=b by multigrid
void gmg_close();

void amg_init(csr_t A);
void amg_apply(int n, double *b, double *x);//solve Ax=b by multigrid
void amg_close();

void Acoo_init(int nx, int ny, double dx, double dy);
void Acoo_apply(int n, double *x, double *y);
void Acoo_close();

void Acsr_init(csr_t *A, int nx, int ny, int nz, double dx, double dy, double dz);
void Acsr_apply(int n, double *x, double *y);
void Acsr_close();

double dotprod(int n, double *a, double *b);
void solve_cg(int n, double *x, double *b, op_t Aop, int niter, double tol, int verb);
void solve_pcg(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol, int verb);
void solve_pbicgstab(int n, double *x, double *b, op_t Aop, op_t invKop, int niter, double tol, int verb);
void solve_gmres_rightpreco(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol, int m, int verb);
void solve_cgnr(int n, double *x, double *b, op_t Aop, op_t Atop, int niter, double tol, int verb);


void Id_apply(int n, double *b, double *x)
{
  memcpy(x, b, n*sizeof(double));
}

  
int main(int argc, char *argv[])
{
  int nx, ny, nz, method, niter, n, nrestart, mgopt, verb;
  double dx, dy, dz, x2, y2, z2, tol;
  int i, j, k, kk;
  double *x, *xt, *b;
  csr_t A;
  FILE *fp;
  
  initargs(argc,argv);

  if(!getparint("verb", &verb)) verb = 1;/* verbosity */
  if(!getparint("nx", &nx)) nx = 128;/* dimension in x */
  if(!getparint("ny", &ny)) ny = 128;/* dimension in y */
  if(!getparint("nz", &nz)) nz = 128;/* dimension in z */
  if(!getparint("method", &method)) method = 0;
  if(!getparint("niter", &niter)) niter = 10;
  if(!getparint("mgopt", &mgopt)) mgopt = 1;
  if(!getpardouble("tol", &tol)) tol = 1e-6;
  if(!getparint("nrestart", &nrestart)) nrestart = 10;/* GMRES restarts every nrestart steps  */
  
  dx = 1./nx;
  dy = 1./ny;
  dz = 1./nz;
  n = (nx+1)*(ny+1)*(nz+1);

  printf("method=%d\n", method);
  printf("niter=%d\n", niter);
  printf("mgopt=%d\n", mgopt);
  
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
  Acsr_init(&A, nx, ny, nz, dx, dy, dz);
  Acsr_apply(n, xt, b);//b=A*xt

  gmg_init(nx, ny, nz, dx, dy, dz);//preconditioner initialization
  gmg_apply(n, b, x);//direct solve Ax=b by multigrid
  gmg_close();

  
  /* output true signal and reconstructed one */
  fp = fopen("result.txt", "w");
  fprintf(fp, "x_true \t x_rec \n");
  for(i=0; i<n; i++) fprintf(fp, "%e \t %e\n", xt[i], x[i]);
  fclose(fp);

  Acsr_close();
  free1double(b);
  free1double(x);
  free1double(xt);

  return 0;
}
