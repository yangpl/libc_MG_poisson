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

void gmg_init(int nx, int ny, double dx, double dy);
void gmg_apply(int n, double *b, double *x);//solve Ax=b by multigrid
void gmg_close();

void amg_init(csr_t *A);
void amg_apply(int n, double *b, double *x);//solve Ax=b by multigrid
void amg_close();


void Acoo_init(int nx, int ny, double dx, double dy);
void Acoo_apply(int n, double *x, double *y);
void Acoo_close();

void Acsr_init(csr_t *A, int nx, int ny, double dx, double dy);
void Acsr_apply(int n, double *x, double *y);
void Acsr_close();

double dotprod(int n, double *a, double *b);
//linear solver using conjugate gradient method
void solve_cg(int n, double *x, double *b, op_t Aop, int niter, double tol);
//linear solver using conjugate gradient method
void solve_pcg(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol);
//linear solver using BiCGStab
void solve_bicgstab(int n, double *x, double *b, op_t Aop, int niter, double tol);
//linear solver using preconditioned BiCGStab
void solve_pbicgstab(int n, double *x, double *b, op_t Aop, op_t invKop, int niter, double tol);
//preconditioned BiCGStab2
void solve_pbicgstab2(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol);
//GMRES without preconditioning
void solve_gmres(int n, double *x, double *b, op_t Aop, int niter, double tol, int m);
//GMRES with left preconditioning
void solve_gmres_leftpreco(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol, int m);
//GMRES with right preconditioning
void solve_gmres_rightpreco(int n, double *x, double *b, op_t Aop, op_t invMop, int niter, double tol, int m);
//CGNR method
void solve_cgnr(int n, double *x, double *b, op_t Aop, op_t Atop, int niter, double tol);


int main(int argc, char *argv[])
{
  int nx, ny, method, niter, n, mgopt;
  double dx, dy, x2, y2;
  int i, j, k;
  double *x, *xt, *b;
  csr_t A;
  FILE *fp;
  
  initargs(argc,argv);

  if(!getparint("nx", &nx)) nx = 128;/* dimension in x */
  if(!getparint("ny", &ny)) ny = 128;/* dimension in y */
  if(!getparint("method", &method)) method = 0;
  if(!getparint("niter", &niter)) niter = 20;
  if(!getparint("mgopt", &mgopt)) mgopt = 1;//0=GMG; 1=AMG
  dx = 1./nx;
  dy = 1./ny;
  n = (nx+1)*(ny+1);

  b = alloc1double(n);
  x = alloc1double(n);
  xt = alloc1double(n);
  
  for(j=0; j<=ny; j++){
    y2 = j*dy;
    y2 = y2*y2;
    for(i=0; i<=nx; i++){
      x2 = i*dx;
      x2 = x2*x2;
      k = i + (nx+1)*j;
      //b[k] = 2.*((1.-6.*x2)*y2*(1-y2) + (1.-6.*y2)*x2*(1.-x2));
      xt[k] = (x2-x2*x2)*(y2*y2-y2);//initialize right hand side
    }
  }
  memset(x, 0, n*sizeof(double));//x=0 as initialization
  Acsr_init(&A, nx, ny, dx, dy);
  Acsr_apply(n, xt, b);//b=A*xt

  printf("method=%d\n", method);
  printf("niter=%d\n", niter);

  if(mgopt==0){
    printf("----- solved by GMG-------\n");
    gmg_init(nx, ny, dx, dy);//preconditioner initialization
    if(method==0)
      gmg_apply(n, b, x);//direct solve Ax=b by multigrid
    else if(method==1)
      solve_pcg(n, x, b, Acsr_apply, gmg_apply, niter, 1e-6);
    else if(method==2)
      solve_pbicgstab(n, x, b, Acsr_apply, gmg_apply, niter, 1e-6);
    else if(method==3)//problematic now using preconditioner
      solve_gmres_rightpreco(n, x, b, Acsr_apply, gmg_apply, niter, 1e-6, 10);
    gmg_close();
  }else{
    printf("----- solved by AMG-------\n");
    amg_init(&A);//preconditioner initialization
    if(method==0)
      amg_apply(n, b, x);//direct solve Ax=b by multigrid
    else if(method==1)
      solve_pcg(n, x, b, Acsr_apply, amg_apply, niter, 1e-6);
    else if(method==2)
      solve_pbicgstab(n, x, b, Acsr_apply, amg_apply, niter, 1e-6);
    else if(method==3)//problematic now using preconditioner
      solve_gmres_rightpreco(n, x, b, Acsr_apply, amg_apply, niter, 1e-6, 10);
    amg_close();
  }
  
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
