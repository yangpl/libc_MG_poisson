/* Demo for the solution of Poisson equation using direct solver - MUMPS
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "sparse.h"
#include "mpi.h"
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654


void Acoo_init(coo_t *A, int nx, int ny, double dx, double dy);
void Acoo_apply(int n, double *x, double *y);
void Acoo_close();
void Acoo_mumps(coo_t *A_);

int main(int argc, char *argv[])
{
  int ierr, iproc;
  int nx, ny, method, niter, n, nrestart, mgopt, verb;
  double dx, dy, x2, y2, tol;
  int i, j, k;
  double *x, *xt, *b;
  coo_t *A;
  FILE *fp;

  initargs(argc,argv);

  if(!getparint("verb", &verb)) verb = 1;/* verbosity */
  if(!getparint("nx", &nx)) nx = 128;/* dimension in x */
  if(!getparint("ny", &ny)) ny = 128;/* dimension in y */
  if(!getparint("method", &method)) method = 0;
  if(!getparint("niter", &niter)) niter = 10;
  if(!getparint("mgopt", &mgopt)) mgopt = 0;
  if(!getpardouble("tol", &tol)) tol = 1e-6;
  if(!getparint("nrestart", &nrestart)) nrestart = 10;/* GMRES restarts every nrestart steps  */

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
      xt[k] = (x2-x2*x2)*(y2-y2*y2);
    }
  }
  memset(x, 0, n*sizeof(double));//x=0 as initialization
  A = malloc(sizeof(coo_t));
  Acoo_init(A, nx, ny, dx, dy);
  Acoo_apply(n, xt, b);//b=A*xt




  int error = 0;  
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  DMUMPS_STRUC_C id;

  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
  id.comm_fortran = USE_COMM_WORLD;
  id.par = 1;
  id.sym = 0;
  id.job = JOB_INIT;
  dmumps_c(&id);

  printf("nrow=%d\n", A->nrow);
  /* Define the problem on the host */
  if (iproc == 0) {
    id.n = A->nrow;
    id.nnz = (MUMPS_INT8) A->nnz;
    id.irn = A->row;
    id.jcn = A->col;
    //the row and column indices called by MUMPS are starting from 1 in Fortran
    for(k=0; k<id.nnz; k++){
      id.irn[k]++;
      id.jcn[k]++;
    }
    id.a = A->val;
    id.rhs = b;
  } 
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
  /* No outputs */
  id.ICNTL(1) = -1;
  id.ICNTL(2) = -1;
  id.ICNTL(3) = -1;
  id.ICNTL(4) = 0;

  /* Call the MUMPS package (analyse, factorization and solve). */
  id.job=6;
  dmumps_c(&id);
  if (id.infog[0]<0) {
    printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n", iproc, id.infog[0], id.infog[1]);
    error = 1;
  }

  /* Terminate instance. */
  id.job = JOB_END;
  dmumps_c(&id);
  if (iproc == 0) {
    if (!error) {
      /* output true signal and reconstructed one */
      fp = fopen("result.txt", "w");
      fprintf(fp, "x_true \t x_rec \n");
      for(i=0; i<n; i++) fprintf(fp, "%e \t %e\n", xt[i], b[i]);
      fclose(fp);
    } else {
      printf("An error has occured, please check error code returned by MUMPS.\n");
    }
  }
  

  Acoo_close();
  free1double(b);
  free1double(x);
  free1double(xt);

  ierr = MPI_Finalize();
  return 0;
}
