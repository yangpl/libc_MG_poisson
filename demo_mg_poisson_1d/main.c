/*
 * Demo for solving 1D possion equation using multigrid method
 * Author: Pengliang Yang, ypl.2100@gmail.com 
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cstd.h"

#define PI 3.14159265

typedef struct {
  int nx; //nx=2^L, total number of grid points should be nx+1
  double dx; //grid spacing
  double *u; // unknown
  double *f; // right hand side
  double *r; // r=f-Au
} grid1d;

int v0; //number of V cycles in F cycle
int v1, v2;// number of pre-/post- smoothing
char *smoother;//GS or RBGS

double dotprod(int n, double *x, double *y)
{
  int i;
  double s =0;

  for(i=0; i<n; i++) s+=x[i]*y[i];

  return s;
}

int maxlevel(int nx)//find L if nx=2^L
{
  int n=nx;
  int i=0;

  while(n!=1){
    n/=2;
    i++;
  }
  return i;
}

//Gauss-Seidel relaxation to update u, boundary points are not updated, homogeneous Dirichlet B.C.
void relax(int nx, double dx, double *f, double *u)
{
  int i;

  if(strcmp(smoother,"RBGS")==0){
    //Red-black/Odd-Even Gauss-Seidel, can be parallelized!
    for(i=1; i<=nx-1; i+=2){//red points, odd indices
      u[i] = 0.5* (u[i+1] + u[i-1] -dx*dx*f[i]);
    }
    for(i=2; i<=nx-1; i+=2){//black points, even indices
      u[i] = 0.5* (u[i+1] + u[i-1] -dx*dx*f[i]);
    }
  }else{//sequential/lexicographic Gauss-Seidel
	//This is order dependent, and cannot be parallelized!
    for(i=1; i<=nx-1; i++){
      u[i] = 0.5* (u[i+1] + u[i-1] -dx*dx*f[i]);
    }
  }
}

//compute residual
void residual(int nx, double dx, double *f, double *u, double *r)
{
  int i;

  for(i=1; i<=nx-1; i++){
    r[i] = f[i] - (u[i+1]-2.*u[i]+u[i-1])/(dx*dx); // r = f- Au
  }
  r[0] = r[nx] = 0.;//boundary is not updated
}

//restriction operator, full weighting
//nxf: number of points on fine grid
//nxc: number of points on coase grid
void restriction(int nxf, double *rf, int nxc, double *fc)
{
  int i;

  for(i=1; i<=nxc-1; i++){
    fc[i] = 0.25* (rf[2*i-1] +2.*rf[2*i] + rf[2*i+1]);
  }
  fc[0] = rf[0];
  fc[nxc] = rf[nxf];//nxf = 2*nxc
}

//prolongation/interpolation operator
//nxf: number of points on fine grid
//nxc: number of points on coase grid
void prolongation(int nxc,  double *rc, int nxf, double *pf)
{
  int i;
    
  for(i=0; i<nxc; i++){
    pf[2*i] = rc[i];
    pf[2*i+1] = 0.5*(rc[i] + rc[i+1]);
  }
  pf[nxf] = rc[nxc]; // nxf = 2*nxc
}


//V cycle: V(v1, v2) from fine grid (2^level*h) to coarse grid (2^{lmax-1}*h)
//level: [0, lmax-1]
//v1: number of pre-smoothing/relaxtion, 0, 1,2
//v2: number of post-smoothing/relaxtion
void v_cycle(int lmax, int level, grid1d *g, double *res)
{
  int i, l;

  if(level==0) {
    //compute residual r^h = f^h-A^h u^h
    residual(g[level].nx, g[level].dx, g[level].f, g[level].u, g[level].r);
    *res = sqrt(dotprod(g[level].nx, g[level].r, g[level].r));
  }


  //------------------------------------------------------------
  //1. from fine to coarse grid step by step
  for(l=level; l<lmax-1; l++){
    //pre-smoothing relaxation: v1 sweeps
    for(i=1; i<=v1; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);

    //compute residual r^h = f^h-A^h u^h
    residual(g[l].nx, g[l].dx,g[l].f, g[l].u, g[l].r);

    //restriction f^2h = I_h^2h r^h
    restriction(g[l].nx, g[l].r, g[l+1].nx, g[l+1].f);

    //initial guess for level 2h = 0
    memset(g[l+1].u, 0, (g[l+1].nx+1)*sizeof(double));	    
  }

  //-------------------------------------------------------------
  //2. direct solve at coarsest grid level, now l=lmax-1
  //pre-smoothing relaxation: v1 sweeps
  for(i=1; i<=v1; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);
  //dirct solve at coarsest grid A u=f
  if(g[l].nx==2) {
    //exactly solve the problem directly at coarsest grid (nx=2, nx+1=3)
    //all points are boundary points except for a single interior point.
    //(u[i-1]-2.*u[i] +u[i+1])/(dx*dx)= f[i]; i=1
    i=1;
    g[l].u[i] = 0.5*( g[l].u[i-1]-g[l].u[i+1] - g[l].dx*g[l].dx*g[l].f[i]);
  }
  //post-smoothing relaxation: v2 sweeps
  for(i=1; i<=v2; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);

  //-------------------------------------------------------------
  //3. from coarse to fine grid step by step
  for(l=lmax-2; l>=level; l--){
    //prolongation  (2h->h)
    prolongation(g[l+1].nx, g[l+1].u, g[l].nx, g[l].r);

    //correcton on level-h
    for(i=0; i<g[l].nx; i++) g[l].u[i] += g[l].r[i];

    //post-smoothing/relaxation: v2 sweeps
    for(i=1; i<=v2; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);
  }
}

void w_cycle(int lmax, int level, grid1d *g, double *res)
{
  int i, l;

  if(level==0) {
    //compute residual r^h = f^h-A^h u^h
    residual(g[level].nx, g[level].dx,g[level].f, g[level].u, g[level].r);
    *res = sqrt(dotprod(g[level].nx, g[level].r, g[level].r));
  }

  //-------------------------------------------------------------
  //1. from fine to coarse grid step by step
  for(l=level; l<=lmax-2; l++){
    //pre-smoothing relaxation: v1 sweeps
    for(i=1; i<=v1; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);

    //compute residual r^h = f^h-A^h u^h
    residual(g[l].nx, g[l].dx,g[l].f, g[l].u, g[l].r);

    //=====================================================
    //start a V cycle based on the new solution at level l
    v_cycle(lmax, l, g, res);
    //=====================================================

    //restriction f^2h = I_h^2h r^h
    restriction(g[l].nx, g[l].r, g[l+1].nx, g[l+1].f);

    //initial guess for level 2h = 0
    memset(g[l+1].u, 0, (g[l+1].nx+1)*sizeof(double));	    
  }

  //------------------------------------------------------------
  //2. direct solve at coarsest grid level, now l=lmax-1
  //pre-smoothing relaxation: v1 sweeps
  for(i=1; i<=v1; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);
  //dirct solve at coarsest grid A u=f
  if(g[l].nx==2) {
    //exactly solve the problem directly at coarsest grid (nx=2, nx+1=3)
    //all points are boundary points except for a single interior point.
    //(u[i-1]-2.*u[i] +u[i+1])/(dx*dx)= f[i]; i=1
    i=1;
    g[l].u[i] = 0.5*( g[l].u[i-1]-g[l].u[i+1] - g[l].dx*g[l].dx*g[l].f[i]);
  }
  //post-smoothing relaxation: v2 sweeps
  for(i=1; i<=v2; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);

  //-------------------------------------------------------------
  //3. from coarse grid to fine grid step by step
  for(l=lmax-2; l>=level; l--){
    //prolongation  (2h->h)
    prolongation(g[l+1].nx, g[l+1].u, g[l].nx, g[l].r);

    //correcton on level-h
    for(i=0; i<g[l].nx; i++) g[l].u[i] += g[l].r[i];

    //=====================================================
    //start a V cycle based on the new solution at level l
    v_cycle(lmax, l, g, res);
    //=====================================================
	
    //post-smoothing/relaxation: v2 sweeps
    for(i=1; i<=v2; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);
  }
}

//full multigrid V cycle (F cycle)
//F cycle: from fine grid (2^level*h) to coarse grid (2^{lmax-1}*h)
//level: [0, lmax-1]
void f_cycle(int lmax, int level, grid1d *g, double *res)
{
  int i, l;

  if(level==0) {
    //return residual at finest grid for computing convergence rate
    //compute residual r^h = f^h-A^h u^h
    residual(g[level].nx, g[level].dx,g[level].f, g[level].u, g[level].r);
    *res = sqrt(dotprod(g[level].nx, g[level].r, g[level].r));
  }

  //-------------------------------------------------------------
  //1. from fine to coarse grid step by step
  for(l=level; l<=lmax-2; l++){
    //pre-smoothing relaxation: v1 sweeps
    for(i=1; i<=v1; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);

    //compute residual r^h = f^h-A^h u^h
    residual(g[l].nx, g[l].dx,g[l].f, g[l].u, g[l].r);

    //restriction f^2h = I_h^2h r^h
    restriction(g[l].nx, g[l].r, g[l+1].nx, g[l+1].f);

    //initial guess for level 2h = 0
    memset(g[l+1].u, 0, (g[l+1].nx+1)*sizeof(double));	     
  }

  //------------------------------------------------------------
  //2. direct solve at coarsest grid level, now l=lmax-1
  //pre-smoothing relaxation: v1 sweeps
  for(i=1; i<=v1; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);
  //dirct solve at coarsest grid A u=f
  if(g[l].nx==2) {
    //exactly solve the problem directly at coarsest grid (nx=2, nx+1=3)
    //all points are boundary points except for a single interior point.
    //(u[i-1]-2.*u[i] +u[i+1])/(dx*dx)= f[i]; i=1
    i=1;
    g[l].u[i] = 0.5*( g[l].u[i-1]-g[l].u[i+1] - g[l].dx*g[l].dx*g[l].f[i]);
  }
  //post-smoothing relaxation: v2 sweeps
  for(i=1; i<=v2; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);

  //-------------------------------------------------------------
  //3. from coarse grid to fine grid step by step
  for(l=lmax-2; l>=level; l--){
    //prolongation  (2h->h)
    prolongation(g[l+1].nx, g[l+1].u, g[l].nx, g[l].r);

    //correcton on level-h
    for(i=0; i<g[l].nx; i++) g[l].u[i] += g[l].r[i];

    //=====================================================
    //start a V cycle based on the new solution at level l
    for(i=1; i<=v0; i++) v_cycle(lmax, l, g, res);//v0 times of V cycle
    //=====================================================
	
    //post-smoothing/relaxation: v2 sweeps
    for(i=1; i<=v2; i++) relax(g[l].nx, g[l].dx, g[l].f, g[l].u);
  }
}
static char *helper=" A multigrid solver for 1D possion equation \n\
    \n\
    v0: number of V cycle in full multigrid method \n\
    v1: number of pre-smoothing/relaxation sweeps \n\
    v2: number of post-smoothing/relaxation sweeps \n\
    nx: number of points in 1D problem=nx+1 including 2 boundary points \n\
    lmax: maximum number of levels to move from fine grid to coarse grid \n\
    tol: convergence tolerance, default=1e-6 \n\
    niter: number of MG iterations \n\
    mgcycle: MG cycle, can be V, W or F \n\
    smoother: GS (Gauss-Seidel) or RBGS (red black Gauss-Seidel)\n";
    

int main(int argc, char* argv[])
{
  int nx, lmax, niter;
  double x0, xL, dx, tol;
  char *mgcycle;
  double *x, *f, *u;

  int i, k, l;
  double res, res_old;
  grid1d *g;

  if(argc==1){ printf("%s",helper); exit(0);}
  initargs(argc,argv);

  if(!getparint("v0", &v0)) v0=1; //number of V-cycle in full multigrid
  if(!getparint("v1", &v1)) v1=1; //number of pre-smoothing/pre-relaxation sweeps
  if(!getparint("v2", &v2)) v2=1; //number of post-smoothing/post-relaxation sweeps
  if(!getparint("nx",&nx))  nx=128;
  if(!getparint("lmax", &lmax)) lmax=maxlevel(nx);
  if(!getpardouble("x0",&x0)) x0=0.;
  if(!getpardouble("xL",&xL)) xL=1.;
  if(!getpardouble("tol",&tol)) tol=1.e-6;
  if(!getparint("niter", &niter)) niter=100;
  if(!(getparstring("mgcycle", &mgcycle))) mgcycle="F";//mgcycle=V,W,F
  if(strcmp(mgcycle,"V")!=0 && strcmp(mgcycle,"W")!=0 && strcmp(mgcycle,"F")!=0) 
    err("mgcycle must be V,W or F!");
  if(!(getparstring("smoother", &smoother))) smoother="GS";//smoother=GS,RBGS
  if(strcmp(smoother,"GS")!=0 && strcmp(smoother,"RBGS")!=0) 
    err("smoother must be GS or RBGS!");
  //if(!sf_getbool("GS",smoother)) smoother=true;

  dx = (xL-x0)/(double)nx;

  warn("nx=%d", nx);
  warn("x0=%f", x0);
  warn("xL=%f", xL);
  warn("tol=%e", tol);
  warn("dx=%f", dx);
  warn("v0=%d", v0);
  warn("v1=%d", v1);
  warn("v2=%d", v2);
  warn("lmax=%d", lmax);

  x=alloc1double(nx+1);
  f=alloc1double(nx+1);
  u=alloc1double(nx+1);
  for(i=0; i<=nx; i++) {
    x[i] = x0 + i*dx;
    f[i] = 0.5*( sin(PI*x[i]) +sin(16.0*PI*x[i]) );//right hand side
  }
  memset(u, 0, (nx+1)*sizeof(double));//initialize unknowns to be 0s

  g = malloc(lmax*sizeof(grid1d));
  g[0].nx = nx;
  g[0].dx = dx;
  g[0].u = alloc1double(nx+1);
  g[0].f = alloc1double(nx+1);
  g[0].r = alloc1double(nx+1);
  memcpy(g[0].u, u, (nx+1)*sizeof(double));
  memcpy(g[0].f, f, (nx+1)*sizeof(double));


  FILE *fp=fopen("iterate.txt","w");
  for(k=0; k<niter; k++){
    for(l=1; l<lmax; l++){
      g[l].nx = g[l-1].nx/2;
      g[l].dx = g[l-1].dx*2;
      g[l].u = alloc1double(g[l].nx+1);
      g[l].f = alloc1double(g[l].nx+1);
      g[l].r = alloc1double(g[l].nx+1);
    }

    if(strcmp(mgcycle,"V")==0) 
      v_cycle(lmax, 0, g, &res);
    else if(strcmp(mgcycle,"W")==0) 
      w_cycle(lmax, 0, g, &res);
    else if(strcmp(mgcycle,"F")==0) 
      f_cycle(lmax, 0, g, &res);

    for(l=1; l<lmax; l++){
      free(g[l].u);
      free(g[l].f);
      free(g[l].r);
    }  

    if (k==0) {
      res_old = res;
      fprintf(fp, "#Iteration \t Misfit \t Convergence\n");
      printf("#Iteration \t Misfit \t Convergence\n");
    }
    fprintf(fp, "%d \t %e \t %e\n", k, res, res/res_old);
    printf("%d \t %e \t %e\n", k, res, res/res_old);
    res_old = res;
    if(res<=tol) break; 	//check convergence
  }
  fclose(fp);

  free(g[0].u);
  free(g[0].f);
  free(g[0].r);

  free(x);
  free(f);

  exit(0);
}

