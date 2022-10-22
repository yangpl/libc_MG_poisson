/* Demo for solving 2D possion equation using multigrid method
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"

int itermax;
int v1;//number of pre-smoothing
int v2;//number of post-smoothing
int lmax;//maximum number of levels for switching between fine and coarse grid
double tol;//tolerance for convergence
double rnorm;//norm of residual vector r

typedef struct{
  int nx;
  int ny;
  double dx;
  double dy;
  double **u;
  double **f;
} gmg_t;
gmg_t *gmg;


//compute dot product between two vectors s=<x,y>
double inner_product(int n, double *x, double *y)
{
  int i;
  double s;

  for(i=0, s=0; i<n; i++) s += x[i]*y[i];
  return s;
}

//relaxation/smoothing step for multigrid
void smoothing(gmg_t *gmg, int lev)
{
  int i, j;
  double a, tmp1, tmp2, _dx2, _dy2, **u, **f;

  u = gmg[lev].u;
  f = gmg[lev].f;
  _dx2 = 1./(gmg[lev].dx*gmg[lev].dx);
  _dy2 = 1./(gmg[lev].dy*gmg[lev].dy);
  a = -2.*_dx2 - 2.*_dy2;

  for(j=1; j<gmg[lev].ny; j++){
    for(i=1; i<gmg[lev].nx; i++){
      tmp1 = (u[j][i+1] + u[j][i-1])*_dx2;
      tmp2 = (u[j+1][i] + u[j-1][i])*_dy2;
      u[j][i] = (f[j][i] - tmp1 - tmp2)/a;
    }
  }
}

//compute residual r=f-Au
void residual(gmg_t *gmg, double **r, int lev)
{
  int i, j;
  double _dx2, _dy2, tmp1, tmp2, **u, **f;
  
  _dx2 = 1./(gmg[lev].dx*gmg[lev].dx);
  _dy2 = 1./(gmg[lev].dy*gmg[lev].dy);  
  u = gmg[lev].u;
  f = gmg[lev].f;
  for(j=1; j<gmg[lev].ny; j++){
    for(i=1; i<gmg[lev].nx; i++){
      tmp1 = (u[j][i+1] - 2.*u[j][i] + u[j][i-1])*_dx2;
      tmp2 = (u[j+1][i] - 2.*u[j][i] + u[j-1][i])*_dy2;
      r[j][i] = f[j][i] - tmp1 - tmp2;
    }
  }

  for(i=0; i<=gmg[lev].nx; i++){
    r[0][i] = 0.;
    r[gmg[lev].ny][i] = 0.;
  }
  for(j=0; j<=gmg[lev].ny; j++) {
    r[j][0] = 0.;
    r[j][gmg[lev].nx] = 0.;
  }
}

//interpolate u from (lev+1) to lev-th level
void prolongation(gmg_t *gmg, double **r, int lev)
{
  double **u;
  int i, j;

  u = gmg[lev+1].u;
  for(j=0; j<gmg[lev+1].ny; j++){
    for(i=0; i<gmg[lev+1].nx; i++){
      r[2*j][2*i] = u[j][i];
      r[2*j][2*i+1] = 0.5*(u[j][i] + u[j][i+1]);
      r[2*j+1][2*i] = 0.5*(u[j][i] + u[j+1][i]);
      r[2*j+1][2*i+1] = 0.25*(u[j][i] + u[j][i+1] + u[j+1][i] + u[j+1][i+1]);
    }
  }

  //handle boundaries
  for(i=0; i<gmg[lev+1].nx; i++) r[gmg[lev].ny][2*i] = u[gmg[lev+1].ny][i];
  for(j=0; j<gmg[lev+1].ny; j++) r[2*j][gmg[lev].nx] = u[j][gmg[lev+1].nx];
}

//restrict r from lev to (lev+1)-th lev: fine to coarse grid
void restriction(gmg_t *gmg, double **r, int lev)
{
  int i, j;
  double tmp1, tmp2, **f;

  //full weighting operator for restriction
  f = gmg[lev+1].f;
  for(j=1; j<gmg[lev+1].ny; j++){
    for(i=1; i<gmg[lev+1].nx; i++){
      tmp1 = r[2*j-1][2*i-1] + r[2*j-1][2*i+1] + r[2*j+1][2*i-1] + r[2*j+1][2*i+1];
      tmp2 = r[2*j][2*i-1] + r[2*j][2*i+1] + r[2*j-1][2*i] + r[2*j+1][2*i];
      f[j][i] = (tmp1 + 2*tmp2 + 4*r[2*j][2*i])/16.;
    }
  }

  //use exactly the same value from fine to coarse grid at boundaries
  for(i=0; i<=gmg[lev+1].nx; i++) {
    f[0][i] = r[0][i];
    f[gmg[lev+1].ny][i] = r[gmg[lev].ny][i];
  }
  for(j=0; j<=gmg[lev+1].ny; j++) {
    f[j][0] = r[j][0];
    f[j][gmg[lev+1].nx] = r[j][gmg[lev].nx];
  }
}

//multigrid V-cycle
void v_cycle(gmg_t *gmg, int lev)
{
  int i, j;
  double **r;
  
  for(i=0; i<v1; i++) smoothing(gmg, lev);//pre-smoothing of u based on u,f at lev-th level

  if(lev<lmax-1){
    r = alloc2double(gmg[lev].nx+1, gmg[lev].ny+1);//residual vector
    residual(gmg, r, lev);//residual r=f-Au at lev-th lev

    if(lev==0){//compute the norm of the residual vector at the beginning of each iteration
      rnorm = sqrt(inner_product((gmg[lev].nx+1)*(gmg[lev].ny+1), &r[0][0], &r[0][0]));
    }
    
    restriction(gmg, r, lev);//restrict r at lev-th lev to gmg[lev+1].f 

    memset(&gmg[lev+1].u[0][0], 0, (gmg[lev+1].nx+1)*(gmg[lev+1].ny+1)*sizeof(double));
    v_cycle(gmg, lev+1);// another v-cycle at (lev+1)-th level

    prolongation(gmg, r, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
    for(j=1; j<gmg[lev].ny; j++)
      for(i=1; i<gmg[lev].nx; i++)
	gmg[lev].u[j][i] += r[j][i];//correct u=u+r at interior without boundaries
    free2double(r);
  }
  //if lev==lmax-1, then nx=ny=2, grid size=3*3, only 1 point at the center is unknwn
  //direct solve is equivalent to smoothing at center point, one post-smoothing will do the joib
  
  for(i=0; i<v2; i++) smoothing(gmg, lev);//post-smoothing
}

//multigrid F-cycle
void f_cycle(gmg_t *gmg, int lev)
{
  int i, j;
  double **r;
    
  if(lev==lmax-1){//coarsest grid, direct solve or smoothing
    memset(&gmg[lev].u[0][0], 0, (gmg[lev].nx+1)*(gmg[lev].ny+1)*sizeof(double));
  }else{
    r = alloc2double(gmg[lev].nx+1, gmg[lev].ny+1);
    residual(gmg, r, lev);//residual r=f-Au at lev-th level
    restriction(gmg, r, lev);//restrict r at lev-th lev to f at (lev+1)-th level

    memset(&gmg[lev+1].u[0][0], 0, (gmg[lev+1].nx+1)*(gmg[lev+1].ny+1)*sizeof(double));
    f_cycle(gmg, lev+1);
    prolongation(gmg, r, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
    for(j=1; j<gmg[lev].ny; j++)
      for(i=1; i<gmg[lev].nx; i++)
	gmg[lev].u[j][i] += r[j][i];//correct u=u+r at interior without boundaries
    free2double(r);
  }
  v_cycle(gmg, lev);// another v-cycle at (lev+1)-th level
}

void gmg_init(int nx, int ny, double dx, double dy)
{
  int nx_;
  int i;

  if(!getparint("itermax", &itermax)) itermax = 10;/* maximum number of iterations */  
  if(!getparint("v1", &v1)) v1 = 1;/* number of pre-smoothing */
  if(!getparint("v2", &v2)) v2 = 1;/* number of post-smoothing */
  if(!getpardouble("tol", &tol)) tol = 1e-6;/* stopping criteria */
  if(!getparint("lmax", &lmax)) {
    lmax = 0;
    nx_ = nx;
    while(nx_>1) {nx_=nx_>>1; lmax++;}
  }
  
  gmg = (gmg_t *)malloc(lmax*sizeof(gmg_t));
  for(i=0; i<lmax; i++){
    gmg[i].nx = nx/(1<<i);//nx/2^i
    gmg[i].ny = ny/(1<<i);//ny/2^i
    gmg[i].dx = dx*(1<<i);//dx*2^i
    gmg[i].dy = dy*(1<<i);//dy*2^i    
    gmg[i].u = (double**) alloc2double(gmg[i].nx+1, gmg[i].ny+1);
    gmg[i].f = (double**) alloc2double(gmg[i].nx+1, gmg[i].ny+1);
  }
  
}

void gmg_close()
{
  int i;
  
  for(i=0; i<lmax; i++){
    free2double(gmg[i].u);
    free2double(gmg[i].f);
  }
  free(gmg);
}

void gmg_apply(int n, double *b, double *x)
{
  int iter;

  memcpy(&gmg[0].f[0][0], b, n*sizeof(double));
  memset(&gmg[0].u[0][0], 0, n*sizeof(double));
  for(iter=0; iter<itermax; iter++) f_cycle(gmg, 0);
  memcpy(x, &gmg[0].u[0][0], n*sizeof(double));
}

