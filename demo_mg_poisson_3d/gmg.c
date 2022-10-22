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
  int nz;
  double dx;
  double dy;
  double dz;
  double ***u;
  double ***f;
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
  int i, j, k;
  double a, tmp1, tmp2, tmp3, _dx2, _dy2, _dz2;
  double ***u, ***f;

  u = gmg[lev].u;
  f = gmg[lev].f;
  _dx2 = 1./(gmg[lev].dx*gmg[lev].dx);
  _dy2 = 1./(gmg[lev].dy*gmg[lev].dy);
  _dz2 = 1./(gmg[lev].dz*gmg[lev].dz);
  a = 2.*_dx2 + 2.*_dy2 + 2.*_dz2;
  for(k=1; k<gmg[lev].nz; k++){
    for(j=1; j<gmg[lev].ny; j++){
      for(i=1; i<gmg[lev].nx; i++){
	tmp1 = (u[k][j][i+1] + u[k][j][i-1])*_dx2;
	tmp2 = (u[k][j+1][i] + u[k][j-1][i])*_dy2;
	tmp3 = (u[k+1][j][i] + u[k-1][j][i])*_dz2;
	u[k][j][i] = (f[k][j][i] + tmp1 + tmp2 + tmp3)/a;
      }
    }
  }

}

//compute residual r=f-Au
void residual(gmg_t *gmg, double ***r, int lev)
{
  int i, j, k;
  double _dx2, _dy2, _dz2, tmp1, tmp2, tmp3;
  double ***u, ***f;
  
  _dx2 = 1./(gmg[lev].dx*gmg[lev].dx);
  _dy2 = 1./(gmg[lev].dy*gmg[lev].dy);
  _dz2 = 1./(gmg[lev].dz*gmg[lev].dz);  
  u = gmg[lev].u;
  f = gmg[lev].f;
  for(k=1; k<gmg[lev].nz; k++){
    for(j=1; j<gmg[lev].ny; j++){
      for(i=1; i<gmg[lev].nx; i++){
	tmp1 = (u[k][j][i+1] - 2.*u[k][j][i] + u[k][j][i-1])*_dx2;
	tmp2 = (u[k][j+1][i] - 2.*u[k][j][i] + u[k][j-1][i])*_dy2;
	tmp3 = (u[k+1][j][i] - 2.*u[k][j][i] + u[k-1][j][i])*_dz2;
	r[k][j][i] = f[k][j][i] + tmp1 + tmp2 + tmp3;
      }
    }
  }

  for(k=0; k<=gmg[lev].nz; k++){
    for(j=0; j<=gmg[lev].ny; j++) {
      r[k][j][0] = 0.;
      r[k][j][gmg[lev].nx] = 0.;
    }
  }

  for(k=0; k<=gmg[lev].nz; k++){
    for(i=0; i<=gmg[lev].nx; i++){
      r[k][0][i] = 0.;
      r[k][gmg[lev].ny][i] = 0.;
    }
  }

  for(j=0; j<=gmg[lev].ny; j++){
    for(i=0; i<=gmg[lev].nx; i++){
      r[0][j][i] = 0.;
      r[gmg[lev].nz][j][i] = 0.;
    }
  }
  
}

//interpolate u from (lev+1) to lev-th level
void prolongation(gmg_t *gmg, double ***r, int lev)
{
  double ***u;
  int i, j, k;

  u = gmg[lev+1].u;
  for(k=0; k<gmg[lev+1].nz; k++){
    for(j=0; j<gmg[lev+1].ny; j++){
      for(i=0; i<gmg[lev+1].nx; i++){
	r[2*k][2*j][2*i] = u[k][j][i];
	r[2*k][2*j][2*i+1] = 0.5*(u[k][j][i] + u[k][j][i+1]);
	r[2*k][2*j+1][2*i] = 0.5*(u[k][j][i] + u[k][j+1][i]);
	//r[2*k][2*j+1][2*i+1] = 0.5*(0.5*(u[k][j][i] + u[k][j][i+1]) + 0.5*(u[k][j+1][i] + u[k][j+1][i+1]));
	r[2*k][2*j+1][2*i+1] = 0.25*(u[k][j][i] + u[k][j][i+1] + u[k][j+1][i] + u[k][j+1][i+1]);
	r[2*k+1][2*j][2*i] = 0.5*(u[k][j][i] + u[k+1][j][i]);
	//r[2*k+1][2*j][2*i+1] = 0.5*(0.5*(u[k][j][i] + u[k][j][i+1]) + 0.5*(u[k+1][j][i] + u[k+1][j][i+1]));
	r[2*k+1][2*j][2*i+1] = 0.25*(u[k][j][i] + u[k][j][i+1] + u[k+1][j][i] + u[k+1][j][i+1]);
	//r[2*k+1][2*j+1][2*i] = 0.5*(0.5*(u[k][j][i] + u[k][j+1][i]) + 0.5*(u[k+1][j][i] + u[k+1][j+1][i]));
	r[2*k+1][2*j+1][2*i] = 0.25*(u[k][j][i] + u[k][j+1][i] + u[k+1][j][i] + u[k+1][j+1][i]);
	r[2*k+1][2*j+1][2*i+1] = 0.125*(  u[k][j][i] + u[k][j][i+1] + u[k][j+1][i] + u[k][j+1][i+1]
					+ u[k+1][j][i] + u[k+1][j][i+1] + u[k+1][j+1][i] + u[k+1][j+1][i+1]);
      }
    }
  }

  //handle boundaries
  for(k=0; k<=gmg[lev+1].nz; k++)  
    for(j=0; j<=gmg[lev+1].ny; j++)
      r[2*k][2*j][gmg[lev].nx] = u[k][j][gmg[lev+1].nx];
  
  for(k=0; k<=gmg[lev+1].nz; k++)
    for(i=0; i<=gmg[lev+1].nx; i++)
      r[2*k][gmg[lev].ny][2*i] = u[k][gmg[lev+1].ny][i];

  for(j=0; j<=gmg[lev+1].ny; j++)
    for(i=0; i<=gmg[lev+1].nx; i++)
      r[gmg[lev].nz][2*j][2*i] = u[gmg[lev+1].nz][j][i];
  
}

//restrict r from lev to (lev+1)-th lev: fine to coarse grid
void restriction(gmg_t *gmg, double ***r, int lev)
{
  int i, j, k;
  double tmp1, tmp2, tmp3, ***f;

  //full weighting operator for restriction
  f = gmg[lev+1].f;
  for(k=1; k<gmg[lev+1].nz; k++){
    for(j=1; j<gmg[lev+1].ny; j++){
      for(i=1; i<gmg[lev+1].nx; i++){
	tmp1 = r[2*k-1][2*j-1][2*i-1] + r[2*k-1][2*j-1][2*i+1] 
	  + r[2*k-1][2*j+1][2*i-1] + r[2*k-1][2*j+1][2*i+1]
	  + r[2*k+1][2*j-1][2*i-1] + r[2*k+1][2*j-1][2*i+1] 
	  + r[2*k+1][2*j+1][2*i-1] + r[2*k+1][2*j+1][2*i+1];
	//tmp1 = r[2*j-1][2*i-1] + r[2*j-1][2*i+1] + r[2*j+1][2*i-1] + r[2*j+1][2*i+1];
	tmp2 = r[2*k-1][2*j][2*i-1] + r[2*k-1][2*j][2*i+1]
	  + r[2*k-1][2*j-1][2*i] + r[2*k-1][2*j+1][2*i]
	  + r[2*k+1][2*j][2*i-1] + r[2*k+1][2*j][2*i+1]
	  + r[2*k+1][2*j-1][2*i] + r[2*k+1][2*j+1][2*i]
	  + r[2*k][2*j-1][2*i-1] + r[2*k][2*j+1][2*i-1]
	  + r[2*k][2*j-1][2*i+1] + r[2*k][2*j+1][2*i+1];
	//tmp2 = r[2*j][2*i-1] + r[2*j][2*i+1] + r[2*j-1][2*i] + r[2*j+1][2*i];
	tmp3 = r[2*k][2*j][2*i-1] + r[2*k][2*j][2*i+1]
	  + r[2*k][2*j-1][2*i] + r[2*k][2*j+1][2*i]
	  + r[2*k-1][2*j][2*i] + r[2*k+1][2*j][2*i];
	f[k][j][i] = (tmp1 + 2*tmp2 + 4*tmp3 + 8*r[2*k][2*j][2*i])/64.;
      }
    }
  }

  //use exactly the same value from fine to coarse grid at boundaries
  for(k=0; k<=gmg[lev+1].nz; k++){
    for(j=0; j<=gmg[lev+1].ny; j++){
      f[k][j][0] = r[k][j][0];
      f[k][j][gmg[lev+1].nx] = r[k][j][gmg[lev].nx];
    }
  }
  for(k=0; k<=gmg[lev+1].nz; k++){
    for(i=0; i<=gmg[lev+1].nx; i++){
      f[k][0][i] = r[k][0][i];
      f[k][gmg[lev+1].ny][i] = r[k][gmg[lev].ny][i];
    }
  }
  for(j=0; j<=gmg[lev+1].ny; j++){
    for(i=0; i<=gmg[lev+1].nx; i++){
      f[0][j][i] = r[0][j][i];
      f[gmg[lev+1].nz][j][i] = r[gmg[lev].nz][j][i];
    }
  }

}

//multigrid V-cycle
void v_cycle(gmg_t *gmg, int lev)
{
  int i, j, k, n;
  double ***r;

  for(i=0; i<v1; i++) smoothing(gmg, lev);//pre-smoothing of u based on u,f at lev-th level

  if(lev<lmax-1){
    r = alloc3double(gmg[lev].nx+1, gmg[lev].ny+1, gmg[lev].nz+1);//residual vector
    residual(gmg, r, lev);//residual r=f-Au at lev-th lev

    if(lev==0){//compute the norm of the residual vector at the beginning of each iteration
      n = (gmg[lev].nx+1)*(gmg[lev].ny+1)*(gmg[lev].nz+1);
      rnorm = sqrt(inner_product(n, &r[0][0][0], &r[0][0][0]));
    }
    
    restriction(gmg, r, lev);//restrict r at lev-th lev to gmg[lev+1].f 

    n = (gmg[lev+1].nx+1)*(gmg[lev+1].ny+1)*(gmg[lev+1].nz+1);
    memset(&gmg[lev+1].u[0][0][0], 0, n*sizeof(double));
    v_cycle(gmg, lev+1);// another v-cycle at (lev+1)-th level

    prolongation(gmg, r, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
    for(k=1; k<gmg[lev].nz; k++)
      for(j=1; j<gmg[lev].ny; j++)
	for(i=1; i<gmg[lev].nx; i++)
	  gmg[lev].u[k][j][i] += r[k][j][i];//correct u=u+r at interior without boundaries
    free3double(r);
  }
  //if lev==lmax-1, then nx=ny=2, grid size=3*3, only 1 point at the center is unknwn
  //direct solve is equivalent to smoothing at center point, one post-smoothing will do the joib
  
  for(i=0; i<v2; i++) smoothing(gmg, lev);//post-smoothing
}

//multigrid F-cycle
void f_cycle(gmg_t *gmg, int lev)
{
  int i, j, k, n;
  double ***r;

  if(lev==lmax-1){//coarsest grid, direct solve or smoothing
    n = (gmg[lev].nx+1)*(gmg[lev].ny+1)*(gmg[lev].nz+1);	
    memset(&gmg[lev].u[0][0][0], 0, n*sizeof(double));
  }else{
    r = alloc3double(gmg[lev].nx+1, gmg[lev].ny+1, gmg[lev].nz+1);
    residual(gmg, r, lev);//residual r=f-Au at lev-th level
    restriction(gmg, r, lev);//restrict r at lev-th lev to f at (lev+1)-th level

    n = (gmg[lev+1].nx+1)*(gmg[lev+1].ny+1)*(gmg[lev+1].nz+1);
    memset(&gmg[lev+1].u[0][0][0], 0, n*sizeof(double));
    f_cycle(gmg, lev+1);
    prolongation(gmg, r, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
    for(k=1; k<gmg[lev].nz; k++)
      for(j=1; j<gmg[lev].ny; j++)
	for(i=1; i<gmg[lev].nx; i++)
	  gmg[lev].u[k][j][i] += r[k][j][i];//correct u=u+r at interior without boundaries
    free3double(r);
  }
  v_cycle(gmg, lev);// another v-cycle at (lev+1)-th level
}

void gmg_init(int nx, int ny, int nz, double dx, double dy, double dz)
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
  
  gmg = (gmg_t*)malloc(lmax*sizeof(gmg_t));
  for(i=0; i<lmax; i++){
    gmg[i].nx = nx/(1<<i);//nx/2^i
    gmg[i].ny = ny/(1<<i);//ny/2^i
    gmg[i].nz = nz/(1<<i);//nz/2^i
    gmg[i].dx = dx*(1<<i);//dx*2^i
    gmg[i].dy = dy*(1<<i);//dy*2^i    
    gmg[i].dz = dz*(1<<i);//dz*2^i    
    gmg[i].u = (double***)alloc3double(gmg[i].nx+1, gmg[i].ny+1, gmg[i].nz+1);
    gmg[i].f = (double***)alloc3double(gmg[i].nx+1, gmg[i].ny+1, gmg[i].nz+1);
  }
  
}

void gmg_close()
{
  int i;
  
  for(i=0; i<lmax; i++){
    free3double(gmg[i].u);
    free3double(gmg[i].f);
  }
  free(gmg);
}

void gmg_apply(int n, double *b, double *x)
{
  int iter;

  memcpy(&gmg[0].f[0][0][0], b, n*sizeof(double));
  memset(&gmg[0].u[0][0][0], 0, n*sizeof(double));

  for(iter=0; iter<itermax; iter++) f_cycle(gmg, 0);

  memcpy(x, &gmg[0].u[0][0][0], n*sizeof(double));
}

