/* Demo for solving 2D possion equation using multigrid method
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"

int cycleopt;
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
  double ***r;
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
//Symmetric Gauss-Seidel for swtiching even and odd indices
void smoothing(gmg_t *gmg, int lev, int iter)
{
  int i, j, k;
  double vol, a, tmp1, tmp2, tmp3, _dx2, _dy2, _dz2;
  double ***u, ***f;

  u = gmg[lev].u;
  f = gmg[lev].f;
  _dx2 = 1./(gmg[lev].dx*gmg[lev].dx);
  _dy2 = 1./(gmg[lev].dy*gmg[lev].dy);
  _dz2 = 1./(gmg[lev].dz*gmg[lev].dz);
  a = 2.*_dx2 + 2.*_dy2 + 2.*_dz2;
  vol = gmg[lev].dx*gmg[lev].dy*gmg[lev].dz;
  for(k=1; k<gmg[lev].nz; k++){
    if(iter%2==1) k = gmg[lev].nz-k;
    for(j=1; j<gmg[lev].ny; j++){
      if(iter%2==1) j = gmg[lev].ny-j;
      for(i=1; i<gmg[lev].nx; i++){
	if(iter%2==1) i = gmg[lev].nx-i;
	tmp1 = (u[k][j][i+1] + u[k][j][i-1])*_dx2;
	tmp2 = (u[k][j+1][i] + u[k][j-1][i])*_dy2;
	tmp3 = (u[k+1][j][i] + u[k-1][j][i])*_dz2;
	u[k][j][i] = (f[k][j][i] + (tmp1 + tmp2 + tmp3)*vol)/(a*vol);
      }
    }
  }

}

//compute residual r=f-Au
void residual(gmg_t *gmg, int lev)
{
  int i, j, k;
  double vol, _dx2, _dy2, _dz2, tmp1, tmp2, tmp3;
  double ***u, ***f, ***r;
  
  _dx2 = 1./(gmg[lev].dx*gmg[lev].dx);
  _dy2 = 1./(gmg[lev].dy*gmg[lev].dy);
  _dz2 = 1./(gmg[lev].dz*gmg[lev].dz);  
  u = gmg[lev].u;
  f = gmg[lev].f;
  r = gmg[lev].r;
  vol = gmg[lev].dx*gmg[lev].dy*gmg[lev].dz;
  for(k=1; k<gmg[lev].nz; k++){
    for(j=1; j<gmg[lev].ny; j++){
      for(i=1; i<gmg[lev].nx; i++){
	tmp1 = (u[k][j][i+1] - 2.*u[k][j][i] + u[k][j][i-1])*_dx2;
	tmp2 = (u[k][j+1][i] - 2.*u[k][j][i] + u[k][j-1][i])*_dy2;
	tmp3 = (u[k+1][j][i] - 2.*u[k][j][i] + u[k-1][j][i])*_dz2;
	r[k][j][i] = f[k][j][i] + vol*(tmp1 + tmp2 + tmp3);
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
void prolongation(gmg_t *gmg, int lev)
{
  double ***e, ***u;
  int i, j, k;
  int ip1, jp1, kp1;
  int ii, jj, kk;
  int iip1, jjp1, kkp1;
    

  e = gmg[lev+1].u;
  u = gmg[lev].u;
  for(k=0; k<gmg[lev+1].nz; k++){
    kp1 = k + 1;
    kk = 2*k;
    kkp1 = kk + 1;
    for(j=0; j<gmg[lev+1].ny; j++){
      jp1 = j + 1;
      jj = 2*j;
      jjp1 = jj + 1;
      for(i=0; i<gmg[lev+1].nx; i++){
	ip1 = i + 1;
	ii = 2*i;
	iip1 = ii + 1;
	
	u[kk][jj][ii] += e[k][j][i];
	u[kk][jj][iip1] += 0.5*(e[k][j][i] + e[k][j][ip1]);
	u[kk][jjp1][ii] += 0.5*(e[k][j][i] + e[k][jp1][i]);
	//u[kk][jjp1][iip1] += 0.5*(0.5*(e[k][j][i] + e[k][j][ip1]) + 0.5*(e[k][jp1][i] + e[k][jp1][ip1]));
	u[kk][jjp1][iip1] += 0.25*(e[k][j][i] + e[k][j][ip1] + e[k][jp1][i] + e[k][jp1][ip1]);
	u[kkp1][jj][ii] += 0.5*(e[k][j][i] + e[kp1][j][i]);
	//u[kkp1][jj][iip1] += 0.5*(0.5*(e[k][j][i] + e[k][j][ip1]) + 0.5*(e[kp1][j][i] + e[kp1][j][ip1]));
	u[kkp1][jj][iip1] += 0.25*(e[k][j][i] + e[k][j][ip1] + e[kp1][j][i] + e[kp1][j][ip1]);
	//u[kkp1][jjp1][ii] += 0.5*(0.5*(e[k][j][i] + e[k][jp1][i]) + 0.5*(e[kp1][j][i] + e[kp1][jp1][i]));
	u[kkp1][jjp1][ii] += 0.25*(e[k][j][i] + e[k][jp1][i] + e[kp1][j][i] + e[kp1][jp1][i]);
	u[kkp1][jjp1][iip1] += 0.125*(  e[k][j][i] + e[k][j][ip1] + e[k][jp1][i] + e[k][jp1][ip1]
					+ e[kp1][j][i] + e[kp1][j][ip1] + e[kp1][jp1][i] + e[kp1][jp1][ip1]);
      }
    }
  }
  
}

//restrict r from lev to (lev+1)-th lev: fine to coarse grid
void restriction(gmg_t *gmg, int lev) //exact adjoint of prolongation	
{
  int i, j, k;
  int ip1, jp1, kp1;
  int ii, jj, kk;
  int iip1, jjp1, kkp1;
  double tmp1, tmp2, tmp3;
  double ***f, ***r;

  //full weighting operator for restriction
  f = gmg[lev+1].f;
  r = gmg[lev].r;
  memset(&f[0][0][0], 0, (gmg[lev+1].nx+1)*(gmg[lev+1].ny+1)*(gmg[lev+1].nz+1)*sizeof(double));
  for(k=0; k<gmg[lev+1].nz; k++){
    kp1 = k + 1;
    kk = 2*k;
    kkp1 = kk + 1;
    for(j=0; j<gmg[lev+1].ny; j++){
      jp1 = j + 1;
      jj = 2*j;
      jjp1 = jj + 1;
      for(i=0; i<gmg[lev+1].nx; i++){
	ip1 = i + 1;
	ii = 2*i;
	iip1 = ii + 1;
	
	f[k][j][i] += r[kk][jj][ii];
	
	f[k][j][i] += 0.5*r[kk][jj][iip1];
	f[k][j][ip1] += 0.5*r[kk][jj][iip1];

	f[k][j][i] += 0.5*r[kk][jjp1][ii];
	f[k][jp1][i] += 0.5*r[kk][jjp1][ii];

	f[k][j][i] += 0.25*r[kk][jjp1][iip1];
	f[k][j][ip1] += 0.25*r[kk][jjp1][iip1];
	f[k][jp1][i] += 0.25*r[kk][jjp1][iip1];
	f[k][jp1][ip1] += 0.25*r[kk][jjp1][iip1];

	f[k][j][i] += 0.5*r[kkp1][jj][ii];
	f[kp1][j][i] += 0.5*r[kkp1][jj][ii];

	f[k][j][i] += 0.25*r[kkp1][jj][iip1];
	f[k][j][ip1] += 0.25*r[kkp1][jj][iip1];
	f[kp1][j][i] += 0.25*r[kkp1][jj][iip1];
	f[kp1][j][ip1] += 0.25*r[kkp1][jj][iip1];

	f[k][j][i] += 0.25*r[kkp1][jjp1][ii];
	f[k][jp1][i] += 0.25*r[kkp1][jjp1][ii];
	f[kp1][j][i] += 0.25*r[kkp1][jjp1][ii];
	f[kp1][jp1][i] += 0.25*r[kkp1][jjp1][ii];

	f[k][j][i] += 0.125*r[kkp1][jjp1][iip1];
	f[k][j][ip1] += 0.125*r[kkp1][jjp1][iip1];
	f[k][jp1][i] += 0.125*r[kkp1][jjp1][iip1];
	f[k][jp1][ip1] += 0.125*r[kkp1][jjp1][iip1];
	f[kp1][j][i] += 0.125*r[kkp1][jjp1][iip1];
	f[kp1][j][ip1] += 0.125*r[kkp1][jjp1][iip1];
	f[kp1][jp1][i] += 0.125*r[kkp1][jjp1][iip1];
	f[kp1][jp1][ip1] += 0.125*r[kkp1][jjp1][iip1];		
      }
    }
  }

}


//multigrid V-cycle
void v_cycle(gmg_t *gmg, int lev)
{
  int i, n;

  //if lev==lmax-1, then nx=ny=2, grid size=3*3, only 1 point at the center is unknwn
  //direct solve is equivalent to smoothing at center point, one post-smoothing will do the joib
  for(i=0; i<v1; i++) smoothing(gmg, lev, i);//pre-smoothing of u based on u,f at lev-th level
  if(lev<lmax-1){
    residual(gmg, lev);//residual r=f-Au at lev-th lev    
    if(cycleopt==1 && lev==0){//compute the norm of the residual vector at the beginning of each iteration
      n = (gmg[lev].nx+1)*(gmg[lev].ny+1)*(gmg[lev].nz+1);
      rnorm = sqrt(inner_product(n, &gmg[lev].r[0][0][0], &gmg[lev].r[0][0][0]));      
      printf("residual=%e\n", rnorm);
    }
    restriction(gmg, lev);//restrict r at lev-th lev to gmg[lev+1].f 

    n = (gmg[lev+1].nx+1)*(gmg[lev+1].ny+1)*(gmg[lev+1].nz+1);
    memset(&gmg[lev+1].u[0][0][0], 0, n*sizeof(double));
    v_cycle(gmg, lev+1);// another v-cycle at (lev+1)-th level

    prolongation(gmg, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
  }
  for(i=0; i<v2; i++) smoothing(gmg, lev, i);//post-smoothing
}

//multigrid F-cycle
void f_cycle(gmg_t *gmg, int lev)
{
  int n;
  
  if(lev==lmax-1){//coarsest grid, direct solve or smoothing
    n = (gmg[lev].nx+1)*(gmg[lev].ny+1)*(gmg[lev].nz+1);	
    memset(&gmg[lev].u[0][0][0], 0, n*sizeof(double));
  }else{
    residual(gmg, lev);//residual r=f-Au at lev-th level    
    if(cycleopt==2 && lev==0){//compute the norm of the residual vector at the beginning of each iteration
      n = (gmg[lev].nx+1)*(gmg[lev].ny+1)*(gmg[lev].nz+1);
      rnorm = sqrt(inner_product(n, &gmg[lev].r[0][0][0], &gmg[lev].r[0][0][0]));      
      printf("residual=%e\n", rnorm);
    }
    restriction(gmg, lev);//restrict r at lev-th lev to f at (lev+1)-th level

    n = (gmg[lev+1].nx+1)*(gmg[lev+1].ny+1)*(gmg[lev+1].nz+1);
    memset(&gmg[lev+1].u[0][0][0], 0, n*sizeof(double));
    f_cycle(gmg, lev+1);
    prolongation(gmg, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
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
  if(!getparint("cycleopt", &cycleopt)) cycleopt = 1;//1=v cycle; 2=f cycle
  if(!getpardouble("tol", &tol)) tol = 1e-6;/* stopping criteria */
  if(!getparint("lmax", &lmax)) {
    lmax = 0;
    nx_ = nx;
    while(nx_>1) {nx_=nx_>>1; lmax++;}
  }
  
  gmg = malloc(lmax*sizeof(gmg_t));
  for(i=0; i<lmax; i++){
    gmg[i].nx = nx/(1<<i);//nx/2^i
    gmg[i].ny = ny/(1<<i);//ny/2^i
    gmg[i].nz = nz/(1<<i);//nz/2^i
    gmg[i].dx = dx*(1<<i);//dx*2^i
    gmg[i].dy = dy*(1<<i);//dy*2^i    
    gmg[i].dz = dz*(1<<i);//dz*2^i    
    gmg[i].u = alloc3double(gmg[i].nx+1, gmg[i].ny+1, gmg[i].nz+1);
    gmg[i].f = alloc3double(gmg[i].nx+1, gmg[i].ny+1, gmg[i].nz+1);
    gmg[i].r = alloc3double(gmg[i].nx+1, gmg[i].ny+1, gmg[i].nz+1);
  }
  
}

void gmg_close()
{
  int i;
  
  for(i=0; i<lmax; i++){
    free3double(gmg[i].u);
    free3double(gmg[i].f);
    free3double(gmg[i].r);
  }
  free(gmg);
}

void gmg_apply(int n, double *b, double *x)
{
  int i, j, k, iter;

  memcpy(&gmg[0].f[0][0][0], b, n*sizeof(double));
  memset(&gmg[0].u[0][0][0], 0, n*sizeof(double));
  for(k=0; k<=gmg[0].nz; k++){
    for(j=0; j<=gmg[0].ny; j++){
      for(i=0; i<=gmg[0].nx; i++){
	gmg[0].f[k][j][i] *= gmg[0].dx*gmg[0].dy*gmg[0].dz;
      }
    }
  }
  for(iter=0; iter<itermax; iter++){
    if(cycleopt==1) v_cycle(gmg, 0);
    if(cycleopt==2) f_cycle(gmg, 0);
  }
  memcpy(x, &gmg[0].u[0][0][0], n*sizeof(double));
}

