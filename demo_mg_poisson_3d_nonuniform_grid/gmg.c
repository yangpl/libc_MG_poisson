/* Demo for solving 3D possion equation using multigrid method
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"

int cycleopt;
int icycle, ncycle;
int v1;//number of pre-smoothing
int v2;//number of post-smoothing
int lmax;//maximum number of levels for switching between fine and coarse grid
double tol;//tolerance for convergence
double rnorm;//norm of residual vector r

typedef struct{
  int n1;
  int n2;
  int n3;
  double *x1, *x2, *x3;//node coordinates along x, y and z
  double *x1s, *x2s, *x3s;//staggered node coordinates along x, y and z
  double *d1, *d2, *d3;//cell sizes centered at integer nodes
  double *d1s, *d2s, *d3s;//cell sizes centered at staggered nodes
  double ***u;
  double ***f;
  double ***r;
  double ***sigma;
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
void gauss_seidel(gmg_t *gmg, int lev, int iter)
{
  int n1, n2, n3;
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  double ***u, ***f;
  double dudx_ip, dudx_im, dudy_jp, dudy_jm, dudz_kp, dudz_km;
  double sigma_ip, sigma_im, sigma_jp, sigma_jm, sigma_kp, sigma_km;
  double Ax, coeff;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  double ***sigma;
  
  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1 = gmg[lev].d1;
  d2 = gmg[lev].d2;
  d3 = gmg[lev].d3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  u = gmg[lev].u;
  f = gmg[lev].f;
  sigma = gmg[lev].sigma;

  for(k=1; k<n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=1; j<n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      for(i=1; i<n1; i++){
	ip1 = i+1;
	im1 = i-1;

	u[k][j][i] = 0;//iniitialize unknowns to be 0, rhs=b-Ax
	
	sigma_ip = 0.25*(sigma[k][j][i] + sigma[k][jm1][i] + sigma[km1][j][i] + sigma[km1][jm1][i]);//sigma_{i+0.5}
	sigma_im = 0.25*(sigma[k][j][im1] + sigma[k][jm1][im1] + sigma[km1][j][im1] + sigma[km1][jm1][im1]);//sigma_{i-0.5}
	dudx_ip = (u[k][j][ip1]-u[k][j][i])/d1s[i];//dudx_{i+0.5}
	dudx_im = (u[k][j][i]-u[k][j][im1])/d1s[im1];//dudx_{i-0.5}
	
	sigma_jp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[km1][j][i] + sigma[km1][j][im1]);//sigma_{j+0.5}
	sigma_jm = 0.25*(sigma[k][jm1][i] + sigma[k][jm1][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//sigma_{j-0.5}
	dudy_jp = (u[k][jp1][i]-u[k][j][i])/d2s[j];//dudy_{j+0.5}
	dudy_jm = (u[k][j][i]-u[k][jm1][i])/d2s[jm1];//dudy_{j-0.5}
	
	sigma_kp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[k][jm1][i] + sigma[k][jm1][im1]);//sigma_{k+0.5}
	sigma_km = 0.25*(sigma[km1][j][i] + sigma[km1][j][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//sigma_{k-0.5}
	dudz_kp = (u[kp1][j][i]-u[k][j][i])/d3s[j];//dudz_{k+0.5}
	dudz_km = (u[k][j][i]-u[km1][j][i])/d3s[jm1];//dudz_{k-0.5}

	Ax = 0;
	Ax += -(sigma_ip*dudx_ip - sigma_im*dudx_im)/d1[i];
	Ax += -(sigma_jp*dudy_jp - sigma_jm*dudy_jm)/d2[j];
	Ax += -(sigma_kp*dudz_kp - sigma_km*dudz_km)/d3[k];

	coeff = 0.;
	coeff += (sigma_ip/d1s[i] + sigma_im/d1s[im1])/d1[i];
	coeff += (sigma_jp/d2s[j] + sigma_jm/d2s[jm1])/d2[j];
	coeff += (sigma_kp/d3s[k] + sigma_km/d3s[km1])/d3[k];

	u[k][j][i] = (f[k][j][i] - Ax)/coeff;
      }
    }
  }

}

/*
  solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
  x - initially contains the input vector v, and returns the solution x. indexed from 0 to n - 1 inclusive
  n - number of equations (length of vector x)
  a - subdiagonal (means it is the diagonal below the main diagonal), indexed from 1 to n - 1 inclusive
  b - the main diagonal, indexed from 0 to n - 1 inclusive
  c - superdiagonal (means it is the diagonal above the main diagonal), indexed from 0 to n - 2 inclusive
     
  Note 1: contents of input vector c will be modified, making this a one-time-use function 
  (scratch space can be allocated instead for this purpose to make it reusable)
  Note 2: We don't check for diagonal dominance, etc.; this is not guaranteed stable
*/
void tridiagonal(int n, double *a, double *b, double *c, double *x)
{
  c[0] = c[0] / b[0];
  x[0] = x[0] / b[0];

  /* loop from 1 to n - 1 inclusive, performing the forward sweep */
  for (int ix = 1; ix < n; ix++) {
    const double m = 1.0f / (b[ix] - a[ix] * c[ix - 1]);
    c[ix] = c[ix] * m;
    x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
  }
    
  /* loop from n - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
  for (int ix = n - 2; ix > 0; ix--)
    x[ix] -= c[ix] * x[ix + 1];
  x[0] -= c[0] * x[1];
}

void gauss_seidel_x(gmg_t *gmg, int lev, int iter)
{
  int n1, n2, n3;
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  double ***u, ***f;
  double dudx_ip, dudx_im, dudy_jp, dudy_jm, dudz_kp, dudz_km;
  double sigma_ip, sigma_im, sigma_jp, sigma_jm, sigma_kp, sigma_km;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  double ***sigma;
  double *a, *b, *c, *r;
  
  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1 = gmg[lev].d1;
  d2 = gmg[lev].d2;
  d3 = gmg[lev].d3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  u = gmg[lev].u;
  f = gmg[lev].f;
  sigma = gmg[lev].sigma;

  a = alloc1double(n1-1);
  b = alloc1double(n1-1);
  c = alloc1double(n1-1);
  r = alloc1double(n1-1);
  for(k=1; k<n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=1; j<n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      //iniitialize unknowns to be 0, rhs=b-Ax
      for(i=1; i<n1; i++) u[k][j][i] = 0;

      for(i=1; i<n1; i++){
	ip1 = i+1;
	im1 = i-1;
	
	sigma_ip = 0.25*(sigma[k][j][i] + sigma[k][jm1][i] + sigma[km1][j][i] + sigma[km1][jm1][i]);//sigma_{i+0.5}
	sigma_im = 0.25*(sigma[k][j][im1] + sigma[k][jm1][im1] + sigma[km1][j][im1] + sigma[km1][jm1][im1]);//sigma_{i-0.5}
	dudx_ip = (u[k][j][ip1]-u[k][j][i])/d1s[i];//dudx_{i+0.5}
	dudx_im = (u[k][j][i]-u[k][j][im1])/d1s[im1];//dudx_{i-0.5}
	
	sigma_jp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[km1][j][i] + sigma[km1][j][im1]);//sigma_{j+0.5}
	sigma_jm = 0.25*(sigma[k][jm1][i] + sigma[k][jm1][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//sigma_{j-0.5}
	dudy_jp = (u[k][jp1][i]-u[k][j][i])/d2s[j];//dudy_{j+0.5}
	dudy_jm = (u[k][j][i]-u[k][jm1][i])/d2s[jm1];//dudy_{j-0.5}
	
	sigma_kp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[k][jm1][i] + sigma[k][jm1][im1]);//sigma_{k+0.5}
	sigma_km = 0.25*(sigma[km1][j][i] + sigma[km1][j][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//sigma_{k-0.5}
	dudz_kp = (u[kp1][j][i]-u[k][j][i])/d3s[j];//dudz_{k+0.5}
	dudz_km = (u[k][j][i]-u[km1][j][i])/d3s[jm1];//dudz_{k-0.5}

	r[im1] = f[k][j][i];
	r[im1] -= -(sigma_ip*dudx_ip - sigma_im*dudx_im)/d1[i];
	r[im1] -= -(sigma_jp*dudy_jp - sigma_jm*dudy_jm)/d2[j];
	r[im1] -= -(sigma_kp*dudz_kp - sigma_km*dudz_km)/d3[k];

	b[im1] = 0;
	b[im1] += (sigma_ip/d1s[i] + sigma_im/d1s[im1])/d1[i];
	b[im1] += (sigma_jp/d2s[j] + sigma_jm/d2s[jm1])/d2[j];
	b[im1] += (sigma_kp/d3s[k] + sigma_km/d3s[km1])/d3[k];

	a[im1] = -sigma_ip/(d1s[i]*d1[i]);
	c[im1] = -sigma_im/(d1s[im1]*d1[i]);
      }
      tridiagonal(n1-1, a, b, c, r);

      for(i=1; i<n1; i++) u[k][j][i] = r[i-1];
    }
  }

  free1double(a);
  free1double(b);
  free1double(c);
  free1double(r);
}

void gauss_seidel_y(gmg_t *gmg, int lev, int iter)
{
  int n1, n2, n3;
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  double ***u, ***f;
  double dudx_ip, dudx_im, dudy_jp, dudy_jm, dudz_kp, dudz_km;
  double sigma_ip, sigma_im, sigma_jp, sigma_jm, sigma_kp, sigma_km;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  double ***sigma;
  double *a, *b, *c, *r;
  
  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1 = gmg[lev].d1;
  d2 = gmg[lev].d2;
  d3 = gmg[lev].d3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  u = gmg[lev].u;
  f = gmg[lev].f;
  sigma = gmg[lev].sigma;

  a = alloc1double(n2-1);
  b = alloc1double(n2-1);
  c = alloc1double(n2-1);
  r = alloc1double(n2-1);
  for(k=1; k<n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(i=1; i<n1; i++){
      ip1 = i+1;
      im1 = i-1;
      //iniitialize unknowns to be 0, rhs=b-Ax
      for(j=1; j<n2; j++) u[k][j][i] = 0;
      
      for(j=1; j<n2; j++){
	jp1 = j+1;
	jm1 = j-1;
	
	sigma_ip = 0.25*(sigma[k][j][i] + sigma[k][jm1][i] + sigma[km1][j][i] + sigma[km1][jm1][i]);//sigma_{i+0.5}
	sigma_im = 0.25*(sigma[k][j][im1] + sigma[k][jm1][im1] + sigma[km1][j][im1] + sigma[km1][jm1][im1]);//sigma_{i-0.5}
	dudx_ip = (u[k][j][ip1]-u[k][j][i])/d1s[i];//dudx_{i+0.5}
	dudx_im = (u[k][j][i]-u[k][j][im1])/d1s[im1];//dudx_{i-0.5}
	
	sigma_jp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[km1][j][i] + sigma[km1][j][im1]);//sigma_{j+0.5}
	sigma_jm = 0.25*(sigma[k][jm1][i] + sigma[k][jm1][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//sigma_{j-0.5}
	dudy_jp = (u[k][jp1][i]-u[k][j][i])/d2s[j];//dudy_{j+0.5}
	dudy_jm = (u[k][j][i]-u[k][jm1][i])/d2s[jm1];//dudy_{j-0.5}
	
	sigma_kp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[k][jm1][i] + sigma[k][jm1][im1]);//sigma_{k+0.5}
	sigma_km = 0.25*(sigma[km1][j][i] + sigma[km1][j][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//sigma_{k-0.5}
	dudz_kp = (u[kp1][j][i]-u[k][j][i])/d3s[j];//dudz_{k+0.5}
	dudz_km = (u[k][j][i]-u[km1][j][i])/d3s[jm1];//dudz_{k-0.5}

	r[jm1] = f[k][j][i];
	r[jm1] -= -(sigma_ip*dudx_ip - sigma_im*dudx_im)/d1[i];
	r[jm1] -= -(sigma_jp*dudy_jp - sigma_jm*dudy_jm)/d2[j];
	r[jm1] -= -(sigma_kp*dudz_kp - sigma_km*dudz_km)/d3[k];

	b[jm1] = 0;
	b[jm1] += (sigma_ip/d1s[i] + sigma_im/d1s[im1])/d1[i];
	b[jm1] += (sigma_jp/d2s[j] + sigma_jm/d2s[jm1])/d2[j];
	b[jm1] += (sigma_kp/d3s[k] + sigma_km/d3s[km1])/d3[k];

	a[jm1] = -sigma_jp/(d2s[j]*d2[j]);
	c[jm1] = -sigma_jm/(d2s[jm1]*d2[j]);
      }
      tridiagonal(n2-1, a, b, c, r);

      for(j=1; j<n2; j++) u[k][j][i] = r[j-1];
    }
  }

  free1double(a);
  free1double(b);
  free1double(c);
  free1double(r);
}


void gauss_seidel_z(gmg_t *gmg, int lev, int iter)
{
  int n1, n2, n3;
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  double ***u, ***f;
  double dudx_ip, dudx_im, dudy_jp, dudy_jm, dudz_kp, dudz_km;
  double sigma_ip, sigma_im, sigma_jp, sigma_jm, sigma_kp, sigma_km;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  double ***sigma;
  double *a, *b, *c, *r;
  
  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1 = gmg[lev].d1;
  d2 = gmg[lev].d2;
  d3 = gmg[lev].d3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  u = gmg[lev].u;
  f = gmg[lev].f;
  sigma = gmg[lev].sigma;
 
  a = alloc1double(n3-1);
  b = alloc1double(n3-1);
  c = alloc1double(n3-1);
  r = alloc1double(n3-1);
  for(j=1; j<n2; j++){
    jp1 = j+1;
    jm1 = j-1;
    for(i=1; i<n1; i++){
      ip1 = i+1;
      im1 = i-1;
      //iniitialize unknowns to be 0, rhs=b-Ax
      for(k=1; k<n3; k++) u[k][j][i] = 0;
      
      for(k=1; k<n3; k++){
	kp1 = k+1;
	km1 = k-1;
	
	sigma_ip = 0.25*(sigma[k][j][i] + sigma[k][jm1][i] + sigma[km1][j][i] + sigma[km1][jm1][i]);//sigma_{i+0.5}
	sigma_im = 0.25*(sigma[k][j][im1] + sigma[k][jm1][im1] + sigma[km1][j][im1] + sigma[km1][jm1][im1]);//sigma_{i-0.5}
	dudx_ip = (u[k][j][ip1]-u[k][j][i])/d1s[i];//dudx_{i+0.5}
	dudx_im = (u[k][j][i]-u[k][j][im1])/d1s[im1];//dudx_{i-0.5}
	
	sigma_jp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[km1][j][i] + sigma[km1][j][im1]);//sigma_{j+0.5}
	sigma_jm = 0.25*(sigma[k][jm1][i] + sigma[k][jm1][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//sigma_{j-0.5}
	dudy_jp = (u[k][jp1][i]-u[k][j][i])/d2s[j];//dudy_{j+0.5}
	dudy_jm = (u[k][j][i]-u[k][jm1][i])/d2s[jm1];//dudy_{j-0.5}
	
	sigma_kp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[k][jm1][i] + sigma[k][jm1][im1]);//sigma_{k+0.5}
	sigma_km = 0.25*(sigma[km1][j][i] + sigma[km1][j][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//sigma_{k-0.5}
	dudz_kp = (u[kp1][j][i]-u[k][j][i])/d3s[j];//dudz_{k+0.5}
	dudz_km = (u[k][j][i]-u[km1][j][i])/d3s[jm1];//dudz_{k-0.5}

	r[km1] = f[k][j][i];
	r[km1] -= -(sigma_ip*dudx_ip - sigma_im*dudx_im)/d1[i];
	r[km1] -= -(sigma_jp*dudy_jp - sigma_jm*dudy_jm)/d2[j];
	r[km1] -= -(sigma_kp*dudz_kp - sigma_km*dudz_km)/d3[k];

	b[km1] = 0;
	b[km1] += (sigma_ip/d1s[i] + sigma_im/d1s[im1])/d1[i];
	b[km1] += (sigma_jp/d2s[j] + sigma_jm/d2s[jm1])/d2[j];
	b[km1] += (sigma_kp/d3s[k] + sigma_km/d3s[km1])/d3[k];

	a[km1] = -sigma_kp/(d3s[k]*d3[k]);
	c[km1] = -sigma_km/(d3s[km1]*d3[k]);
      }
      tridiagonal(n3-1, a, b, c, r);

      for(k=1; k<n3; k++) u[k][j][i] = r[k-1];
    }
  }

  free1double(a);
  free1double(b);
  free1double(c);
  free1double(r);
}

void smoothing(gmg_t *gmg, int lev, int iter)
{
  /*
  gauss_seidel_x(gmg, lev, iter);
  gauss_seidel_y(gmg, lev, iter);
  gauss_seidel_z(gmg, lev, iter);
  */
  gauss_seidel(gmg, lev, iter);
}

//compute residual r=f-Au
void residual(gmg_t *gmg, int lev)
{
  int n1, n2, n3;
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  double ***u, ***f, ***r;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  double ***sigma;
  double sigma_ip, sigma_im, sigma_jp, sigma_jm, sigma_kp, sigma_km;
  double dudx_ip, dudx_im, dudy_jp, dudy_jm, dudz_kp, dudz_km;
  
  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;  
  d1 = gmg[lev].d1;
  d2 = gmg[lev].d2;
  d3 = gmg[lev].d3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  u = gmg[lev].u;
  f = gmg[lev].f;
  r = gmg[lev].r;
  sigma = gmg[lev].sigma;

  for(k=1; k<n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=1; j<n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      for(i=1; i<n1; i++){
	ip1 = i+1;
	im1 = i-1;
	r[k][j][i] = f[k][j][i];
	
	sigma_ip = 0.25*(sigma[k][j][i] + sigma[k][jm1][i] + sigma[km1][j][i] + sigma[km1][jm1][i]);//epsilon_{i+0.5}
	sigma_im = 0.25*(sigma[k][j][im1] + sigma[k][jm1][im1] + sigma[km1][j][im1] + sigma[km1][jm1][im1]);//epsilon_{i-0.5}
	dudx_ip = (u[k][j][ip1]-u[k][j][i])/d1s[i];//dudx_{i+0.5}
	dudx_im = (u[k][j][i]-u[k][j][im1])/d1s[im1];//dudx_{i-0.5}
	r[k][j][i] -= -(sigma_ip*dudx_ip - sigma_im*dudx_im)/d1[i];
		
	sigma_jp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[km1][j][i] + sigma[km1][j][im1]);//epsilon_{j+0.5}
	sigma_jm = 0.25*(sigma[k][jm1][i] + sigma[k][jm1][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//epsilon_{j-0.5}
	dudy_jp = (u[k][jp1][i]-u[k][j][i])/d2s[j];//dudy_{j+0.5}
	dudy_jm = (u[k][j][i]-u[k][jm1][i])/d2s[jm1];//dudy_{j-0.5}
	r[k][j][i] -= -(sigma_jp*dudy_jp - sigma_jm*dudy_jm)/d2[j];
	
	sigma_kp = 0.25*(sigma[k][j][i] + sigma[k][j][im1] + sigma[k][jm1][i] + sigma[k][jm1][im1]);//epsilon_{k+0.5}
	sigma_km = 0.25*(sigma[km1][j][i] + sigma[km1][j][im1] + sigma[km1][jm1][i] + sigma[km1][jm1][im1]);//epsilon_{k-0.5}
	dudz_kp = (u[kp1][j][i]-u[k][j][i])/d3s[j];//dudz_{k+0.5}
	dudz_km = (u[k][j][i]-u[km1][j][i])/d3s[jm1];//dudz_{k-0.5}
	r[k][j][i] -= -(sigma_kp*dudz_kp - sigma_km*dudz_km)/d3[k];
      }
    }
  }

  for(k=0; k<=gmg[lev].n3; k++){
    for(j=0; j<=gmg[lev].n2; j++) {
      r[k][j][0] = 0.;
      r[k][j][gmg[lev].n1] = 0.;
    }
  }

  for(k=0; k<=gmg[lev].n3; k++){
    for(i=0; i<=gmg[lev].n1; i++){
      r[k][0][i] = 0.;
      r[k][gmg[lev].n2][i] = 0.;
    }
  }

  for(j=0; j<=gmg[lev].n2; j++){
    for(i=0; i<=gmg[lev].n1; i++){
      r[0][j][i] = 0.;
      r[gmg[lev].n3][j][i] = 0.;
    }
  }

}

//interpolate u from (lev+1) to lev-th level
void prolongation(gmg_t *gmg, int lev)
{
  double ***e, ***u;
  int i, j, k;

  e = gmg[lev+1].u;
  u = gmg[lev].u;
  for(k=0; k<gmg[lev+1].n3; k++){
    for(j=0; j<gmg[lev+1].n2; j++){
      for(i=0; i<gmg[lev+1].n1; i++){
	u[2*k][2*j][2*i] += e[k][j][i];
	u[2*k][2*j][2*i+1] += 0.5*(e[k][j][i] + e[k][j][i+1]);
	u[2*k][2*j+1][2*i] += 0.5*(e[k][j][i] + e[k][j+1][i]);
	//u[2*k][2*j+1][2*i+1] += 0.5*(0.5*(e[k][j][i] + e[k][j][i+1]) + 0.5*(e[k][j+1][i] + e[k][j+1][i+1]));
	u[2*k][2*j+1][2*i+1] += 0.25*(e[k][j][i] + e[k][j][i+1] + e[k][j+1][i] + e[k][j+1][i+1]);
	u[2*k+1][2*j][2*i] += 0.5*(e[k][j][i] + e[k+1][j][i]);
	//u[2*k+1][2*j][2*i+1] += 0.5*(0.5*(e[k][j][i] + e[k][j][i+1]) + 0.5*(e[k+1][j][i] + e[k+1][j][i+1]));
	u[2*k+1][2*j][2*i+1] += 0.25*(e[k][j][i] + e[k][j][i+1] + e[k+1][j][i] + e[k+1][j][i+1]);
	//u[2*k+1][2*j+1][2*i] += 0.5*(0.5*(e[k][j][i] + e[k][j+1][i]) + 0.5*(e[k+1][j][i] + e[k+1][j+1][i]));
	u[2*k+1][2*j+1][2*i] += 0.25*(e[k][j][i] + e[k][j+1][i] + e[k+1][j][i] + e[k+1][j+1][i]);
	u[2*k+1][2*j+1][2*i+1] += 0.125*(  e[k][j][i] + e[k][j][i+1] + e[k][j+1][i] + e[k][j+1][i+1]
					+ e[k+1][j][i] + e[k+1][j][i+1] + e[k+1][j+1][i] + e[k+1][j+1][i+1]);
      }
    }
  }
  
}

//restrict r from lev to (lev+1)-th lev: fine to coarse grid
void restriction(gmg_t *gmg, int lev)
{
  int i, j, k;
  double tmp1, tmp2, tmp3;
  double ***f, ***r;

  //full weighting operator for restriction
  f = gmg[lev+1].f;
  r = gmg[lev].r;
  for(k=1; k<gmg[lev+1].n3; k++){
    for(j=1; j<gmg[lev+1].n2; j++){
      for(i=1; i<gmg[lev+1].n1; i++){
	tmp1 = r[2*k-1][2*j-1][2*i-1] + r[2*k-1][2*j-1][2*i+1] 
	  + r[2*k-1][2*j+1][2*i-1] + r[2*k-1][2*j+1][2*i+1]
	  + r[2*k+1][2*j-1][2*i-1] + r[2*k+1][2*j-1][2*i+1] 
	  + r[2*k+1][2*j+1][2*i-1] + r[2*k+1][2*j+1][2*i+1];
	tmp2 = r[2*k-1][2*j][2*i-1] + r[2*k-1][2*j][2*i+1]
	  + r[2*k-1][2*j-1][2*i] + r[2*k-1][2*j+1][2*i]
	  + r[2*k+1][2*j][2*i-1] + r[2*k+1][2*j][2*i+1]
	  + r[2*k+1][2*j-1][2*i] + r[2*k+1][2*j+1][2*i]
	  + r[2*k][2*j-1][2*i-1] + r[2*k][2*j+1][2*i-1]
	  + r[2*k][2*j-1][2*i+1] + r[2*k][2*j+1][2*i+1];
	tmp3 = r[2*k][2*j][2*i-1] + r[2*k][2*j][2*i+1]
	  + r[2*k][2*j-1][2*i] + r[2*k][2*j+1][2*i]
	  + r[2*k-1][2*j][2*i] + r[2*k+1][2*j][2*i];
	f[k][j][i] = (tmp1 + 2*tmp2 + 4*tmp3 + 8*r[2*k][2*j][2*i])/64.;
      }
    }
  }

  //use exactly the same value from fine to coarse grid at boundaries
  for(k=0; k<=gmg[lev+1].n3; k++){
    for(j=0; j<=gmg[lev+1].n2; j++){
      f[k][j][0] = r[2*k][j][0];
      f[k][j][gmg[lev+1].n1] = r[2*k][2*j][gmg[lev].n1];
    }
  }
  for(k=0; k<=gmg[lev+1].n3; k++){
    for(i=0; i<=gmg[lev+1].n1; i++){
      f[k][0][i] = r[2*k][0][2*i];
      f[k][gmg[lev+1].n2][i] = r[2*k][gmg[lev].n2][2*i];
    }
  }
  for(j=0; j<=gmg[lev+1].n2; j++){
    for(i=0; i<=gmg[lev+1].n1; i++){
      f[0][j][i] = r[0][2*j][2*i];
      f[gmg[lev+1].n3][j][i] = r[gmg[lev].n3][2*j][2*i];
    }
  }

}


//multigrid V-cycle
void v_cycle(gmg_t *gmg, int lev)
{
  int i, n;

  for(i=0; i<v1; i++) smoothing(gmg, lev, i);//pre-smoothing of u based on u,f at lev-th level
  if(lev<lmax-1){
    residual(gmg, lev);//residual r=f-Au at lev-th lev    
    if(cycleopt==1 && lev==0){//compute the norm of the residual vector at the beginning of each iteration
      n = (gmg[lev].n1+1)*(gmg[lev].n2+1)*(gmg[lev].n3+1);
      rnorm = sqrt(inner_product(n, &gmg[lev].r[0][0][0], &gmg[lev].r[0][0][0]));      
      printf("icycle=%d, residual=%e\n", icycle, rnorm);
    }
    restriction(gmg, lev);//restrict r at lev-th lev to gmg[lev+1].f 

    n = (gmg[lev+1].n1+1)*(gmg[lev+1].n2+1)*(gmg[lev+1].n3+1);
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
    n = (gmg[lev].n1+1)*(gmg[lev].n2+1)*(gmg[lev].n3+1);	
    memset(&gmg[lev].u[0][0][0], 0, n*sizeof(double));
  }else{
    residual(gmg, lev);//residual r=f-Au at lev-th level    
    if(cycleopt==2 && lev==0){//compute the norm of the residual vector at the beginning of each iteration
      n = (gmg[lev].n1+1)*(gmg[lev].n2+1)*(gmg[lev].n3+1);
      rnorm = sqrt(inner_product(n, &gmg[lev].r[0][0][0], &gmg[lev].r[0][0][0]));      
      printf("icycle=%d, residual=%e\n", icycle, rnorm);
    }
    restriction(gmg, lev);//restrict r at lev-th lev to f at (lev+1)-th level

    n = (gmg[lev+1].n1+1)*(gmg[lev+1].n2+1)*(gmg[lev+1].n3+1);
    memset(&gmg[lev+1].u[0][0][0], 0, n*sizeof(double));
    f_cycle(gmg, lev+1);
    prolongation(gmg, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
  }
  v_cycle(gmg, lev);// another v-cycle at (lev+1)-th level
}


//use x1[] to derive x1s[], d1[], d1s[]
void generate_xs_dx(int n1, double *x1, double *x1s, double *d1, double *d1s)
{
  int i;
  
  for(i=0; i<n1; i++) {
    d1s[i] = x1[i+1] - x1[i];
    x1s[i] = 0.5*(x1[i] + x1[i+1]);
  }
  x1s[n1] = x1[n1] + 0.5*d1s[n1-1];
  d1s[n1] = x1s[n1] - x1[n1];

  for(i=1; i<=n1; i++) d1[i] = x1s[i] - x1s[i-1];
  d1[0] = x1s[0] - x1[0];
}


void gmg_init(int n1, int n2, int n3, double d1, double d2, double d3)
{
  int n1_;
  int lev, i, j, k;
  int ii, jj, kk;
  int iip1, jjp1, kkp1;
  double vol, val;
  
  if(!getparint("ncycle", &ncycle)) ncycle = 30;/* maximum number of iterations */  
  if(!getparint("v1", &v1)) v1 = 1;/* number of pre-smoothing */
  if(!getparint("v2", &v2)) v2 = 1;/* number of post-smoothing */
  if(!getparint("cycleopt", &cycleopt)) cycleopt = 1;//1=v cycle; 2=f cycle
  if(!getpardouble("tol", &tol)) tol = 1e-6;/* stopping criteria */
  if(!getparint("lmax", &lmax)) {
    lmax = 0;
    n1_ = n1;
    while(n1_>1) {n1_=n1_>>1; lmax++;}
  }
  
  gmg = malloc(lmax*sizeof(gmg_t));
  for(lev=0; lev<lmax; lev++){
    printf("lev=%d\n", lev);
    
    gmg[lev].n1 = n1/(1<<lev);//n1/2^i
    gmg[lev].n2 = n2/(1<<lev);//n2/2^i
    gmg[lev].n3 = n3/(1<<lev);//n3/2^i
    gmg[lev].x1 = alloc1double(gmg[lev].n1+1);
    gmg[lev].x2 = alloc1double(gmg[lev].n2+1);
    gmg[lev].x3 = alloc1double(gmg[lev].n3+1);
    gmg[lev].x1s = alloc1double(gmg[lev].n1+1);
    gmg[lev].x2s = alloc1double(gmg[lev].n2+1);
    gmg[lev].x3s = alloc1double(gmg[lev].n3+1);
    gmg[lev].d1 = alloc1double(gmg[lev].n1+1);
    gmg[lev].d2 = alloc1double(gmg[lev].n2+1);
    gmg[lev].d3 = alloc1double(gmg[lev].n3+1);
    gmg[lev].d1s = alloc1double(gmg[lev].n1+1);
    gmg[lev].d2s = alloc1double(gmg[lev].n2+1);
    gmg[lev].d3s = alloc1double(gmg[lev].n3+1);
    for(i=0; i<=gmg[lev].n1; i++) gmg[lev].x1[i] = (lev==0)?i*d1:gmg[lev-1].x1[2*i];
    for(j=0; j<=gmg[lev].n2; j++) gmg[lev].x2[j] = (lev==0)?j*d2:gmg[lev-1].x2[2*j];
    for(k=0; k<=gmg[lev].n3; k++) gmg[lev].x3[k] = (lev==0)?k*d3:gmg[lev-1].x3[2*k];
    generate_xs_dx(gmg[lev].n1, gmg[lev].x1, gmg[lev].x1s, gmg[lev].d1, gmg[lev].d1s);
    generate_xs_dx(gmg[lev].n2, gmg[lev].x2, gmg[lev].x2s, gmg[lev].d2, gmg[lev].d2s);
    generate_xs_dx(gmg[lev].n3, gmg[lev].x3, gmg[lev].x3s, gmg[lev].d3, gmg[lev].d3s);

    gmg[lev].sigma = alloc3double(gmg[lev].n1, gmg[lev].n2, gmg[lev].n3);
    if(lev==0){
      for(k=0; k<gmg[lev].n3; k++){
	for(j=0; j<gmg[lev].n2; j++){
	  for(i=0; i<gmg[lev].n1; i++){
	    gmg[lev].sigma[k][j][i] = 1.;
	  }
	}
      }
    }else{//restrict model from fine grid to coarse grid
      float maxval1 = 0;
      float maxval2 = 0;
      float maxval3 = 0.;
      for(k=0; k<gmg[lev].n3; k++){
	kk = 2*k;
	kkp1 = 2*k+1;
	for(j=0; j<gmg[lev].n2; j++){
	  jj = 2*j;
	  jjp1 = 2*j+1;
	  for(i=0; i<gmg[lev].n1; i++){
	    ii = 2*i;
	    iip1 = 2*i+1;

	    maxval1 = MAX(maxval1, MAX(MAX(gmg[lev].d1[i], gmg[lev].d2[j]), gmg[lev].d3[k]));
	    maxval2 = MAX(maxval2, MAX(MAX(gmg[lev].d1s[i], gmg[lev].d2s[j]), gmg[lev].d3s[k]));
	    
	    vol = gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kk]
	      + gmg[lev-1].d1s[iip1]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kk]
	      + gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jjp1]*gmg[lev-1].d3s[kk]
	      + gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kk]
	      + gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kkp1]
	      + gmg[lev-1].d1s[iip1]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kkp1]
	      + gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jjp1]*gmg[lev-1].d3s[kkp1]
	      + gmg[lev-1].d1s[iip1]*gmg[lev-1].d2s[jjp1]*gmg[lev-1].d3s[kkp1];
	    val = gmg[lev-1].sigma[kk][jj][ii]*gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kk]
	      + gmg[lev-1].sigma[kk][jj][iip1]*gmg[lev-1].d1s[iip1]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kk]
	      + gmg[lev-1].sigma[kk][jjp1][ii]*gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jjp1]*gmg[lev-1].d3s[kk]
	      + gmg[lev-1].sigma[kk][jjp1][iip1]*gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kk]
	      + gmg[lev-1].sigma[kkp1][jj][ii]*gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kkp1]
	      + gmg[lev-1].sigma[kkp1][jj][iip1]*gmg[lev-1].d1s[iip1]*gmg[lev-1].d2s[jj]*gmg[lev-1].d3s[kkp1]
	      + gmg[lev-1].sigma[kkp1][jjp1][ii]*gmg[lev-1].d1s[ii]*gmg[lev-1].d2s[jjp1]*gmg[lev-1].d3s[kkp1]
	      + gmg[lev-1].sigma[kkp1][jjp1][iip1]*gmg[lev-1].d1s[iip1]*gmg[lev-1].d2s[jjp1]*gmg[lev-1].d3s[kkp1];
	    gmg[lev].sigma[k][j][i] = val/vol;
	    maxval3 = MAX(maxval3, gmg[lev].sigma[k][j][i]);
	  }
	}
      }
    }
    gmg[lev].u = alloc3double(gmg[lev].n1+1, gmg[lev].n2+1, gmg[lev].n3+1);
    gmg[lev].f = alloc3double(gmg[lev].n1+1, gmg[lev].n2+1, gmg[lev].n3+1);
    gmg[lev].r = alloc3double(gmg[lev].n1+1, gmg[lev].n2+1, gmg[lev].n3+1);
  }
  
}

void gmg_close()
{
  int lev;
  
  for(lev=0; lev<lmax; lev++){
    free1double(gmg[lev].x1);
    free1double(gmg[lev].x2);
    free1double(gmg[lev].x3);
    free1double(gmg[lev].x1s);
    free1double(gmg[lev].x2s);
    free1double(gmg[lev].x3s);
    free1double(gmg[lev].d1);
    free1double(gmg[lev].d2);
    free1double(gmg[lev].d3);
    free1double(gmg[lev].d1s);
    free1double(gmg[lev].d2s);
    free1double(gmg[lev].d3s);
    free3double(gmg[lev].sigma);
    free3double(gmg[lev].u);
    free3double(gmg[lev].f);
    free3double(gmg[lev].r);
  }
  free(gmg);
}

void gmg_apply(int n, double *b, double *x)
{
  memcpy(&gmg[0].f[0][0][0], b, n*sizeof(double));
  memset(&gmg[0].u[0][0][0], 0, n*sizeof(double));
  for(icycle=0; icycle<ncycle; icycle++){
    if(cycleopt==1) v_cycle(gmg, 0);
    if(cycleopt==2) f_cycle(gmg, 0);
  }
  memcpy(x, &gmg[0].u[0][0][0], n*sizeof(double));
}

