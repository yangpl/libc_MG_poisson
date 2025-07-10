#include "cstd.h"
#include "sparse.h"

double dx, dy, dz;
int nx, ny, nz;

//build matrix in coo format for laplace operator A=-Dx^2 - Dy^2 - Dz^2
void matrix_init(int nx_, int ny_, int nz_, double dx_, double dy_, double dz_)
{
  nx = nx_;
  ny = ny_;
  nz = nz_;
  dx = dx_;
  dy = dy_;
  dz = dz_;
}

//compute y=Ax by applying operator A in coo format
void matrix_apply(int n, double *x, double *y)
{
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  double s1, s2, s3;

#define id(i,j,k) (i+(nx+1)*(j + (ny+1)*(k)))
  
  memset(y, 0, n*sizeof(double));
  for(k=0; k<=nz; ++k){
    kp1 = k+1;
    km1 = k-1;
    for(j=0; j<=ny; ++j){
      jp1 = j+1;
      jm1 = j-1;
      for(i=0; i<=nx; ++i){
	ip1 = i+1;
	im1 = i-1;
	
	/* y[id(i,j,k)] += (x[id(ip1,j,k)]-2*x[id(i,j,k)]+x[id(im1,j,k)])/(dx*dx); */
	/* y[id(i,j,k)] += (x[id(i,jp1,k)]-2*x[id(i,j,k)]+x[id(i,jm1,k)])/(dy*dy); */
	/* y[id(i,j,k)] += (x[id(i,j,kp1)]-2*x[id(i,j,k)]+x[id(i,j,km1)])/(dz*dz); */

	s1 = -2*x[id(i,j,k)];
	if(ip1<=nx) s1 += x[id(ip1,j,k)];
	if(im1>=0) s1 += x[id(im1,j,k)];
	s1 /= (dx*dx);

	s2 = -2*x[id(i,j,k)];
	if(jp1<=ny) s2 += x[id(i,jp1,k)];
	if(jm1>=0) s2 += x[id(i,jm1,k)];
	s2 /= (dy*dy);

	s3 = -2*x[id(i,j,k)];
	if(kp1<=nz) s3 += x[id(i,j,kp1)];
	if(km1>=0) s3 += x[id(i,j,km1)];
	s3 /= (dz*dz);

	y[id(i,j,k)] += s1 + s2 + s3;
      }
    }
  }
#undef id
  
}
