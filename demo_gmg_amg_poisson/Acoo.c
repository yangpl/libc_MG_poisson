#include "cstd.h"
#include "sparse.h"

coo_t Acoo;

//build matrix in coo format for laplace operator A=-Dx^2 - Dy^2
void Acoo_init(int nx, int ny, double dx, double dy)
{
  int i, j, k;
  double _dx2, _dy2;
  
  //count the number of nonzeros in sparse banded matrix A
  k = 0;
  for(j=0; j<=ny; j++){
    for(i=0; i<=nx; i++){
      //the diagonal element A(i,j)
      k++;
      
      //the off-diagonal element
      if(j-1>=0){//element A(i,j-1)
	k++;	
      }
      if(i-1>=0){//element A(i-1,j)
	k++;	
      }
      if(i+1<=nx){//element A(i+1,j)
	k++;
      }
      if(j+1<=ny){//element A(i,j+1)
	k++;
      }      
    }
  }
  Acoo.nnz = k;//number of non-zeros
  Acoo.nrow = (nx+1)*(ny+1);
  Acoo.ncol = (nx+1)*(ny+1);
  
  
  Acoo.row = alloc1int(Acoo.nnz);
  Acoo.col = alloc1int(Acoo.nnz);
  Acoo.val = alloc1double(Acoo.nnz);

  _dx2 = 1./(dx*dx);
  _dy2 = 1./(dy*dy);
  k = 0;
  for(j=0; j<=ny; j++){
    for(i=0; i<=nx; i++){
      //the diagonal element A(i,j)
      Acoo.val[k] = -2.*(_dx2 + _dy2);
      Acoo.row[k] = i + (nx+1)*j;
      Acoo.col[k] = i + (nx+1)*j;
      k++;
      
      //the off-diagonal element
      if(i-1>=0){//element A(i-1,j)
	Acoo.val[k] = _dx2;
	Acoo.row[k] = i + (nx+1)*j;
	Acoo.col[k] = (i-1) + (nx+1)*j;
	k++;	
      }
      if(i+1<=nx){//element A(i+1,j)
	Acoo.val[k] = _dx2;
	Acoo.row[k] = i + (nx+1)*j;
	Acoo.col[k] = (i+1) + (nx+1)*j;
	k++;
      }
      if(j-1>=0){//element A(i,j-1)
	Acoo.val[k] = _dy2;
	Acoo.row[k] = i + (nx+1)*j;
	Acoo.col[k] = i + (nx+1)*(j-1);
	k++;	
      }
      if(j+1<=ny){//element A(i,j+1)
	Acoo.val[k] = _dy2;
	Acoo.row[k] = i + (nx+1)*j;
	Acoo.col[k] = i + (nx+1)*(j+1);
	k++;
      }      
    }
  }
}

void Acoo_close()
{
  free1int(Acoo.row);
  free1int(Acoo.col);
  free1double(Acoo.val);
}

//compute y=Ax by applying operator A in coo format
void Acoo_apply(int n, double *x, double *y)
{
  int i, j, k;

  memset(y, 0, n*sizeof(double));  
  for(k=0; k<Acoo.nnz; k++){
    i = Acoo.row[k];
    j = Acoo.col[k];
    y[i] += Acoo.val[k]*x[j];//y_i += a_ij*x_j
  }
}
