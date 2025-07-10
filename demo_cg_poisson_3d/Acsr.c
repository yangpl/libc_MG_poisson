#include "cstd.h"
#include "sparse.h"

csr_t *A;

#define id(i,j,k) (i + (nx+1)*(j + (ny+1)*(k)))

//build matrix in coo format for laplace operator A=-Dx^2 - Dy^2 - Dz^2
void Acsr_init(csr_t *A_, int nx, int ny, int nz, double dx, double dy, double dz)
{
  int i, j, k, kk, row_ind;
  double _dx2, _dy2, _dz2;
  
  //count the number of nonzeros in sparse banded matrix A
  kk = 0;
  for(k=0; k<=nz; k++){
    for(j=0; j<=ny; j++){
      for(i=0; i<=nx; i++){
	//the diagonal element A(i,j)
	kk++;
      
	//the off-diagonal element
	if(i-1>=0){//element A(i-1,j,k)
	  kk++;	
	}
	if(i+1<=nx){//element A(i+1,j,k)
	  kk++;
	}
	if(j-1>=0){//element A(i,j-1,k)
	  kk++;	
	}
	if(j+1<=ny){//element A(i,j+1,k)
	  kk++;
	}
	if(k-1>=0){//element A(i,j,k-1)
	  kk++;	
	}
	if(k+1<=nz){//element A(i,j,k+1)
	  kk++;
	}
      
      }
    }
  }
  A = A_;
  A->nnz = kk;//number of non-zeros
  A->nrow = (nx+1)*(ny+1)*(nz+1);
  A->ncol = (nx+1)*(ny+1)*(nz+1);
  
  A->row_ptr = alloc1int(A->nrow+1);
  A->col_ind = alloc1int(A->nnz);
  A->val = alloc1double(A->nnz);

  _dx2 = 1./(dx*dx);
  _dy2 = 1./(dy*dy);
  _dz2 = 1./(dz*dz);
  kk = 0;
  A->row_ptr[0] = kk;
  for(k=0; k<=nz; k++){
    for(j=0; j<=ny; j++){
      for(i=0; i<=nx; i++){
	//the diagonal element A(i,j,k)
	A->val[kk] = 2.*(_dx2 + _dy2 + _dz2);
	A->col_ind[kk] = id(i,j,k);
	kk++;
      
	//the off-diagonal element
	if(i-1>=0){//element A(i-1,j,k)
	  A->val[kk] = -_dx2;
	  A->col_ind[kk] = id(i-1,j,k);
	  kk++;	
	}
	if(i+1<=nx){//element A(i+1,j,k)
	  A->val[kk] = -_dx2;
	  A->col_ind[kk] = id(i+1,j,k);
	  kk++;
	}
	if(j-1>=0){//element A(i,j-1,k)
	  A->val[kk] = -_dy2;
	  A->col_ind[kk] = id(i,j-1,k);
	  kk++;	
	}
	if(j+1<=ny){//element A(i,j+1,k)
	  A->val[kk] = -_dy2;
	  A->col_ind[kk] = id(i,j+1,k);
	  kk++;
	}
	if(k-1>=0){//element A(i,j,k-1)
	  A->val[kk] = -_dz2;
	  A->col_ind[kk] = id(i,j,k-1);
	  kk++;	
	}
	if(k+1<=nz){//element A(i,j,k+1)
	  A->val[kk] = -_dz2;
	  A->col_ind[kk] = id(i,j,k+1);
	  kk++;
	}
      
	row_ind = id(i,j,k);
	A->row_ptr[row_ind+1] = kk;
      }
    }
  }
  
}


void Acsr_close()
{
  free1int(A->row_ptr);
  free1int(A->col_ind);
  free1double(A->val);
}

//compute y=Ax by applying operator A in coo format
void Acsr_apply(int n, double *x, double *y)
{
  int i, j, k;

  for(i=0; i<A->nrow; i++){
    y[i] = 0;
    for(k=A->row_ptr[i]; k<A->row_ptr[i+1]; k++){
      j = A->col_ind[k];
      y[i] += A->val[k]*x[j];
    }
  }

}
