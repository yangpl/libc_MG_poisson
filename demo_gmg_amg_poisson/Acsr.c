#include "cstd.h"
#include "sparse.h"

csr_t *A;

//================================================================
//build matrix in coo format for laplace operator A=-Dx^2 - Dy^2
void Acsr_init(csr_t *A_, int nx, int ny, double dx, double dy)
{
  int i, j, k, row_ind;
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
  A = A_;
  A->nnz = k;//number of non-zeros
  A->nrow = (nx+1)*(ny+1);
  A->ncol = (nx+1)*(ny+1);
  
  A->row_ptr = alloc1int(A->nrow+1);
  A->col_ind = alloc1int(A->nnz);
  A->val = alloc1double(A->nnz);

  _dx2 = 1./(dx*dx);
  _dy2 = 1./(dy*dy);
  k = 0;
  for(j=0; j<=ny; j++){
    for(i=0; i<=nx; i++){
      //the diagonal element A(i,j)
      A->val[k] = -2.*(_dx2 + _dy2);
      A->col_ind[k] = i + (nx+1)*j;
      k++;
      
      //the off-diagonal element
      if(i-1>=0){//element A(i-1,j)
	A->val[k] = _dx2;
	A->col_ind[k] = (i-1) + (nx+1)*j;
	k++;	
      }
      if(i+1<=nx){//element A(i+1,j)
	A->val[k] = _dx2;
	A->col_ind[k] = (i+1) + (nx+1)*j;
	k++;
      }
      if(j-1>=0){//element A(i,j-1)
	A->val[k] = _dy2;
	A->col_ind[k] = i + (nx+1)*(j-1);
	k++;	
      }
      if(j+1<=ny){//element A(i,j+1)
	A->val[k] = _dy2;
	A->col_ind[k] = i + (nx+1)*(j+1);
	k++;
      }
      
      row_ind = i + (nx+1)*j;
      A->row_ptr[row_ind+1] = k;
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
