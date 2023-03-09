#include "cstd.h"
#include "sparse.h"

void gaussian_elimination_ikj(int n, double **a, double *x, double *b)
{
  int i, j, k;
  double s, *y = x;

  //A=LU
  for(i=1; i<n; i++){
    for(k=0; k<i; k++){
      a[i][k] /= a[k][k];
      for(j=k+1; j<n; j++){
	a[i][j] -= a[i][k]*a[k][j];
      }
    }
  }

  //Ly=b
  for(i=0; i<n; i++){
    s = 0;
    for(j=0; j<i; j++){
      s += a[i][j]*y[j];
    }
    y[i] = b[i] - s;
  }

  //Ux=y
  for(i=n-1; i>=0; i--){
    s = 0;
    for(j=i+1; j<n; j++){
      s += a[i][j]*x[j];
    }
    x[i] = (y[i]-s)/a[i][i];
  }
  
}

//incomplete LU(0)
void ilu0(csr_t A_, double *x, double *b)
{
  int i, j, k, k1, k2, k3;
  int *kii;
  double *y, tmp;
  csr_t A;

  A.nrow = A_.nrow;
  A.ncol = A_.ncol;
  A.nnz = A_.nnz;
  A.row_ptr = alloc1int(A.nrow+1);
  A.col_ind = alloc1int(A.nnz);
  A.val = alloc1double(A.nnz);
  memcpy(A.row_ptr, A_.row_ptr, (A.nrow+1)*sizeof(int));
  memcpy(A.col_ind, A_.col_ind, A.nnz*sizeof(int));
  memcpy(A.val, A_.val, A.nnz*sizeof(double));
  
  kii = alloc1int(A.nrow);//a working array for bookkeeping nnz index of aii
  for(i=0; i<A.nrow; i++){
    for(k1=A.row_ptr[i]; k1<A.row_ptr[i+1]; k1++){
      k = A.col_ind[k1];
      if(k==i) {
	kii[i] = k1; //bookkeep nnz index of diagonal element
	break;
      }
    }
  }

  for(i=1; i<A.nrow; i++){
    for(k1=A.row_ptr[i]; k1<kii[i]; k1++){
      k = A.col_ind[k1];//k=1,...i-1
      A.val[k1] /= A.val[kii[k]];//a_ik /= a_kk, a_kk=A.val[kii[k]]
      for(k2=kii[k]+1; k2<A.row_ptr[k+1]; k2++){
	j = A.col_ind[k2];//a_kj, j=k+1,...n
	tmp = A.val[k1]*A.val[k2];//a_ik*a_kj
	for(k3=A.row_ptr[i]; k3<A.row_ptr[i+1]; k3++){
	  if(j==A.col_ind[k3]) A.val[k3] -= tmp;//a_ij-=a_ik*a_kj
	}//end for k3
      }//end for k2
    }//end for k1
  }//end for i

  //A=LU, now we solve: Ly=b and Ux = y
  y = x;
  for(i=0; i<A.nrow; i++){
    tmp = 0;
    for(k1=A.row_ptr[i]; k1<kii[i]; k1++){
      j = A.col_ind[k1];
      tmp += A.val[k1]*y[j];//y_i = b_i-\sum_{j<i} a_ij*y_j
    }
    y[i] = b[i]-tmp;
  }

  for(i=A.nrow-1; i>=0; i--){
    tmp = 0.;
    for(k1=kii[i]+1; k1<A.row_ptr[i+1]; k1++){
      j = A.col_ind[k1];
      tmp += A.val[k1]*x[j];
    }
    x[i] = (y[i]-tmp)/A.val[kii[i]];//a_ii
  }

  free1int(kii);
  free1int(A.row_ptr);
  free1int(A.col_ind);
  free1double(A.val);
  
}


void gauss_seidel(int n, csr_t A, double *x, double *b)
{
  int i, j, k;
  double aii;

  aii = 1;
  for(i=0; i<A.nrow; i++){
    x[i] = b[i];
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      if(j!=i) x[i] -= A.val[k]*x[j];// x[i] -= a_ij*x_j, (j!=i)
      else aii = A.val[k];
    }
    x[i] /= aii;
  }
}

//symmetric successive over relaxation
void ssor(int n, csr_t A, double *x, double *b)
{
  int i, j, k;
  double aii;

  aii = 1;
  for(i=0; i<A.nrow; i++){
    x[i] = b[i];
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      if(j!=i) x[i] -= A.val[k]*x[j];// x[i] -= a_ij*x_j, (j!=i)
      else aii = A.val[k];
    }
    x[i] /= aii;
  }

  for(i=A.nrow-1; i>=0; i--){
    x[i] = b[i];
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      if(j!=i) x[i] -= A.val[k]*x[j];// x[i] -= a_ij*x_j, (j!=i)
      else aii = A.val[k];
    }
    x[i] /= aii;
  }
  
}


void sor(int n, csr_t A, double *x, double *b, double w)
{
  int i, j, k;
  double aii, xi;

  aii = 1;
  for(i=0; i<A.nrow; i++){
    xi = x[i];
    x[i] = b[i];
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      if(j!=i) x[i] -= A.val[k]*x[j];// x[i] -= a_ij*x_j, (j!=i)
      else aii = A.val[k];
    }
    x[i] /= aii;
    x[i] = (1.-w)*xi + w*x[i];
  }
}

int main()
{
  int i, j, k, n;
  csr_t A;
  double **B;

  n = 3;
  B = alloc2double(n, n);
  
  B[0][0] = 4;
  B[0][1] = 2;
  B[0][2] = 3;
  B[1][0] = 3;
  B[1][1] = -5;
  B[1][2] = 2;
  B[2][0] = -2;
  B[2][1] = 3;
  B[2][2] = 8;
  

  A.nrow = n;
  A.ncol = n;
  A.nnz = n*n;
  A.row_ptr = alloc1int(A.nrow+1);
  A.col_ind = alloc1int(A.nnz);
  A.val = alloc1double(A.nnz);

  A.row_ptr[0] = 0;
  k = 0;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      A.col_ind[k] = j;

      /* if(i==j) A.val[k] = 4; */
      /* else if(abs(i-j)==1) A.val[k] = -1; */
      /* else     A.val[k] = 0; */
      A.val[k] = B[i][j];
      k++;
    }
    A.row_ptr[i+1] = k;
  }

  printf("-------A---------\n");
  for(i=0; i<A.nrow; i++){
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      printf("%g \t", A.val[k]);
    }
    printf("\n");
  }

  double *x, *b;
  x = alloc1double(n);
  b = alloc1double(n);

  x[0] = -1;
  x[1] = 3;
  x[2] = 2;
  printf("------b-----\n");  
  for(i=0; i<A.nrow; i++){
    b[i] = 0;
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];
      b[i] += A.val[k]*x[j];
    }
    printf("%g\n", b[i]);	
  }

  for(i=0; i<n; i++) x[i] = 0;
  //for(j=0; j<100; j++)  gauss_seidel(n, A, x, b);
  //for(j=0; j<100; j++)  sor(n, A, x, b, 0.5);
  for(j=0; j<100; j++) ssor(n, A, x, b);

  printf("-----x---\n");
  for(i=0; i<n; i++) printf("%g\n", x[i]);

  for(i=0; i<n; i++) x[i] = 0;
  //gaussian_elimination_ikj(n, B, x, b);  
  ilu0(A, x, b);
  
  printf("-----x---\n");
  for(i=0; i<n; i++) printf("%g\n", x[i]);
  
}
