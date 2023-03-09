#include "cstd.h"
#include "sparse.h"

csr_t Alu;
int *kii;

void ilu0_init(csr_t A_)
{
  int i, j, k, k1, k2, k3;
  double tmp;

  Alu.nrow = A_.nrow;
  Alu.ncol = A_.ncol;
  Alu.nnz = A_.nnz;
  Alu.row_ptr = alloc1int(Alu.nrow+1);
  Alu.col_ind = alloc1int(Alu.nnz);
  Alu.val = alloc1double(Alu.nnz);
  memcpy(Alu.row_ptr, A_.row_ptr, (Alu.nrow+1)*sizeof(int));
  memcpy(Alu.col_ind, A_.col_ind, Alu.nnz*sizeof(int));
  memcpy(Alu.val, A_.val, Alu.nnz*sizeof(double));
  
  kii = alloc1int(Alu.nrow);//a working array for bookkeeping nnz index of aii
  for(i=0; i<Alu.nrow; i++){
    for(k1=Alu.row_ptr[i]; k1<Alu.row_ptr[i+1]; k1++){
      k = Alu.col_ind[k1];
      if(k==i) {
	kii[i] = k1; //bookkeep nnz index of diagonal element
	break;
      }
    }
  }

  for(i=1; i<Alu.nrow; i++){
    for(k1=Alu.row_ptr[i]; k1<kii[i]; k1++){
      k = Alu.col_ind[k1];//k=1,...i-1
      Alu.val[k1] /= Alu.val[kii[k]];//a_ik /= a_kk, a_kk=Alu.val[kii[k]]
      for(k2=kii[k]+1; k2<Alu.row_ptr[k+1]; k2++){
	j = Alu.col_ind[k2];//a_kj, j=k+1,...n
	tmp = Alu.val[k1]*Alu.val[k2];//a_ik*a_kj
	for(k3=Alu.row_ptr[i]; k3<Alu.row_ptr[i+1]; k3++){
	  if(j==Alu.col_ind[k3]) Alu.val[k3] -= tmp;//a_ij-=a_ik*a_kj
	}//end for k3
      }//end for k2
    }//end for k1
  }//end for i

}

void ilu0_close()
{
  free1int(kii);
  free1int(Alu.row_ptr);
  free1int(Alu.col_ind);
  free1double(Alu.val);
}


//incomplete LU(0): Ax=b, A=LU
void ilu0_apply(int n, double *b, double *x)
{
  int i, j, k1;
  double s, *y;

  //A=LU, now we solve: Ly=b and Ux = y
  y = x;
  for(i=0; i<Alu.nrow; i++){
    s = 0;
    for(k1=Alu.row_ptr[i]; k1<kii[i]; k1++){
      j = Alu.col_ind[k1];
      s += Alu.val[k1]*y[j];//s=\sum_{j<i} a_ij*y_j
    }
    y[i] = b[i] - s;//y_i = b_i-\sum_{j<i} a_ij*y_j
  }

  for(i=Alu.nrow-1; i>=0; i--){
    s = 0.;
    for(k1=kii[i]+1; k1<Alu.row_ptr[i+1]; k1++){
      j = Alu.col_ind[k1];
      s += Alu.val[k1]*x[j];
    }
    x[i] = (y[i]-s)/Alu.val[kii[i]];//a_ii
  }

}

