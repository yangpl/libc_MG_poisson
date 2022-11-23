#include "cstd.h"
#include "sparse.h"

icsr_t icsr_transpose(icsr_t S)
{
  int i, j, k, m;
  int *pos;
  icsr_t St;

  St.nrow = S.ncol;
  St.ncol = S.nrow;
  St.nnz = S.nnz;
  St.row_ptr = alloc1int(St.nrow+1);
  St.col_ind = alloc1int(St.nnz);
  St.val = alloc1int(St.nnz);
  
  memset(St.row_ptr, 0, (St.nrow+1)*sizeof(int));
  for(i=0; i<S.nrow; i++){
    for(k=S.row_ptr[i]; k<S.row_ptr[i+1]; k++){
      j = S.col_ind[k];//a_ij=S.val[k]
      St.row_ptr[j+1]++;//first, count the number of nonzeros in each column of A
    }
  }
  //then, add the total number of nonzeros before i-th column
  for(j=0; j<St.nrow; j++) St.row_ptr[j+1] += St.row_ptr[j];
  
  // construct an array of size n to record current available position in each column of A
  pos = alloc1int(St.nrow);
  memset(pos, 0, St.nrow*sizeof(int));
  for(i=0; i<S.nrow; i++){
    for(k=S.row_ptr[i]; k<S.row_ptr[i+1]; k++){
      j = S.col_ind[k];//a_ij=S.val[k]
      //add number of nonzeros before current column and number of nonzeros in this column
      m = St.row_ptr[j] + pos[j];//row_ptr[j] stores number of nonzeros before j-th row of St
      St.val[m] = S.val[k];//so the m-th element in St.val[] is: a_ij=S.val[k]
      St.col_ind[m] = i;
      pos[j]++;//increase the counter
    }
  }
  free1int(pos);

  return St;
}

void icsr_close(icsr_t *S)
{
  free1int(S->row_ptr);
  free1int(S->col_ind);
  free1int(S->val);
}

csr_t csr_transpose(csr_t A)
{
  int i, j, k, m;
  int *pos;
  csr_t At;

  At.nrow = A.ncol;
  At.ncol = A.nrow;
  At.row_ptr = alloc1int(At.nrow+1);

  memset(At.row_ptr, 0, (At.nrow+1)*sizeof(int));
  for(i=0; i<A.nrow; i++){
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];//a_ij=A.val[k]
      //printf("j=%d ncol=%d\n", j, A.ncol);
      At.row_ptr[j+1]++;//first, count the number of nonzeros in each column of A (each row of At)
    }
  }
  //then, add the total number of nonzeros before i-th column
  for(j=0; j<A.ncol; j++) At.row_ptr[j+1] += At.row_ptr[j];

  At.nnz = At.row_ptr[At.nrow];
  At.val = alloc1double(At.nnz);
  At.col_ind = alloc1int(At.nnz);
  
  // construct an array of size n to record current available position in each column of A
  pos = alloc1int(At.ncol);
  memset(pos, 0, At.ncol*sizeof(int));
  for(i=0; i<A.nrow; i++){
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];//a_ij=A.val[k]
      //add number of nonzeros before current column and number of nonzeros in this column
      m = At.row_ptr[j] + pos[j];//row_ptr[j]=number of nonzeros in previous (j-1)-row of At, (j-1)-column of A
      At.val[m] = A.val[k];//so the m-th element in At.val[] is: a_ij=A.val[k]
      At.col_ind[m] = i;
      pos[j]++;//increase the counter
    }
  }
  free1int(pos);

  return At;
}


void csr_close(csr_t *A)
{
  free1int(A->row_ptr);
  free1int(A->col_ind);
  free1double(A->val);
}
