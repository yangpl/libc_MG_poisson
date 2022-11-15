#include "cstd.h"
#include "sparse.h"

csr_t csr_transpose(csr_t A)
{

  int i, j, k, m;
  int *pos;
  csr_t At;

  At.nrow = A.ncol;
  At.ncol = A.nrow;
  At.nnz = A.nnz;
  At.row_ptr = alloc1int(At.nrow+1);
  At.col_ind = alloc1int(At.nnz);
  At.val = alloc1double(At.nnz);

  memset(At.row_ptr, 0, (At.nrow+1)*sizeof(int));
  for(i=0; i<A.nrow; i++){
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];//a_ij=A.val[k]
      At.row_ptr[j+1]++;//first, count the number of nonzeros in each column of A
    }
  }
  //then, add the total number of nonzeros before i-th column
  //printf("%d\n", At.row_ptr[0]);
  for(j=0; j<At.nrow; j++){
    At.row_ptr[j+1] += At.row_ptr[j];
    //printf("%d\n", At.row_ptr[j+1]);
  }

  
  // construct an array of size n to record current available position in each column of A
  pos = alloc1int(At.nrow);
  memset(pos, 0, At.nrow*sizeof(int));
  for(i=0; i<A.nrow; i++){
    for(k=A.row_ptr[i]; k<A.row_ptr[i+1]; k++){
      j = A.col_ind[k];//a_ij=A.val[k]
      //add number of nonzeros before current column and number of nonzeros in this column
      m = At.row_ptr[j] + pos[j];//row_ptr[j] stores number of nonzeros before j-th row of At
      At.val[m] = A.val[k];//so the m-th element in At.val[] is: a_ij=A.val[k]
      At.col_ind[m] = i;
      pos[j]++;//increase the counter
    }
  }
  free1int(pos);

  return At;
}

//C_ik=sum_j A_ij*B_jk, CSR sparse matrix multiplication
csr_t build_AB(csr_t A, csr_t B)
{
  int i, j, k, k1, k2, m;
  csr_t C;
  int *pos;

  C.nrow = A.nrow;
  C.ncol = B.ncol;
  C.row_ptr = alloc1int(C.nrow+1);

  pos = alloc1int(C.ncol);
  for(i=0; i<C.ncol; i++) pos[i] = -1;
  m = 0;
  for(i=0; i<A.nrow; i++){
    C.row_ptr[i] = m;
    for(k1=A.row_ptr[i]; k1<A.row_ptr[i+1]; k1++){
      j = A.col_ind[k1];
      for(k2=B.row_ptr[j]; k2<B.row_ptr[j+1]; k2++){
	k = B.col_ind[k2];
	if(pos[k]<C.row_ptr[i]){//it is a new column because the pos index is not up to date
	  printf("k=%d pos[n]=%d \t", k, pos[k]);
	  pos[k] = m;//T_mn is the k-th element of matrix T
	  printf("after: pos[n]=%d\n", pos[k]);
	  m++;//increase global element counter for matrix T
	}
      }
    }
  }
  C.row_ptr[C.nrow] = m;
  C.nnz = m;
  printf("nnz=%d\n", m);
  C.col_ind = alloc1int(C.nnz);
  C.val = alloc1double(C.nnz);
  memset(C.val, 0, C.nnz*sizeof(double));

  for(i=0; i<C.ncol; i++) pos[i] = -1;
  m = 0;
  for(i=0; i<A.nrow; i++){
    C.row_ptr[i] = m;
    for(k1=A.row_ptr[i]; k1<A.row_ptr[i+1]; k1++){
      j = A.col_ind[k1];//a_ij=A.val[k1]
      for(k2=B.row_ptr[j]; k2<B.row_ptr[j+1]; k2++){
	k = B.col_ind[k2];//b_jk=B.val[k2]
	if(pos[k]<C.row_ptr[i]){
	  //this is a new column because the pos index is not up to date
	  pos[k] = m;//T_mn is the k-th element of matrix T
	  C.col_ind[m] = k;
	  m++;//increase global element counter for matrix T
	}
	C.val[pos[k]] += A.val[k1]*B.val[k2];
      }
    }
  }  

  free1int(pos);

  return C;
}


int main()
{
  int i, j, k;
  int n, B[4][4];
  csr_t A, At, C;

  n = 4;
  A.nrow = n;
  A.ncol = n;
  A.nnz = n*n;
  A.row_ptr = alloc1int(A.nrow+1);
  A.col_ind = alloc1int(A.nnz);
  A.val = alloc1double(A.nnz);

  B[0][0] = 1;
  B[0][1] = 2;
  B[0][2] = 3;
  B[0][3] = 4;
  B[1][0] = 5;
  B[1][1] = 6;
  B[1][2] = 7;
  B[1][3] = 8;
  B[2][0] = 9;
  B[2][1] = 10;
  B[2][2] = 11;
  B[2][3] = 12;
  B[3][0] = 13;
  B[3][1] = 14;
  B[3][2] = 15;
  B[3][3] = 16;

  A.row_ptr[0] = 0;
  k = 0;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      A.col_ind[k] = j;
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

  
  At = csr_transpose(A);
  printf("-------At---------\n");
  for(i=0; i<At.nrow; i++){
    for(k=At.row_ptr[i]; k<At.row_ptr[i+1]; k++){
      j = At.col_ind[k];
      printf("%g \t", At.val[k]);
    }
    printf("\n");
  }

  C = build_AB(A, A);
  printf("----------C=A*A------\n");
  for(i=0; i<C.nrow; i++){
    for(k=C.row_ptr[i]; k<C.row_ptr[i+1]; k++){
      j = C.col_ind[k];
      printf("%g \t", C.val[k]);
    }
    printf("\n");
  }
  

}
