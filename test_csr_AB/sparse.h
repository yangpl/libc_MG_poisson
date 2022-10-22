#ifndef __sparse_h__
#define __sparse_h__


typedef struct{
  int nrow;
  int ncol;
  int nnz;//number of non-zeros
  int *row_ptr;//row pointer
  int *col_ind;//column indices
  double *val;//values of the matrix A
} csr_t;//compressed sparse row format (compressed row storage, CRS)


#endif
