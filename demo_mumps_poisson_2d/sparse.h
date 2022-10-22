#ifndef __sparse_h__
#define __sparse_h__

typedef struct{
  int nrow, ncol;//matrix size m-rows, n-columns
  int nnz;//number of non-zeros
  int *row;//row indices
  int *col;//column indices
  double *val;//values of the matrix A
} coo_t;//COOrdinate format



typedef struct{
  int nrow;
  int ncol;
  int nnz;
  int *row_ptr;
  int *col_ind;
  double *val;
}csr_t;

typedef struct{
  int nrow;
  int ncol;
  int nnz;
  int *row_ptr;
  int *col_ind;
  int *val;
}icsr_t;//csr for integer array

#endif
