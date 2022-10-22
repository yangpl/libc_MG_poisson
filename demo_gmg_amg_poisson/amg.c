#include "cstd.h"
#include "sparse.h"
#include "sxamg.h"


/* preconditioner */
SX_AMG amg;
SX_AMG_PARS pars;

void amg_init(csr_t *A)
{
  SX_MAT A_;

  A_.num_rows = A->nrow;
  A_.num_cols = A->ncol;
  A_.num_nnzs = A->nnz;
  A_.Ap = A->row_ptr;
  A_.Aj = A->col_ind;
  A_.Ax = A->val;
  
  sx_amg_pars_init(&pars);
  pars.maxit = 1;
  sx_amg_setup(&amg, &A_, &pars);
}

void amg_close()
{
  sx_amg_data_destroy(&amg);
}

void amg_apply(int n, double *x, double *y)//preconditioner, y=invA*x
{

  SX_VEC z, r;
  z.n = n;
  z.d = y;
  r.n = n;
  r.d = x;
  sx_blas_vec_set(&z, 0.);
  sx_solver_amg_solve(&amg, &z, &r);

}
