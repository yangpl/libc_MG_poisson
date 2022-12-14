
#include "sxamg.h"
SX_MAT laplacian_7pt_bound(const SX_INT Nx, const SX_INT Ny, const SX_INT Nz);

int main(void)
{
    SX_AMG_PARS pars;
    SX_MAT A;
    SX_VEC b, x;
    SX_AMG mg;
    SX_RTN rtn;
    SX_INT ncx = 23, ncy = 13, ncz = 24;
    SX_INT nglobal = 0;

    /* create distributed matrix */
    nglobal = ncx * ncy * ncz;
    A = laplacian_7pt_bound(ncx, ncy, ncz);
    sx_printf("sx: problem size: %d, %d x %d x %d.\n", nglobal, ncx, ncy, ncz);

    /* pars */
    sx_amg_pars_init(&pars);
    pars.maxit = 1000;
    pars.verb = 2;
    
    /* print info */
    sx_printf("\nA: m = %"dFMT", n = %"dFMT", nnz = %"dFMT"\n", A.num_rows,
            A.num_cols, A.num_nnzs);

    sx_amg_pars_print(&pars);

    // Step 1: AMG setup phase
    sx_amg_setup(&mg, &A, &pars);

    // Step 2: AMG solve phase
    b = sx_vec_create(A.num_rows);
    sx_vec_set_value(&b, 1.0);
    
    x = sx_vec_create(A.num_rows);
    sx_vec_set_value(&x, 1.0);

    /* solve */
    rtn = sx_solver_amg_solve(&mg, &x, &b);

    sx_printf("AMG residual: %"fFMTg"\n", rtn.ares);
    sx_printf("AMG relative residual: %"fFMTg"\n", rtn.rres);
    sx_printf("AMG iterations: %"dFMT"\n", rtn.nits);

    sx_mat_destroy(&A);
    sx_vec_destroy(&x);
    sx_vec_destroy(&b);
    sx_amg_data_destroy(&mg);
    
    return 0;
}
