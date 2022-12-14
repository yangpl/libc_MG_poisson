#include "sxamg.h"
#include <assert.h>
SX_MAT laplacian_7pt_bound(const SX_INT Nx, const SX_INT Ny, const SX_INT Nz);

int main(void)
{
  SX_MAT A;
  SX_INT ncx = 123, ncy = 133, ncz = 124;
  SX_INT k, i;
  SX_VEC x, b;

  assert(ncx > 0);
  assert(ncy > 0);
  assert(ncz > 0);

  /* create distributed matrix */
  k = ncx * ncy * ncz;
  A = laplacian_7pt_bound(ncx, ncy, ncz);
  sx_printf("sx: problem size: %d, %d x %d x %d.\n", k, ncx, ncy, ncz);
  sx_printf("A: m = %"dFMT", n = %"dFMT", nnz = %"dFMT"\n", A.num_rows,
            A.num_cols, A.num_nnzs);

  /* create vectors */
  x = sx_vec_create(k);
  b = sx_vec_create(k);
  /* init x and b */
  for (i = 0; i < k; i++) {
    sx_vec_set_entry(&b, i, 1);
    sx_vec_set_entry(&x, i, 0);
  }

  /* solve Ax = b using CG with AMG preconditioner */
  {
    SX_INT maxits = 2000;  /* maximal iterations */
    SX_FLT tol = 1e-6;     /* stop tolerance */
    SX_FLT err, err0;
    SX_VEC r;     /* residual */
    SX_VEC p;
    SX_VEC z, q;
    SX_FLT rho1, rho0 = 0., alpha, beta;
    SX_FLT stm;

    /* preconditioner */
    SX_AMG mg;
    SX_AMG_PARS pars;

    /* pars */
    sx_amg_pars_init(&pars);
    pars.maxit = 1;
    sx_amg_setup(&mg, &A, &pars);

    /* create vector */
    r = sx_vec_create(k);
    z = sx_vec_create(k);
    q = sx_vec_create(k);
    p = sx_vec_create(k);

    /* initial residual */
    sx_blas_mv_amxpbyz(-1., &A, &x, 1., &b, &r);
    err0 = sx_blas_vec_norm2(&r);

    sx_printf("\nsx: solver: CG, preconditioner: AMG\n");
    sx_printf("Convergence settings: relative residual: %"fFMTe
	      ", maximal iterations: %"dFMT"\n\n", tol, maxits);
    sx_printf("Initial residual: %"fFMTe"\n\n", err0);

    /* timer */
    stm = sx_get_time();

    for (i = 0; i < maxits; i++) {
      /* supposed to solve preconditioning system Mz = r */
#if 0       /* no pc */
      sx_blas_vec_copy(&r, &z);
#else       /* amg pc */
      sx_blas_vec_set(&z, 0.);
      sx_solver_amg_solve(&mg, &z, &r);
#endif

      /* rho = <r, z> */
      rho1 = sx_blas_vec_dot(&r, &z);

      if (i == 0) {
	/* p = z */
	sx_blas_vec_copy(&z, &p);
      }
      else {
	beta = rho1 / rho0;

	/* update p */
	sx_blas_vec_axpby(1, &z, beta, &p);
      }

      /* save rho */
      rho0 = rho1;

      /* update q */
      sx_blas_mv_mxy(&A, &p, &q);

      /* compute alpha */
      alpha = rho1 / sx_blas_vec_dot(&p, &q);

      /* update x */
      sx_blas_vec_axpy(alpha, &p, &x);

      /* update r */
      sx_blas_vec_axpy(-alpha, &q, &r);

      /* check convergence */
      err = sx_blas_vec_norm2(&r);

      sx_printf("itr: %6"dFMT",     residual: %"fFMTe", relative error: %"fFMTe"\n",
		i + 1, err, err / err0);

      if (err / err0 <= tol) break;
    }

    if (i == maxits) {
      sx_printf("\ncg failed to converge.\n\n");
    }
    else {
      stm = sx_get_time() - stm;

      sx_printf("\ncg converged: absolute residual: %"fFMTe
		", total iterations: %"dFMT", time: %g s\n\n", err, i + 1, stm);
    }

    sx_vec_destroy(&r);
    sx_vec_destroy(&p);
    sx_vec_destroy(&z);
    sx_vec_destroy(&q);

    sx_amg_data_destroy(&mg);
  }

  /* release memory */
  sx_mat_destroy(&A);
  sx_vec_destroy(&x);
  sx_vec_destroy(&b);

  return 0;
}
