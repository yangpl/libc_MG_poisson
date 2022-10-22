#include "sxamg.h"
#include <assert.h>

SX_MAT laplacian_7pt_bound(const SX_INT Nx, const SX_INT Ny, const SX_INT Nz)
{
  SX_MAT A;
  SX_INT nz = 0;
  SX_INT i, j, k;

  A.num_rows = Nx * Ny * Nz;
  A.num_cols = Nx * Ny * Nz;
  A.num_nnzs = 7 * Nx * Ny * Nz;

  A.Ap = sx_malloc(sizeof(SX_INT) * (A.num_rows + 1));
  A.Aj = sx_malloc(sizeof(SX_INT) * A.num_nnzs);
  A.Ax = sx_malloc(sizeof(SX_FLT) * A.num_nnzs);

  A.Ap[0] = 0;
  for (i = 0; i < A.num_rows; i++) A.Ap[i + 1] = 0;

  for (i = 0; i < Nz; i++) {
    for (j = 0; j < Ny; j++) {
      for (k = 0; k < Nx; k++) {
	SX_INT indx = Nx * Ny * i + Nx * j + k;

	if (i > 0) {
	  A.Aj[nz] = indx - Nx * Ny;
	  A.Ax[nz] = -1;
	  nz++;
	}

	if (j > 0) {
	  A.Aj[nz] = indx - Nx;
	  A.Ax[nz] = -1;
	  nz++;
	}

	if (k > 0) {
	  A.Aj[nz] = indx - 1;
	  A.Ax[nz] = -1;
	  nz++;
	}

	A.Aj[nz] = indx;
	A.Ax[nz] = 6;
	nz++;

	if (k < Nx - 1) {
	  A.Aj[nz] = indx + 1;
	  A.Ax[nz] = -1;
	  nz++;
	}

	if (j < Ny - 1) {
	  A.Aj[nz] = indx + Nx;
	  A.Ax[nz] = -1;
	  nz++;
	}

	if (i < Nz - 1) {
	  A.Aj[nz] = indx + Nx * Ny;
	  A.Ax[nz] = -1;
	  nz++;
	}

	A.Ap[indx + 1] = nz;
      }
    }
  }

  A.num_nnzs = nz;
  return A;
}

