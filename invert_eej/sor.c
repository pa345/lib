/* sor.c
 * 
 * Copyright (C) 2018 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>

/* this module contains routines related to the successive over relaxation method */

/*
pde_sor()
  Perform one iteration of the successive over relaxation method for the square
system A x = b

Inputs: omega - relaxation parameter
        tol   - convergence tolerance
        A     - sparse matrix
        b     - right hand side vector
        x     - (input/output) on input, initial guess for solution
                               on output, final solution

Return: GSL_SUCCESS if converged (||b - A x|| <= tol * ||b||)
        GSL_CONTINUE if not converged
        error code for input parameter error
*/

int
pde_sor(const double omega, const double tol, const gsl_spmatrix * A, const gsl_vector * b, gsl_vector * x)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR("matrix must be square", GSL_ENOTSQR);
    }
  else if (b->size != N)
    {
      GSL_ERROR("right hand side vector does not match matrix", GSL_EBADLEN);
    }
  else if (x->size != N)
    {
      GSL_ERROR("solution vector does not match matrix", GSL_EBADLEN);
    }
  else if (omega < 0.0 || omega > 2.0)
    {
      GSL_ERROR("relaxation parameter omega must be between 0 and 2", GSL_EDOM);
    }
  else if (A->sptype != GSL_SPMATRIX_CRS)
    {
      GSL_ERROR("SOR currently supports only compressed row format", GSL_EINVAL);
    }
  else
    {
      int status;
      const double normb = gsl_blas_dnrm2(b);
      const size_t *p = A->p;
      const size_t *colidx = A->i;
      const double *data = A->data;
      gsl_vector *r = gsl_vector_alloc(N);
      double normr;
      size_t i;

      for (i = 0; i < N; ++i)
        {
          int diag_idx = -1; /* index of a_{ii} */
          double sigma = 0.0;
          double *xi = gsl_vector_ptr(x, i);
          double bi = gsl_vector_get(b, i);
          size_t k;

          /* loop over elements of row i */
          for (k = p[i]; k < p[i + 1]; ++k)
            {
              if (i == colidx[k])
                diag_idx = (int) k;
              else
                sigma += data[k] * gsl_vector_get(x, colidx[k]);
            }

          if (diag_idx < 0 || data[diag_idx] == 0.0)
            {
              GSL_ERROR("diagonal element is zero", GSL_EDOM);
            }

          *xi += omega * ((bi - sigma) / data[diag_idx] - (*xi));
        }

      /* convergence test: || b - A x || <= tol * || b || */

      gsl_vector_memcpy(r, b);
      gsl_spblas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);
      normr = gsl_blas_dnrm2(r);

      if (normr <= tol * normb)
        status = GSL_SUCCESS;
      else
        status = GSL_CONTINUE;

      gsl_vector_free(r);

      return status;
    }
}
