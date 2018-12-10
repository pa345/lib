/*
 * mfield_error.c
 *
 * Routines for computing uncertainties in model parameters
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

#include "mfield.h"
#include "mfield_error.h"

#include "lapack_wrapper.h"

/*
mfield_calc_uncertainties()
  Print uncertainties for internal field coefficients to a file

Inputs: filename - output file
        covar    - covariance matrix
        w        - workspace
*/

int
mfield_print_uncertainties(const char * filename, const gsl_matrix * covar, mfield_workspace * w)
{
  int s = GSL_SUCCESS;
  gsl_vector_const_view d = gsl_matrix_const_diagonal(covar);
  FILE *fp;
  size_t n;

  fp = fopen(filename, "w");

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic order m\n", n++);
  fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
  fprintf(fp, "# Field %zu: g(n,m) (nT)\n", n++);
  fprintf(fp, "# Field %zu: g(n,m) error (nT)\n", n++);
  fprintf(fp, "# Field %zu: uncertainty in MF g(n,m) (dimensionless)\n", n++);
  fprintf(fp, "# Field %zu: d/dt g(n,m) (nT)\n", n++);
  fprintf(fp, "# Field %zu: d/dt g(n,m) error (nT)\n", n++);
  fprintf(fp, "# Field %zu: uncertainty in SV g(n,m) (dimensionless)\n", n++);
  fprintf(fp, "# Field %zu: (d/dt)^2 g(n,m) (nT)\n", n++);
  fprintf(fp, "# Field %zu: (d/dt)^2 g(n,m) error (nT)\n", n++);
  fprintf(fp, "# Field %zu: uncertainty in SA g(n,m) (dimensionless)\n", n++);

  for (n = 1; n <= w->nmax; ++n)
    {
      int m, ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm = mfield_get_mf(w->c, cidx, w);
          double dgnm = mfield_get_sv(w->c, cidx, w);
          double ddgnm = mfield_get_sa(w->c, cidx, w);
          double err_gnm = mfield_get_mf(&d.vector, cidx, w);
          double err_dgnm = mfield_get_sv(&d.vector, cidx, w);
          double err_ddgnm = mfield_get_sa(&d.vector, cidx, w);

          fprintf(fp, "%5d %5zu %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                  m,
                  n,
                  gnm,
                  err_gnm,
                  fabs(err_gnm / gnm),
                  dgnm,
                  err_dgnm,
                  fabs(err_dgnm / dgnm),
                  ddgnm,
                  err_ddgnm,
                  fabs(err_ddgnm / ddgnm));
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/*
mfield_covariance()
  Compute covariance matrix:

C = [J^T W J]^{-1}

and

W = diag(w_i) = diag(1 / sigma^2)  with
*/

int
mfield_covariance(gsl_matrix * covar, mfield_workspace *w)
{
  int s = GSL_SUCCESS;

  if (w->lls_solution == 1)
    {
      /*
       * for linear systems, mfield_calc_nonlinear_multilarge() stores the
       * Cholesky factor in w->choleskyL
       */

      gsl_matrix_tricpy('L', 1, covar, w->choleskyL);

      /* invert (J^T J) matrix using Cholesky factor */
      lapack_cholesky_invert(covar);

      /* copy lower to upper triangle */
      gsl_matrix_transpose_tricpy('L', 0, covar, covar);
    }
  else
    {
      s = gsl_multilarge_nlinear_covar(covar, w->nlinear_workspace_p);
    }

  return s;
}
