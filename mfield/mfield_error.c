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
#include <gsl/gsl_linalg.h>
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
#if 0 /*XXX*/
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
#endif

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

      gsl_matrix_tricpy(CblasLower, CblasNonUnit, covar, w->choleskyL);
      s += gsl_linalg_cholesky_invert(covar);
    }
  else
    {
      s = gsl_multilarge_nlinear_covar(covar, w->nlinear_workspace_p);
    }

  /* copy lower to upper triangle */
  gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, covar, covar);

  return s;
}

/*
mfield_correlation()
  The correlation matrix is computed from the covariance matrix:

Corr = D^{-1/2} Covar D^{-1/2}

where D = diag(diag(Covar))

Inputs: C - on input, covariance matrix
            on output, correlation matrix
*/

int
mfield_correlation(gsl_matrix * C, mfield_workspace *w)
{
  int s = 0;
  const size_t N = C->size1;
  gsl_vector_const_view diag = gsl_matrix_const_diagonal(C);
  gsl_vector * d = gsl_vector_alloc(N);
  size_t i;

  gsl_vector_memcpy(d, &diag.vector);

  /* compute D^{-1/2} C */
  for (i = 0; i < N; ++i)
    {
      double di = gsl_vector_get(d, i);
      gsl_vector_view v = gsl_matrix_row(C, i);
      gsl_blas_dscal(1.0 / sqrt(di), &v.vector);
    }

  /* compute C D^{-1/2} */
  for (i = 0; i < N; ++i)
    {
      double di = gsl_vector_get(d, i);
      gsl_vector_view v = gsl_matrix_column(C, i);
      gsl_blas_dscal(1.0 / sqrt(di), &v.vector);
    }

  gsl_vector_free(d);

  return s;
}
