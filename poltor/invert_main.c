/*
 * invert_main.c
 *
 * Invert vector magnetic field measurements for toroidal currents
 * flowing in a thin shell at 110km altitude
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <libconfig.h>

#include <flow/flow.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_track.h>
#include <mainlib/ml_magdata_list.h>

#include "green.h"
#include "lapack_wrapper.h"

#include "poltor.h"

/* define to use QD coordinates in spherical harmonic transforms */
#define POLTOR_QD_HARMONICS        1

static int
parse_config_file(const char *filename, poltor_parameters *poltor_params)
{
  int s;
  config_t cfg;
  double fval;
  int ival;

  config_init(&cfg);

  s = config_read_file(&cfg, filename);
  if (s != CONFIG_TRUE)
    {
      fprintf(stderr, "parse_config_file: %s:%d - %s\n",
              config_error_file(&cfg),
              config_error_line(&cfg),
              config_error_text(&cfg));
      config_destroy(&cfg);
      return -1;
    }

  if (config_lookup_int(&cfg, "nmax_int", &ival))
    poltor_params->nmax_int = (size_t) ival;
  if (config_lookup_int(&cfg, "mmax_int", &ival))
    poltor_params->mmax_int = (size_t) ival;
  if (config_lookup_int(&cfg, "nmax_ext", &ival))
    poltor_params->nmax_ext = (size_t) ival;
  if (config_lookup_int(&cfg, "mmax_ext", &ival))
    poltor_params->mmax_ext = (size_t) ival;
  if (config_lookup_int(&cfg, "nmax_sh", &ival))
    poltor_params->nmax_sh = (size_t) ival;
  if (config_lookup_int(&cfg, "mmax_sh", &ival))
    poltor_params->mmax_sh = (size_t) ival;
  if (config_lookup_int(&cfg, "nmax_tor", &ival))
    poltor_params->nmax_tor = (size_t) ival;
  if (config_lookup_int(&cfg, "mmax_tor", &ival))
    poltor_params->mmax_tor = (size_t) ival;

  if (config_lookup_float(&cfg, "R", &fval))
    poltor_params->R = fval;
  if (config_lookup_float(&cfg, "b", &fval))
    poltor_params->b = fval;

  if (config_lookup_int(&cfg, "shell_J", &ival))
    poltor_params->shell_J = (size_t) ival;

  if (config_lookup_int(&cfg, "max_iter", &ival))
    poltor_params->max_iter = (size_t) ival;

  if (config_lookup_int(&cfg, "use_weights", &ival))
    poltor_params->use_weights = ival;
  if (config_lookup_int(&cfg, "regularize", &ival))
    poltor_params->regularize = ival;
  if (config_lookup_int(&cfg, "solar_flux_factor", &ival))
    poltor_params->solar_flux_factor = ival;

  if (config_lookup_float(&cfg, "alpha_int", &fval))
    poltor_params->alpha_int = fval;
  if (config_lookup_float(&cfg, "alpha_sh", &fval))
    poltor_params->alpha_sh = fval;
  if (config_lookup_float(&cfg, "alpha_tor", &fval))
    poltor_params->alpha_tor = fval;

  if (config_lookup_float(&cfg, "weight_X", &fval))
    poltor_params->weight_X = fval;
  if (config_lookup_float(&cfg, "weight_Y", &fval))
    poltor_params->weight_Y = fval;
  if (config_lookup_float(&cfg, "weight_Z", &fval))
    poltor_params->weight_Z = fval;
  if (config_lookup_float(&cfg, "weight_F", &fval))
    poltor_params->weight_F = fval;
  if (config_lookup_float(&cfg, "weight_DX_NS", &fval))
    poltor_params->weight_DX_NS = fval;
  if (config_lookup_float(&cfg, "weight_DY_NS", &fval))
    poltor_params->weight_DY_NS = fval;
  if (config_lookup_float(&cfg, "weight_DZ_NS", &fval))
    poltor_params->weight_DZ_NS = fval;
  if (config_lookup_float(&cfg, "weight_DF_NS", &fval))
    poltor_params->weight_DF_NS = fval;

  if (config_lookup_int(&cfg, "fit_X", &ival))
    poltor_params->fit_X = ival;
  if (config_lookup_int(&cfg, "fit_Y", &ival))
    poltor_params->fit_Y = ival;
  if (config_lookup_int(&cfg, "fit_Z", &ival))
    poltor_params->fit_Z = ival;
  if (config_lookup_int(&cfg, "fit_F", &ival))
    poltor_params->fit_F = ival;
  if (config_lookup_int(&cfg, "fit_DX_NS", &ival))
    poltor_params->fit_DX_NS = ival;
  if (config_lookup_int(&cfg, "fit_DY_NS", &ival))
    poltor_params->fit_DY_NS = ival;
  if (config_lookup_int(&cfg, "fit_DZ_NS", &ival))
    poltor_params->fit_DZ_NS = ival;
  if (config_lookup_int(&cfg, "fit_DF_NS", &ival))
    poltor_params->fit_DF_NS = ival;

#if 0
  if (config_lookup_int(&cfg, "fit_DX_EW", &ival))
    data_params->fit_DX_EW = ival;
  if (config_lookup_int(&cfg, "fit_DY_EW", &ival))
    data_params->fit_DY_EW = ival;
  if (config_lookup_int(&cfg, "fit_DZ_EW", &ival))
    data_params->fit_DZ_EW = ival;
  if (config_lookup_int(&cfg, "fit_DF_EW", &ival))
    data_params->fit_DF_EW = ival;

  if (config_lookup_int(&cfg, "fit_Z_highlat", &ival))
    data_params->fit_Z_highlat = ival;
  if (config_lookup_int(&cfg, "fit_F_highlat", &ival))
    data_params->fit_F_highlat = ival;
  if (config_lookup_int(&cfg, "fit_DZ_NS_highlat", &ival))
    data_params->fit_DZ_NS_highlat = ival;
  if (config_lookup_int(&cfg, "fit_DF_NS_highlat", &ival))
    data_params->fit_DF_NS_highlat = ival;
  if (config_lookup_int(&cfg, "fit_DZ_EW_highlat", &ival))
    data_params->fit_DZ_EW_highlat = ival;
  if (config_lookup_int(&cfg, "fit_DF_EW_highlat", &ival))
    data_params->fit_DF_EW_highlat = ival;
#endif

  if (config_lookup_int(&cfg, "synth_data", &ival))
    poltor_params->synth_data = ival;

  config_destroy(&cfg);

  return 0;
}

int
initial_guess(gsl_vector_complex *c, poltor_workspace *w)
{
  int s = 0;
  const poltor_parameters *params = &(w->params);

  if (params->synth_data)
    {
#if 0
      poltor_synth_init(c, w);
#else
      /* set the (1,0) element to 1, since a 0 vector causes problems with the scalar Jacobian,
       * which divides by || B_model || */
      gsl_vector_complex_set_zero(c);
      gsl_vector_complex_set(c, 1, GSL_COMPLEX_ONE);
#endif
    }
  else
    {
      /* set the (1,0) element to 1, since a 0 vector causes problems with the scalar Jacobian,
       * which divides by || B_model || */
      gsl_vector_complex_set_zero(c);
      gsl_vector_complex_set(c, 1, GSL_COMPLEX_ONE);
    }

  return s;
}

/*
set_flags()
  The preproc program will have calculated the magdata struct,
but won't set the individual X/Y/Z fitting flags, since that should
be done here. However it would set the MAGDATA_FLG_DZ_NS flag if
gradient information is available for each point.

So go through data, and based on the fit settings, re-assign
magdata flags
*/

int
set_flags(const poltor_parameters *params, magdata_list *list)
{
  int s = 0;
  size_t i, j;

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* don't fit X/Y data at high latitudes */
          if (fabs(mptr->qdlat[j]) > 60.0)
            {
              mptr->flags[j] &= ~(MAGDATA_FLG_X|MAGDATA_FLG_Y);
              mptr->flags[j] &= ~(MAGDATA_FLG_DX_NS|MAGDATA_FLG_DY_NS);
            }

          if (!params->fit_X)
            mptr->flags[j] &= ~MAGDATA_FLG_X;

          if (!params->fit_Y)
            mptr->flags[j] &= ~MAGDATA_FLG_Y;

          if (!params->fit_Z)
            mptr->flags[j] &= ~MAGDATA_FLG_Z;

          if (!params->fit_F)
            mptr->flags[j] &= ~MAGDATA_FLG_F;

          if (!params->fit_DX_NS)
            mptr->flags[j] &= ~MAGDATA_FLG_DX_NS;

          if (!params->fit_DY_NS)
            mptr->flags[j] &= ~MAGDATA_FLG_DY_NS;

          if (!params->fit_DZ_NS)
            mptr->flags[j] &= ~MAGDATA_FLG_DZ_NS;

          if (!params->fit_DF_NS)
            mptr->flags[j] &= ~MAGDATA_FLG_DF_NS;

          if (!(mptr->flags[j] & (MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z | MAGDATA_FLG_F |
                                  MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DF_NS |
                                  MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW | MAGDATA_FLG_DF_EW)))
            {
              /* this point is no longer used in the fitting */
              mptr->flags[j] |= MAGDATA_FLG_DISCARD;
            }
        }
    }

  return s;
}

/*
replace_qdlat()
  Replace theta with theta_q
*/

int
replace_qdlat(magdata_list *list)
{
  int s = 0;
  size_t i, j;

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          mptr->theta[j] = M_PI / 2.0 - mptr->qdlat[j] * M_PI / 180.0;
          mptr->theta_ns[j] = M_PI / 2.0 - mptr->qdlat_ns[j] * M_PI / 180.0;
        }
    }

  return s;
}

int
print_error_source(FILE *fp, const size_t nmax, const size_t mmax,
                   const size_t offset, gsl_vector_complex *err,
                   gsl_vector_complex *c)
{
  int s = 0;
  size_t n;

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic order m\n", n++);
  fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
  fprintf(fp, "# Field %zu: uncertainty in c(n,m) (nT)\n", n++);
  fprintf(fp, "# Field %zu: Re[c(n,m)] (nT)\n", n++);
  fprintf(fp, "# Field %zu: Im[c(n,m)] (nT)\n", n++);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_idx(n, m, mmax) + offset;
          gsl_complex cnm = gsl_vector_complex_get(c, cidx);
          gsl_complex err_cnm = gsl_vector_complex_get(err, cidx);

          fprintf(fp, "%5d %5zu %20.4e %20.4e %20.4e\n", m, n, GSL_REAL(err_cnm), GSL_REAL(cnm), GSL_IMAG(cnm));
        }

      fprintf(fp, "\n");
    }

  fprintf(fp, "\n\n");

  return s;
}

/* w->JHJ must contain covariance matrix on entry */
int
print_uncertanties(const char *filename, poltor_workspace *w)
{
  int s = 0;
  gsl_matrix_complex *JHJ = w->JHJ;
  gsl_vector_complex_view diag = gsl_matrix_complex_diagonal(JHJ);
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_uncertanties: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  if (w->p_pint > 0)
    print_error_source(fp, w->nmax_int, w->mmax_int, w->pint_offset, &diag.vector, w->c);

  if (w->p_pext > 0)
    print_error_source(fp, w->nmax_ext, w->mmax_ext, w->pext_offset, &diag.vector, w->c);

  if (w->p_psh > 0)
    print_error_source(fp, w->nmax_sh, w->mmax_sh, w->psh_offset, &diag.vector, w->c);

  if (w->p_tor > 0)
    print_error_source(fp, w->nmax_tor, w->mmax_tor, w->tor_offset, &diag.vector, w->c);

  fclose(fp);

  return s;
}

int
print_correlation(const char *filename, poltor_workspace *w)
{
  int s = 0;
  const size_t p = w->p;
  const size_t nmax = GSL_MIN(w->nmax_int, w->nmax_sh);
  const size_t mmax = GSL_MIN(w->mmax_int, w->mmax_sh);
  gsl_matrix_complex *B;
  gsl_vector_complex_view d;
  size_t n;
  size_t i, j;

  B = gsl_matrix_complex_calloc(p, p);
  gsl_matrix_complex_memcpy(B, w->JHJ);
  d = gsl_matrix_complex_diagonal(B);

  /* covariance matrix is in upper triangle, so fill in lower half of B */
  for (j = 1; j < p; ++j)
    {
      for (i = 0; i < j; ++i)
        {
          gsl_complex z = gsl_matrix_complex_get(B, i, j);
          GSL_SET_IMAG(&z, -GSL_IMAG(z));
          gsl_matrix_complex_set(B, j, i, z); 
        }
    }

  /* compute correlation matrix: B = diag(C)^{-1/2} C diag(C)^{-1/2} */
  for (i = 0; i < n; ++i)
    {
      gsl_complex di = gsl_vector_complex_get(&d.vector, i);
      gsl_vector_complex_view ri = gsl_matrix_complex_row(B, i);
      gsl_vector_complex_view ci = gsl_matrix_complex_column(B, i);
      gsl_complex z;

      GSL_SET_COMPLEX(&z, 1.0 / sqrt(GSL_REAL(di)), 0.0);
      gsl_vector_complex_scale(&ri.vector, z);
      gsl_vector_complex_scale(&ci.vector, z);
    }

  printherm_octave(B, filename);

#if 0
  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax);
      int m;

      for (m = -ni; m <= ni; ++m)
        {
          size_t idx1 = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          size_t idx2 = poltor_jnmidx(0, n, m, w);
          gsl_complex z = gsl_matrix_complex_get(B, idx1, idx2);

          fprintf(fp, "%3zu %3d %e\n",
                  n,
                  m,
                  gsl_complex_abs(z));
        }

      fprintf(fp, "\n\n");
    }
#endif

  gsl_matrix_complex_free(B);

  return s;
} /* print_correlation() */

int
print_coefficients(poltor_workspace *w)
{
  const poltor_parameters *params = &(w->params);

  if (params->synth_data)
    {
      /* print and check synthetic coefficients */
      poltor_synth_print(w);
    }
  else
    {
      const double b = w->b;
      const double R = w->R;
      const size_t nmax_int = GSL_MIN(3, w->nmax_int);
      const size_t mmax_int = GSL_MIN(3, w->mmax_int);
      const size_t nmax_ext = GSL_MIN(2, w->nmax_ext);
      const size_t mmax_ext = GSL_MIN(2, w->mmax_ext);
      const size_t nmax_sh = GSL_MIN(2, w->nmax_sh);
      const size_t mmax_sh = GSL_MIN(2, w->mmax_sh);
      const size_t nmax_tor = GSL_MIN(3, w->nmax_tor);
      const size_t mmax_tor = GSL_MIN(3, w->mmax_tor);
      size_t n;

      /* print internal poloidal coefficients */
      fprintf(stderr, "Internal poloidal coefficients:\n");
      for (n = 1; n <= nmax_int; ++n)
        {
          int ni = (int) GSL_MIN(n, mmax_int);
          int m;

          /*
           * only need to print positive m coefficients since
           * c_{n,-m} = c_{nm}
           */
          for (m = 0; m <= ni; ++m)
            {
              size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
              gsl_complex coef = poltor_get(cidx, w);
              double gnm = GSL_REAL(coef);
              double qnm = -(2.0*n + 1.0) / (double)n * gnm *
                           pow(R / b, n + 3.0);

              fprintf(stderr, "g(%2zu,%2d) = %12g q(%2zu,%2d) = %12g [nT]\n",
                      n, m, gnm,
                      n, m, qnm);
            }
        }

      /* print external poloidal coefficients */
      fprintf(stderr, "External poloidal coefficients:\n");
      for (n = 1; n <= nmax_ext; ++n)
        {
          int ni = (int) GSL_MIN(n, mmax_ext);
          int m;

          /*
           * only need to print positive m coefficients since
           * c_{n,-m} = c_{nm}
           */
          for (m = 0; m <= ni; ++m)
            {
              size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
              gsl_complex coef = poltor_get(cidx, w);
              double knm = GSL_REAL(coef);

              fprintf(stderr, "k(%2zu,%2d) = %12g [nT]\n",
                      n,
                      m,
                      knm);
            }
        }

      /* print shell poloidal coefficients */
      {
        char buf[2048];
        char *bufptr = buf;
        int offset;
        size_t j;

        fprintf(stderr, "Shell poloidal coefficients:\n");
        for (n = 1; n <= nmax_sh; ++n)
          {
            int ni = (int) GSL_MIN(n, mmax_sh);
            int m;

            /*
             * only need to print positive m coefficients since
             * c_{n,-m} = c_{nm}
             */
            for (m = 0; m <= ni; ++m)
              {
                bufptr = buf;
                for (j = 0; j <= w->shell_J; ++j)
                  {
                    size_t cidx = poltor_jnmidx(j, n, m, w);
                    gsl_complex coef = poltor_get(cidx, w);
                    double qnm = GSL_REAL(coef);

                    sprintf(bufptr, "%12g%s %n",
                            qnm,
                            (j < w->shell_J) ? "," : "",
                            &offset);
                    bufptr += offset;
                  }

                fprintf(stderr, "q(%2zu,%2d) = %s [nT]\n", n, m, buf);
              }
          }
      }

      /* print toroidal coefficients */
      fprintf(stderr, "Shell toroidal coefficients:\n");
      for (n = 1; n <= nmax_tor; ++n)
        {
          int ni = (int) GSL_MIN(n, mmax_tor);
          int m;

          /*
           * only need to print positive m coefficients since
           * c_{n,-m} = c_{nm}
           */
          for (m = 0; m <= ni; ++m)
            {
              size_t cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
              gsl_complex coef = poltor_get(cidx, w);
              double phinm = GSL_REAL(coef);

              fprintf(stderr, "phi(%2zu,%2d) = %12g [nT]\n",
                      n,
                      m,
                      phinm);
            }
        }
    }

  return 0;
}

int
print_Lcurve(const char *filename, poltor_workspace *w)
{
  int s = 0;
#if 0 /*XXX*/
  FILE *fp;
  const size_t p = w->p;
  double rnorm, Lnorm;
  gsl_vector_complex_view v = gsl_vector_complex_subvector(w->rhs, 0, p);
  size_t i;

  fp = fopen(filename, "a");
  if (!fp)
    {
      fprintf(stderr, "print_Lcurve: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  /* construct A and b, and calculate chi^2 = ||b - A c||^2 */
  poltor_build_ls(0, w);
  rnorm = sqrt(w->chisq);

  /* compute v = L c; L is stored in w->L by poltor_solve() */
  for (i = 0; i < p; ++i)
    {
      gsl_complex ci = gsl_vector_complex_get(w->c, i);
      double li = gsl_vector_get(w->L, i);
      gsl_complex val = gsl_complex_mul_real(ci, li);

      gsl_vector_complex_set(&v.vector, i, val);
    }

  /* compute || L c || */
  Lnorm = gsl_blas_dznrm2(&v.vector);

  fprintf(fp, "%.12e %.12e %.6e %.6e %.6e\n",
          log(rnorm),
          log(Lnorm),
          w->alpha_int,
          w->alpha_sh,
          w->alpha_tor);

  printcv_octave(w->residuals, "r");
  printcv_octave(w->c, "c");
  printv_octave(w->L, "L");

  fclose(fp);
#endif/*XXX*/

  return s;
} /* print_Lcurve() */

int
poltor_print_residual_stat(const char *component_str, const gsl_rstat_workspace *rstat_p)
{
  if (component_str == NULL)
    {
      /* print header */
      fprintf(stderr, "%12s %10s %12s %12s %12s\n",
              "", "N", "mean (nT)", "sigma (nT)", "rms (nT)");

    }
  else
    {
      const size_t n = gsl_rstat_n(rstat_p);

      if (n > 0)
        {
          fprintf(stderr, "%12s %10zu %12.2f %12.2f %12.2f\n",
                  component_str,
                  n,
                  gsl_rstat_mean(rstat_p),
                  gsl_rstat_sd(rstat_p),
                  gsl_rstat_rms(rstat_p));
        }
    }

  return GSL_SUCCESS;
}

int
poltor_print_residuals(const char *prefix, const size_t iter, poltor_workspace *w)
{
  int s = 0;
  char buf[2048];
  FILE *fp[12];
  const size_t n = 12; /* number of components to write to disk */
  const char *fmtstr = "%ld %9.4f %8.4f %8.4f %9.4f %.3e %.3e %.3e %.4f %.4f %.4f\n";
  const char *fmtstr_F = "%ld %9.4f %8.4f %8.4f %9.4f %.3e %.3e %.3e %.4f %.4f %.4f\n";
  const char *fmtstr_grad = "%ld %9.4f %8.4f %8.4f %9.4f %.3e %.3e %.3e %.4f %.4f %.4f\n";
  const double qdlat_cutoff = 55.0; /* cutoff latitude for high/low statistics */
  size_t i;
  magdata_list *list = w->data;
  size_t idx = 0;
  gsl_rstat_workspace *rstat_x = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_y = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_z = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_f = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_lowz = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_highz = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_lowf = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_highf = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_dx_ns = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_dy_ns = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_low_dz_ns = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_high_dz_ns = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_df_ns = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_dx_ew = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_dy_ew = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_low_dz_ew = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_high_dz_ew = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_df_ew = gsl_rstat_alloc();

  fprintf(stderr, "\n");

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);
      size_t j, k;

      gsl_rstat_reset(rstat_x);
      gsl_rstat_reset(rstat_y);
      gsl_rstat_reset(rstat_z);
      gsl_rstat_reset(rstat_f);
      gsl_rstat_reset(rstat_lowz);
      gsl_rstat_reset(rstat_highz);
      gsl_rstat_reset(rstat_lowf);
      gsl_rstat_reset(rstat_highf);
      gsl_rstat_reset(rstat_dx_ns);
      gsl_rstat_reset(rstat_dy_ns);
      gsl_rstat_reset(rstat_low_dz_ns);
      gsl_rstat_reset(rstat_high_dz_ns);
      gsl_rstat_reset(rstat_df_ns);
      gsl_rstat_reset(rstat_dx_ew);
      gsl_rstat_reset(rstat_dy_ew);
      gsl_rstat_reset(rstat_low_dz_ew);
      gsl_rstat_reset(rstat_high_dz_ew);
      gsl_rstat_reset(rstat_df_ew);

      sprintf(buf, "%s/res%zu_X_iter%zu.dat", prefix, i, iter);
      fp[0] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_Y_iter%zu.dat", prefix, i, iter);
      fp[1] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_Z_iter%zu.dat", prefix, i, iter);
      fp[2] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_F_iter%zu.dat", prefix, i, iter);
      fp[3] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_DX_NS_iter%zu.dat", prefix, i, iter);
      fp[4] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_DY_NS_iter%zu.dat", prefix, i, iter);
      fp[5] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_DZ_NS_iter%zu.dat", prefix, i, iter);
      fp[6] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_DF_NS_iter%zu.dat", prefix, i, iter);
      fp[7] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_DX_EW_iter%zu.dat", prefix, i, iter);
      fp[8] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_DY_EW_iter%zu.dat", prefix, i, iter);
      fp[9] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_DZ_EW_iter%zu.dat", prefix, i, iter);
      fp[10] = fopen(buf, "w");

      sprintf(buf, "%s/res%zu_DF_EW_iter%zu.dat", prefix, i, iter);
      fp[11] = fopen(buf, "w");

      /* header line */
      fprintf(fp[0], "# X vector residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[1], "# Y vector residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[2], "# Z vector residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[3], "# F scalar residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[4], "# DX gradient (N/S) vector residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[5], "# DY gradient (N/S) vector residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[6], "# DZ gradient (N/S) vector residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[7], "# DZ gradient (N/S) scalar residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[8], "# DX gradient (E/W) vector residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[9], "# DY gradient (E/W) vector residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[10], "# DZ gradient (E/W) vector residuals (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[11], "# DF gradient (E/W) scalar residuals (satellite %zu, iteration %zu)\n", i, iter);

      for (j = 0; j < n; ++j)
        {
          k = 1;
          fprintf(fp[j], "# Field %zu: timestamp (UT seconds since 1970-01-01)\n", k++);
          fprintf(fp[j], "# Field %zu: longitude (degrees)\n", k++);
          fprintf(fp[j], "# Field %zu: geocentric latitude (degrees)\n", k++);
          fprintf(fp[j], "# Field %zu: QD latitude (degrees)\n", k++);
          fprintf(fp[j], "# Field %zu: geocentric altitude (km)\n", k++);
          fprintf(fp[j], "# Field %zu: spatial weight factor\n", k++);
          fprintf(fp[j], "# Field %zu: robust weight factor\n", k++);
          fprintf(fp[j], "# Field %zu: total weight factor\n", k++);
        }

      fprintf(fp[0], "# Field %zu: X vector measurement (observation minus a priori model) (nT)\n", k);
      fprintf(fp[1], "# Field %zu: Y vector measurement (observation minus a priori model) (nT)\n", k);
      fprintf(fp[2], "# Field %zu: Z vector measurement (observation minus a priori model) (nT)\n", k);
      fprintf(fp[3], "# Field %zu: F scalar measurement (no subtraction of prior model) (nT)\n", k);
      fprintf(fp[4], "# Field %zu: DX N/S vector difference (a priori models subtracted) (nT)\n", k);
      fprintf(fp[5], "# Field %zu: DY N/S vector difference (a priori models subtracted) (nT)\n", k);
      fprintf(fp[6], "# Field %zu: DZ N/S vector difference (a priori models subtracted) (nT)\n", k);
      fprintf(fp[7], "# Field %zu: DF N/S scalar difference (no subtraction of prior models) (nT)\n", k);
      fprintf(fp[8], "# Field %zu: X vector measurement (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y vector measurement (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z vector measurement (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F scalar measurement (nT)\n", k);
      ++k;

      fprintf(fp[0], "# Field %zu: X fitted model (nT)\n", k);
      fprintf(fp[1], "# Field %zu: Y fitted model (nT)\n", k);
      fprintf(fp[2], "# Field %zu: Z fitted model (nT)\n", k);
      fprintf(fp[3], "# Field %zu: F total model || B_prior + B_fitted || (nT)\n", k);
      fprintf(fp[4], "# Field %zu: DX N/S fitted model (nT)\n", k);
      fprintf(fp[5], "# Field %zu: DY N/S fitted model (nT)\n", k);
      fprintf(fp[6], "# Field %zu: DZ N/S fitted model (nT)\n", k);
      fprintf(fp[7], "# Field %zu: DF N/S difference of fitted models || B_fitted_NS + B_prior_NS || - || B_fitted + B_prior || (nT)\n", k);

      fprintf(fp[8], "# Field %zu: X a priori model (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y a priori model (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z a priori model (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F a priori model (nT)\n", k);
      ++k;

      fprintf(fp[0], "# Field %zu: X residual (nT)\n", k);
      fprintf(fp[1], "# Field %zu: Y residual (nT)\n", k);
      fprintf(fp[2], "# Field %zu: Z residual (nT)\n", k);
      fprintf(fp[3], "# Field %zu: scalar residual (nT)\n", k);
      fprintf(fp[4], "# Field %zu: DX N/S residual (nT)\n", k);
      fprintf(fp[5], "# Field %zu: DY N/S residual (nT)\n", k);
      fprintf(fp[6], "# Field %zu: DZ N/S residual (nT)\n", k);
      fprintf(fp[7], "# Field %zu: DF N/S residual (nT)\n", k);

      fprintf(fp[8], "# Field %zu: X fitted model (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y fitted model (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z fitted model (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F fitted model (nT)\n", k);
      ++k;

      fprintf(fp[8], "# Field %zu: X vector measurement at E/W gradient point (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y vector measurement at E/W gradient point (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z vector measurement at E/W gradient point (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F scalar measurement at E/W gradient point (nT)\n", k);
      ++k;

      fprintf(fp[8], "# Field %zu: X a priori model at E/W gradient point (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y a priori model at E/W gradient point (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z a priori model at E/W gradient point (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F a priori model at E/W gradient point (nT)\n", k);
      ++k;

      fprintf(fp[8], "# Field %zu: X fitted model at E/W gradient point (nT)\n", k);
      fprintf(fp[9], "# Field %zu: Y fitted model at E/W gradient point (nT)\n", k);
      fprintf(fp[10], "# Field %zu: Z fitted model at E/W gradient point (nT)\n", k);
      fprintf(fp[11], "# Field %zu: F fitted model at E/W gradient point (nT)\n", k);
      ++k;

      fprintf(fp[8], "# Field %zu: DX E/W residual (nT)\n", k);
      fprintf(fp[9], "# Field %zu: DY E/W residual (nT)\n", k);
      fprintf(fp[10], "# Field %zu: DZ E/W residual (nT)\n", k);
      ++k;

      for (j = 0; j < mptr->n; ++j)
        {
          time_t unix_time = satdata_epoch2timet(mptr->t[j]);
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          double lon = phi * 180.0 / M_PI;
          double lat = 90.0 - theta * 180.0 / M_PI;
          double alt = r - mptr->R;
          double qdlat = mptr->qdlat[j];
          double B_obs[4], B_obs_grad[4];     /* observations in NEC */
          double B_prior[4], B_prior_grad[4]; /* prior models (main, crust, external) */
          double B_fit[4], B_grad_fit[4];     /* fitted (ionosphere) field model */
          double B_model[4], B_model_grad[4]; /* B_fit + B_prior */
          double res[4], res_grad[4];         /* residuals */

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute B_prior = B_main + B_crust + B_ext */
          magdata_prior(j, B_prior, mptr);

          B_obs[0] = mptr->Bx_nec[j];
          B_obs[1] = mptr->By_nec[j];
          B_obs[2] = mptr->Bz_nec[j];
          B_obs[3] = mptr->F[j];

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              /* evaluate model */
              poltor_eval_B(mptr->t[j], r, theta, phi, B_fit, w);

              if (MAGDATA_ExistDX_NS(mptr->flags[j]) || MAGDATA_ExistDY_NS(mptr->flags[j]) ||
                  MAGDATA_ExistDZ_NS(mptr->flags[j]) || MAGDATA_ExistDF_NS(mptr->flags[j]) ||
                  MAGDATA_ExistDX_EW(mptr->flags[j]) || MAGDATA_ExistDY_EW(mptr->flags[j]) ||
                  MAGDATA_ExistDZ_EW(mptr->flags[j]) || MAGDATA_ExistDF_EW(mptr->flags[j]))
                {
                  /* compute B_prior_grad = B_main + B_crust + B_ext */
                  magdata_prior_grad(j, B_prior_grad, mptr);

                  B_obs_grad[0] = mptr->Bx_nec_ns[j];
                  B_obs_grad[1] = mptr->By_nec_ns[j];
                  B_obs_grad[2] = mptr->Bz_nec_ns[j];
                  B_obs_grad[3] = mptr->F_ns[j];

                  poltor_eval_B(mptr->t_ns[j], mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j], B_grad_fit, w);
                }

              /* calculate vector residuals */
              for (k = 0; k < 3; ++k)
                {
                  B_model[k] = B_fit[k] + B_prior[k];
                  B_model_grad[k] = B_grad_fit[k] + B_prior_grad[k];

                  res[k] = B_obs[k] - B_model[k];
                  res_grad[k] = (B_obs_grad[k] - B_model_grad[k]) - res[k];
                }

              B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);
              B_model_grad[3] = gsl_hypot3(B_model_grad[0], B_model_grad[1], B_model_grad[2]);

              res[3] = B_obs[3] - B_model[3];
              res_grad[3] = (B_obs_grad[3] - B_model_grad[3]) - res[3];

              for (k = 0; k < 4; ++k)
                res_grad[k] *= mptr->satdir[j];
            }

          if ((j > 0) && (mptr->flags[j] & MAGDATA_FLG_TRACK_START))
            {
              for (k = 0; k < n; ++k)
                fprintf(fp[k], "\n\n");
            }

          if (MAGDATA_ExistX(mptr->flags[j]))
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double ws = gsl_vector_get(w->wts_spatial, idx);
                  double wr = gsl_vector_get(w->wts_robust, idx);
                  double wf = gsl_vector_get(w->wts_final, idx);

                  fprintf(fp[0], fmtstr, unix_time, lon, lat, qdlat, alt, ws, wr, wf, B_obs[0] - B_prior[0], B_fit[0], res[0]);
                  gsl_rstat_add(res[0], rstat_x);
                }

              ++idx;
            }

          if (MAGDATA_ExistY(mptr->flags[j]))
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double ws = gsl_vector_get(w->wts_spatial, idx);
                  double wr = gsl_vector_get(w->wts_robust, idx);
                  double wf = gsl_vector_get(w->wts_final, idx);
                  fprintf(fp[1], fmtstr, unix_time, lon, lat, qdlat, alt, ws, wr, wf, B_obs[1] - B_prior[1], B_fit[1], res[1]);
                  gsl_rstat_add(res[1], rstat_y);
                }

              ++idx;
            }

          if (MAGDATA_ExistZ(mptr->flags[j]))
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double ws = gsl_vector_get(w->wts_spatial, idx);
                  double wr = gsl_vector_get(w->wts_robust, idx);
                  double wf = gsl_vector_get(w->wts_final, idx);
                  fprintf(fp[2], fmtstr, unix_time, lon, lat, qdlat, alt, ws, wr, wf, B_obs[2] - B_prior[2], B_fit[2], res[2]);

                  gsl_rstat_add(res[2], rstat_z);

                  if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                    gsl_rstat_add(res[2], rstat_lowz);
                  else
                    gsl_rstat_add(res[2], rstat_highz);
                }

              ++idx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            {
              double ws = gsl_vector_get(w->wts_spatial, idx);
              double wr = gsl_vector_get(w->wts_robust, idx);
              double wf = gsl_vector_get(w->wts_final, idx);

              fprintf(fp[3], fmtstr_F, unix_time, lon, lat, qdlat, alt, ws, wr, wf, B_obs[3], B_model[3], res[3]);
              gsl_rstat_add(res[3], rstat_f);

              if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                gsl_rstat_add(res[3], rstat_lowf);
              else
                gsl_rstat_add(res[3], rstat_highf);

              ++idx;
            }

          if (MAGDATA_ExistDX_NS(mptr->flags[j]))
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double ws = gsl_vector_get(w->wts_spatial, idx);
                  double wr = gsl_vector_get(w->wts_robust, idx);
                  double wf = gsl_vector_get(w->wts_final, idx);
                  fprintf(fp[4], fmtstr_grad, unix_time, lon, lat, qdlat, alt, ws, wr, wf, (B_obs_grad[0] - B_prior_grad[0]) - (B_obs[0] - B_prior[0]), B_grad_fit[0] - B_fit[0], res_grad[0]);
                  gsl_rstat_add(res_grad[0], rstat_dx_ns);
                }

              ++idx;
            }

          if (MAGDATA_ExistDY_NS(mptr->flags[j]))
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double ws = gsl_vector_get(w->wts_spatial, idx);
                  double wr = gsl_vector_get(w->wts_robust, idx);
                  double wf = gsl_vector_get(w->wts_final, idx);
                  fprintf(fp[5], fmtstr_grad, unix_time, lon, lat, qdlat, alt, ws, wr, wf, (B_obs_grad[1] - B_prior_grad[1]) - (B_obs[1] - B_prior[1]), B_grad_fit[1] - B_fit[1], res_grad[1]);
                  gsl_rstat_add(res_grad[1], rstat_dy_ns);
                }

              ++idx;
            }

          if (MAGDATA_ExistDZ_NS(mptr->flags[j]))
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double ws = gsl_vector_get(w->wts_spatial, idx);
                  double wr = gsl_vector_get(w->wts_robust, idx);
                  double wf = gsl_vector_get(w->wts_final, idx);

                  fprintf(fp[6], fmtstr_grad, unix_time, lon, lat, qdlat, alt, ws, wr, wf, (B_obs_grad[2] - B_prior_grad[2]) - (B_obs[2] - B_prior[2]), B_grad_fit[2] - B_fit[2], res_grad[2]);

                  if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                    gsl_rstat_add(res_grad[2], rstat_low_dz_ns);
                  else
                    gsl_rstat_add(res_grad[2], rstat_high_dz_ns);
                }

              ++idx;
            }

          if (MAGDATA_ExistDF_NS(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            {
              double ws = gsl_vector_get(w->wts_spatial, idx);
              double wr = gsl_vector_get(w->wts_robust, idx);
              double wf = gsl_vector_get(w->wts_final, idx);
              fprintf(fp[7], fmtstr_grad, unix_time, lon, lat, qdlat, alt, ws, wr, wf, B_obs_grad[3] - B_obs[3], B_model_grad[3] - B_model[3], res_grad[3]);
              gsl_rstat_add(res_grad[3], rstat_df_ns);

              ++idx;
            }

          if (MAGDATA_ExistDX_EW(mptr->flags[j]))
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double ws = gsl_vector_get(w->wts_spatial, idx);
                  double wr = gsl_vector_get(w->wts_robust, idx);
                  double wf = gsl_vector_get(w->wts_final, idx);
                  fprintf(fp[8], fmtstr_grad, unix_time, lon, lat, qdlat, alt, ws, wr, wf, (B_obs_grad[0] - B_prior_grad[0]) - (B_obs[0] - B_prior[0]), B_grad_fit[0] - B_fit[0], res_grad[0]);
                  gsl_rstat_add(res_grad[0], rstat_dx_ew);
                }

              ++idx;
            }

          if (MAGDATA_ExistDY_EW(mptr->flags[j]))
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double ws = gsl_vector_get(w->wts_spatial, idx);
                  double wr = gsl_vector_get(w->wts_robust, idx);
                  double wf = gsl_vector_get(w->wts_final, idx);
                  fprintf(fp[9], fmtstr_grad, unix_time, lon, lat, qdlat, alt, ws, wr, wf, (B_obs_grad[1] - B_prior_grad[1]) - (B_obs[1] - B_prior[1]), B_grad_fit[1] - B_fit[1], res_grad[1]);
                  gsl_rstat_add(res_grad[1], rstat_dx_ew);
                }

              ++idx;
            }

          if (MAGDATA_ExistDZ_EW(mptr->flags[j]))
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double ws = gsl_vector_get(w->wts_spatial, idx);
                  double wr = gsl_vector_get(w->wts_robust, idx);
                  double wf = gsl_vector_get(w->wts_final, idx);
                  fprintf(fp[10], fmtstr_grad, unix_time, lon, lat, qdlat, alt, ws, wr, wf, (B_obs_grad[2] - B_prior_grad[2]) - (B_obs[2] - B_prior[2]), B_grad_fit[2] - B_fit[2], res_grad[2]);

                  if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                    gsl_rstat_add(res_grad[2], rstat_low_dz_ew);
                  else
                    gsl_rstat_add(res_grad[2], rstat_high_dz_ew);
                }

              ++idx;
            }

          if (MAGDATA_ExistDF_EW(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            {
              double ws = gsl_vector_get(w->wts_spatial, idx);
              double wr = gsl_vector_get(w->wts_robust, idx);
              double wf = gsl_vector_get(w->wts_final, idx);
              /*fprintf(fp[11], fmtstr_grad, unix_time, lon, lat, qdlat, alt, ws, wr, wf, B[3], B_model[3], B_fit[3], B_grad[3], B_grad_model[3], B_grad_fit[3]);
              gsl_rstat_add(B[3] - B_model[3] - B_fit[3] - (B_grad[3] - B_grad_model[3] - B_grad_fit[3]), rstat_df_ew);*/

              ++idx;
            }
        }

      fprintf(stderr, "=== FIT STATISTICS SATELLITE %zu ===\n", i);

      /* print header */
      poltor_print_residual_stat(NULL, NULL);

      poltor_print_residual_stat("X", rstat_x);
      poltor_print_residual_stat("Y", rstat_y);
      poltor_print_residual_stat("Z", rstat_z);
      poltor_print_residual_stat("F", rstat_f);
      poltor_print_residual_stat("low Z", rstat_lowz);
      poltor_print_residual_stat("high Z", rstat_highz);
      poltor_print_residual_stat("low F", rstat_lowf);
      poltor_print_residual_stat("high F", rstat_highf);

      poltor_print_residual_stat("N/S DX", rstat_dx_ns);
      poltor_print_residual_stat("N/S DY", rstat_dy_ns);
      poltor_print_residual_stat("low N/S DZ", rstat_low_dz_ns);
      poltor_print_residual_stat("high N/S DZ", rstat_high_dz_ns);
      poltor_print_residual_stat("N/S DF", rstat_df_ns);

      poltor_print_residual_stat("E/W DX", rstat_dx_ew);
      poltor_print_residual_stat("E/W DY", rstat_dy_ew);
      poltor_print_residual_stat("low E/W DZ", rstat_low_dz_ew);
      poltor_print_residual_stat("high E/W DZ", rstat_high_dz_ew);
      poltor_print_residual_stat("E/W DF", rstat_df_ew);
    }

  assert(idx == w->n);

  for (i = 0; i < n; ++i)
    fclose(fp[i]);

  gsl_rstat_free(rstat_x);
  gsl_rstat_free(rstat_y);
  gsl_rstat_free(rstat_z);
  gsl_rstat_free(rstat_f);
  gsl_rstat_free(rstat_lowf);
  gsl_rstat_free(rstat_highf);
  gsl_rstat_free(rstat_lowz);
  gsl_rstat_free(rstat_highz);
  gsl_rstat_free(rstat_dx_ns);
  gsl_rstat_free(rstat_dy_ns);
  gsl_rstat_free(rstat_low_dz_ns);
  gsl_rstat_free(rstat_high_dz_ns);
  gsl_rstat_free(rstat_df_ns);
  gsl_rstat_free(rstat_dx_ew);
  gsl_rstat_free(rstat_dy_ew);
  gsl_rstat_free(rstat_low_dz_ew);
  gsl_rstat_free(rstat_high_dz_ew);
  gsl_rstat_free(rstat_df_ew);

  return s;
}

#if 0
int
print_residuals(const char *filename, poltor_workspace *w)
{
#if 0/*XXX*/
  size_t i, j;
  FILE *fp;
  magdata *data = w->data;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_residuals: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
  fprintf(fp, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: X observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: DX observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: DY observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: DZ observation (nT)\n", i++);
  fprintf(fp, "# Field %zu: X residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: DX residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: DY residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: DZ residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal Z model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal Z model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal Z model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal DX model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal DY model (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal DZ model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal DX model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal DY model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal DZ model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal DX model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal DY model (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell toroidal DZ model (nT)\n", i++);
  fprintf(fp, "# Field %zu: vector along-track gradient available (1 or 0)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      double B_obs[4], B_model[3], B_int[4], B_sh[4], B_ext[4], B_tor[4], B_res[3];
      double dB_model[3], B_int2[4], B_sh2[4], B_ext2[4], B_tor2[4];
      double dB_obs[4] = { 0.0, 0.0, 0.0, 0.0 };
      double dB_res[3] = { 0.0, 0.0, 0.0 };
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double lt = get_localtime(unix_time, data->phi[i]);
      double theta = poltor_theta(i, w);

      /* store observation vector */
      magdata_residual(i, B_obs, data);

      /* compute individual magnetic field models at this point */
      poltor_eval_B_all(data->r[i], theta, data->phi[i], B_int, B_ext, B_sh, B_tor, w);

      if (data->flags[i] & MAGDATA_FLG_DZ_NS)
        {
          double theta_ns = poltor_theta_ns(i, w);

          magdata_residual_dB_ns(i, dB_obs, data);

          /* evaluate model at along-track point */
          poltor_eval_B_all(data->r_ns[i], theta_ns, data->phi_ns[i], B_int2, B_ext2, B_sh2, B_tor2, w);
        }

      /* compute total model and residual vector */
      for (j = 0; j < 3; ++j)
        {
          B_model[j] = B_int[j] + B_ext[j] + B_sh[j] + B_tor[j];
          B_res[j] = B_obs[j] - B_model[j];

          if (data->flags[i] & MAGDATA_FLG_DZ_NS)
            {
              dB_model[j] = B_int2[j] + B_ext2[j] + B_sh2[j] + B_tor2[j] - B_model[j];
              dB_res[j] = dB_obs[j] - dB_model[j];
            }
        }

      fprintf(fp, "%f %.2f %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d\n",
              satdata_epoch2year(data->t[i]),
              lt,
              data->r[i] - w->R,
              data->phi[i] * 180.0 / M_PI,
              90.0 - data->theta[i] * 180.0 / M_PI,
              data->qdlat[i],
              B_obs[0],
              B_obs[1],
              B_obs[2],
              data->satdir[i] * dB_obs[0],
              data->satdir[i] * dB_obs[1],
              data->satdir[i] * dB_obs[2],
              B_res[0],
              B_res[1],
              B_res[2],
              data->satdir[i] * dB_res[0],
              data->satdir[i] * dB_res[1],
              data->satdir[i] * dB_res[2],
              B_int[0],
              B_int[1],
              B_int[2],
              B_sh[0],
              B_sh[1],
              B_sh[2],
              B_tor[0],
              B_tor[1],
              B_tor[2],
              data->satdir[i] * (B_int2[0] - B_int[0]),
              data->satdir[i] * (B_int2[1] - B_int[1]),
              data->satdir[i] * (B_int2[2] - B_int[2]),
              data->satdir[i] * (B_sh2[0] - B_sh[0]),
              data->satdir[i] * (B_sh2[1] - B_sh[1]),
              data->satdir[i] * (B_sh2[2] - B_sh[2]),
              data->satdir[i] * (B_tor2[0] - B_tor[0]),
              data->satdir[i] * (B_tor2[1] - B_tor[1]),
              data->satdir[i] * (B_tor2[2] - B_tor[2]),
              data->flags[i] & MAGDATA_FLG_DZ_NS ? 1 : 0);
    }

  fclose(fp);
#endif

  return 0;
} /* print_residuals() */

#endif/*XXX*/

int
print_chisq(const char *filename, poltor_workspace *w)
{
  FILE *fp;
  poltor_parameters *params = &(w->params);

  fp = fopen(filename, "a");
  if (!fp)
    {
      fprintf(stderr, "print_chisq: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fprintf(fp, "%zu %zu %zu %zu %zu %zu %e\n",
          params->nmax_int,
          params->mmax_int,
          params->nmax_tor,
          params->mmax_tor,
          params->nmax_sh,
          params->mmax_sh,
          w->chisq / w->dof);

  fclose(fp);

  return 0;
} /* print_chisq() */

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options] input1.dat input2.dat ...\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --alpha_int | -c alpha_int          - damping parameter for internal coefficients\n");
  fprintf(stderr, "\t --alpha_sh | -d alpha_sh            - damping parameter for shell poloidal coefficients\n");
  fprintf(stderr, "\t --alpha_tor | -j alpha_tor          - damping parameter for shell toroidal coefficients\n");
  fprintf(stderr, "\t --print_residuals | -r              - print residuals for each iteration\n");
  fprintf(stderr, "\t --coef_file | -o file               - coefficient output file\n");
  fprintf(stderr, "\t --chisq_file | -p file              - chi^2 output file (for L-curve analysis)\n");
  fprintf(stderr, "\t --lls_file | -l lls_file            - LS system file (matrix and rhs)\n");
  fprintf(stderr, "\t --lcurve_file | -k lcurve_file      - output file for L-curve data\n");
  fprintf(stderr, "\t --maxit | -n maxit                  - number of robust iterations\n");
  fprintf(stderr, "\t --print_data | -u                   - print data file\n");
}

int
main(int argc, char *argv[])
{
  int status;
  double alpha_int = -1.0;
  double alpha_sh = -1.0;
  double alpha_tor = -1.0;
  const double R = R_EARTH_KM;
  const double d = R + 350.0;   /* radius of current shell for gravity/diamag */
  char *config_file = "PT.cfg";
  char *output_prefix = "output";
  char *spectrum_file = "poltor.s";
  char *corr_file = "corr.dat";
  char *output_file = NULL;
  char *chisq_file = NULL;
  char *lls_file = NULL;
  char *Lcurve_file = NULL;
  magdata_list *mlist = NULL;
  poltor_workspace *poltor_p;
  poltor_parameters params;
  struct timeval tv0, tv1;
  int print_data = 0;
  int print_residuals = 0; /* print residuals for each iteration */
  size_t maxit = 0;        /* maximum iterations */
  double rmin, rmax;       /* min/max radii of sources (km) */
  int nsource;             /* number of different satellites */

  poltor_init_params(&params);

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "print_residuals", no_argument, NULL, 'r' },
          { "output_file", required_argument, NULL, 'o' },
          { "chisq_file", required_argument, NULL, 'p' },
          { "lls_file", required_argument, NULL, 'l' },
          { "lcurve_file", required_argument, NULL, 'k' },
          { "alpha_int", required_argument, NULL, 'c' },
          { "alpha_sh", required_argument, NULL, 'd' },
          { "alpha_tor", required_argument, NULL, 'j' },
          { "maxit", required_argument, NULL, 'n' },
          { "print_data", no_argument, NULL, 'u' },
          { "config_file", no_argument, NULL, 'C' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "c:C:d:j:k:l:o:p:n:ru", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'C':
            config_file = optarg;
            break;

          case 'c':
            alpha_int = atof(optarg);
            break;

          case 'd':
            alpha_sh = atof(optarg);
            break;

          case 'j':
            alpha_tor = atof(optarg);
            break;

          case 'r':
            print_residuals = 1;
            break;

          case 'k':
            Lcurve_file = optarg;
            break;

          case 'o':
            output_file = optarg;
            break;

          case 'p':
            chisq_file = optarg;
            break;

          case 'l':
            lls_file = optarg;
            break;

          case 'n':
            maxit = (size_t) atoi(optarg);
            break;

          case 'u':
            print_data = 1;
            break;

          default:
            break;
        }
    }

  nsource = argc - optind;
  if (nsource == 0)
    {
      print_help(argv);
      exit(1);
    }

  mlist = magdata_list_alloc(nsource);
  if (!mlist)
    {
      print_help(argv);
      exit(1);
    }

  /* read in magdata files */
  while (optind < argc)
    {
      size_t ndata;

      fprintf(stderr, "main: reading %s...", argv[optind]);
      gettimeofday(&tv0, NULL);
      ndata = magdata_list_add(argv[optind], mlist);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%zu data total, %g seconds)\n",
              ndata, time_diff(tv0, tv1));

      ++optind;
    }

  /* parse configuration file */
  fprintf(stderr, "main: parsing configuration file %s...", config_file);
  status = parse_config_file(config_file, &params);
  fprintf(stderr, "done (status = %d)\n", status);
  if (status)
    exit(1);

  /* check if command line arguments override config file values */
  if (alpha_int >= 0.0)
    params.alpha_int = alpha_int;
  if (alpha_sh >= 0.0)
    params.alpha_sh = alpha_sh;
  if (alpha_tor >= 0.0)
    params.alpha_tor = alpha_tor;
  if (maxit > 0)
    params.max_iter = maxit;

  if (params.synth_data)
    {
      params.nmax_int = 10;
      params.mmax_int = 10;
      params.nmax_ext = 3;
      params.mmax_ext = 2;
      params.nmax_sh = 4;
      params.mmax_sh = 4;
      params.nmax_tor = 5;
      params.mmax_tor = 3;

      /* turn off weights for synthetic test */
      params.use_weights = 0;

      /* turn off regularization */
      params.regularize = 0;
    }

  params.mmax_int = GSL_MIN(params.mmax_int, params.nmax_int);
  params.mmax_ext = GSL_MIN(params.mmax_ext, params.nmax_ext);
  params.mmax_sh = GSL_MIN(params.mmax_sh, params.nmax_sh);
  params.mmax_tor = GSL_MIN(params.mmax_tor, params.nmax_tor);

  fprintf(stderr, "main: nmax_int  = %zu\n", params.nmax_int);
  fprintf(stderr, "main: mmax_int  = %zu\n", params.mmax_int);
  fprintf(stderr, "main: nmax_ext  = %zu\n", params.nmax_ext);
  fprintf(stderr, "main: mmax_ext  = %zu\n", params.mmax_ext);
  fprintf(stderr, "main: nmax_sh   = %zu\n", params.nmax_sh);
  fprintf(stderr, "main: mmax_sh   = %zu\n", params.mmax_sh);
  fprintf(stderr, "main: nmax_tor  = %zu\n", params.nmax_tor);
  fprintf(stderr, "main: mmax_tor  = %zu\n", params.mmax_tor);
  fprintf(stderr, "main: alpha_int = %g\n", params.alpha_int);
  fprintf(stderr, "main: alpha_sh  = %g\n", params.alpha_sh);
  fprintf(stderr, "main: alpha_tor = %g\n", params.alpha_tor);
  fprintf(stderr, "main: shell_J   = %zu\n", params.shell_J);
  fprintf(stderr, "main: R         = %g [km]\n", params.R);
  fprintf(stderr, "main: b         = %g [km]\n", params.b);

  if (Lcurve_file)
    fprintf(stderr, "main: L-curve file  = %s\n", Lcurve_file);

  /* XXX */
  {
    fprintf(stderr, "main: replacing theta with thetaq...");
    replace_qdlat(mlist);
    fprintf(stderr, "done\n");
  }

  /*
   * re-compute flags for fitting components / gradient, etc;
   * must be called before magdata_init()
   */
  set_flags(&params, mlist);

#if 0/*XXX*/
  fprintf(stderr, "main: initializing spatial weighting histogram...");
  gettimeofday(&tv0, NULL);
  magdata_init(mlist);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* re-compute weights, nvec, nres based on flags update */
  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mlist);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
#endif

#if 0/*XXX*/
  if (params.synth_data)
    {
      fprintf(stderr, "main: setting unit spatial weights...");
      magdata_unit_weights(mlist);
      fprintf(stderr, "done\n");
    }
#endif

  magdata_list_rminmax(mlist, &rmin, &rmax);

  fprintf(stderr, "main: satellite rmin = %.1f (%.1f) [km]\n", rmin, rmin - R);
  fprintf(stderr, "main: satellite rmax = %.1f (%.1f) [km]\n", rmax, rmax - R);

  params.d = d;
  params.rmin = GSL_MAX(rmin, R + 250.0);
  params.rmax = GSL_MIN(rmax, R + 450.0);
  params.data = mlist;

#if POLTOR_QD_HARMONICS
  params.flags = POLTOR_FLG_QD_HARMONICS;
#else
  params.flags = 0;
#endif

  poltor_p = poltor_alloc(&params);

  fprintf(stderr, "main: poltor rmin = %.1f (%.1f) [km]\n",
          params.rmin, params.rmin - R);
  fprintf(stderr, "main: poltor rmax = %.1f (%.1f) [km]\n",
          params.rmax, params.rmax - R);

  if (params.synth_data)
    {
      fprintf(stderr, "main: replacing with synthetic data...");
      gettimeofday(&tv0, NULL);
      poltor_synth(poltor_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
    }

  if (print_data)
    {
      fprintf(stderr, "main: writing data to %s...", output_prefix);
      magdata_list_print(output_prefix, mlist);
      fprintf(stderr, "done\n");

      fprintf(stderr, "main: writing data map to %s...", output_prefix);
      magdata_list_map(output_prefix, mlist);
      fprintf(stderr, "done\n");
    }

  if (lls_file)
    {
      /* use previously computed LS system from file */
      fprintf(stderr, "main: loading LS system from %s...", lls_file);
      lls_complex_load(lls_file, poltor_p->lls_workspace_p);
      fprintf(stderr, "done\n");

      /* solve LS system */
      poltor_solve(poltor_p);
    }
  else
    {
      size_t maxiter = params.max_iter;
      size_t iter = 0;
      char buf[2048];
      gsl_vector_complex *cprev = gsl_vector_complex_alloc(poltor_p->p);
      gsl_vector_complex *c = poltor_p->c; /* model coefficients */
      size_t *res_cnt = poltor_p->res_cnt;

      fprintf(stderr, "main: X residuals      = %zu\n", res_cnt[MAGDATA_LIST_IDX_X]);
      fprintf(stderr, "main: Y residuals      = %zu\n", res_cnt[MAGDATA_LIST_IDX_Y]);
      fprintf(stderr, "main: Z residuals      = %zu\n", res_cnt[MAGDATA_LIST_IDX_Z]);
      fprintf(stderr, "main: F residuals      = %zu\n", res_cnt[MAGDATA_LIST_IDX_F]);
      fprintf(stderr, "main: DX N/S residuals = %zu\n", res_cnt[MAGDATA_LIST_IDX_DX_NS]);
      fprintf(stderr, "main: DY N/S residuals = %zu\n", res_cnt[MAGDATA_LIST_IDX_DY_NS]);
      fprintf(stderr, "main: DZ N/S residuals = %zu\n", res_cnt[MAGDATA_LIST_IDX_DZ_NS]);
      fprintf(stderr, "main: DF N/S residuals = %zu\n", res_cnt[MAGDATA_LIST_IDX_DF_NS]);
      fprintf(stderr, "main: total residuals  = %zu\n", poltor_p->n);
      fprintf(stderr, "main: ncoeff           = %zu\n", poltor_p->p);
      fprintf(stderr, "main: nres / ncoeff    = %.1f\n", (double) poltor_p->n / (double) poltor_p->p);

      if (params.synth_data == 1 && poltor_p->lls_solution == 1)
        maxiter = 1;

      /* construct initial guess vector */
      initial_guess(c, poltor_p);

      while (iter++ < maxiter)
        {
          fprintf(stderr, "main: ROBUST ITERATION %zu/%zu\n", iter, maxiter);

          /* save current coefficient vector */
          gsl_vector_complex_memcpy(cprev, c);

          poltor_calc_nonlinear(iter, c, poltor_p);

          /* compute || c_{k+1} - c_k || */
          gsl_vector_complex_sub(cprev, c);
          fprintf(stderr, "main: || c_{k+1} - c_k || = %.12e\n", gsl_blas_dznrm2(cprev));

          /* output spectrum for this iteration */
          sprintf(buf, "%s/%s.iter%zu", output_prefix, spectrum_file, iter);
          fprintf(stderr, "main: printing spectrum to %s...", buf);
          poltor_print_spectrum(buf, c, poltor_p);
          fprintf(stderr, "done\n");

          /* output coefficient file for this iteration */
          if (output_file)
            {
              sprintf(buf, "%s/%s.iter%zu", output_prefix, output_file, iter);
              fprintf(stderr, "main: writing output coefficients to %s...", buf);
              poltor_write(buf, poltor_p);
              fprintf(stderr, "done\n");
            }

          /*
           * JHJ contains the Cholesky factor of J^H W J after calling poltor_calc_nonlinear(),
           * invert the matrix to find the covariance matrix
           */
          fprintf(stderr, "main: calculating model covariance matrix...");
          lapack_complex_cholesky_invert(poltor_p->JHJ);
          fprintf(stderr, "done\n");

          sprintf(buf, "%s/error.iter%zu", output_prefix, iter);
          fprintf(stderr, "main: writing coefficient uncertainties to %s...", buf);
          print_uncertanties(buf, poltor_p);
          fprintf(stderr, "done\n");

          sprintf(buf, "corr_iter%zu", iter);
          fprintf(stderr, "main: writing correlation matrix to %s...", buf);
          print_correlation(buf, poltor_p);
          fprintf(stderr, "done\n");

          if (print_residuals)
            {
              fprintf(stderr, "main: printing residuals to %s...", output_prefix);
              poltor_print_residuals(output_prefix, iter, poltor_p);
              fprintf(stderr, "done\n");
            }
        }

      gsl_vector_complex_free(cprev);
    }

  print_coefficients(poltor_p);

  fprintf(stderr, "main: printing spectrum to %s...", spectrum_file);
  poltor_print_spectrum(spectrum_file, poltor_p->c, poltor_p);
  fprintf(stderr, "done\n");

  if (Lcurve_file)
    {
      fprintf(stderr, "main: writing L-curve data to %s...", Lcurve_file);
      print_Lcurve(Lcurve_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  if (output_file)
    {
      fprintf(stderr, "main: writing output coefficients to %s...", output_file);
      poltor_write(output_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  if (chisq_file)
    {
      fprintf(stderr, "main: printing chisq/dof to %s...", chisq_file);
      print_chisq(chisq_file, poltor_p);
      fprintf(stderr, "done\n");
    }

  magdata_list_free(mlist);
  poltor_free(poltor_p);

  return 0;
} /* main() */
