/*
 * poltor_synth.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <mainlib/ml_magdata.h>
#include <mainlib/ml_green.h>

#include "poltor.h"

static int
poltor_synth_fill_random(const size_t nmax, const size_t mmax, gsl_vector_complex *v, gsl_rng *r)
{
  int s = 0;
  size_t n;

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;

      for (m = 0; m <= M; ++m)
        {
          size_t cidx = green_idx(n, m, mmax);
          double a = 2.0 * gsl_rng_uniform(r) - 1.0; /* [-1,1] */
          double b;
          
          if (m == 0)
            b = 0.0;
          else
            b = 2.0 * gsl_rng_uniform(r) - 1.0;

          gsl_vector_complex_set(v, cidx, gsl_complex_rect(a, b));

          /* g_{n,-m} = conj(g_{nm}) */
          if (m > 0)
            {
              cidx = green_idx(n, -m, mmax);
              gsl_vector_complex_set(v, cidx, gsl_complex_rect(a, -b));
            }
        }
    }

  return s;
}

/* initialize coefficient vector with random coefficients in [-1,1] */
int
poltor_synth_init(gsl_vector_complex *g, poltor_workspace *w)
{
  int s = 0;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_vector_complex_view v;

  gsl_vector_complex_set_zero(g);

  if (w->p_pint > 0)
    {
      v = gsl_vector_complex_subvector(g, w->pint_offset, w->p_pint);
      poltor_synth_fill_random(w->nmax_int, w->mmax_int, &v.vector, r);
    }

  if (w->p_pext > 0)
    {
      v = gsl_vector_complex_subvector(g, w->pext_offset, w->p_pext);
      poltor_synth_fill_random(w->nmax_ext, w->mmax_ext, &v.vector, r);
    }

  if (w->p_psh > 0)
    {
      size_t j;

      for (j = 0; j <= w->shell_J; ++j)
        {
          v = gsl_vector_complex_subvector(g, w->psh_offset + j * w->nnm_sh, w->nnm_sh);
          poltor_synth_fill_random(w->nmax_sh, w->mmax_sh, &v.vector, r);
        }
    }

  if (w->p_tor > 0)
    {
      v = gsl_vector_complex_subvector(g, w->tor_offset, w->p_tor);
      poltor_synth_fill_random(w->nmax_tor, w->mmax_tor, &v.vector, r);
    }

  gsl_rng_free(r);

  return s;
} /* poltor_synth_init() */

/* calculate synthetic field value for single point */
static int
poltor_synth_calc(const double r, const double theta, const double phi,
                  gsl_vector_complex *g, double B[3], poltor_workspace *w)
{
  int s = 0;
  const double invsint = 1.0 / sin(theta);
  complex double X = 0.0, Y = 0.0, Z = 0.0;
  double ratio;
  green_complex_workspace *green_p = w->green_p[0];
  complex double *Ynm = green_p->Ynm;
  complex double *dYnm = green_p->dYnm;
  size_t n, j;
  int m;

  /* compute Ynm and d/dtheta Ynm */
  green_complex_Ynm_deriv(theta, phi, green_p);

  /* internal potential field */
  ratio = w->R / r;
  for (n = 1; n <= w->nmax_int; ++n)
    {
      int ni = (int) GSL_MIN(n, w->mmax_int);
      double rterm = pow(ratio, n + 2.0);

      for (m = -ni; m <= ni; ++m)
        {
          size_t mabs = (size_t) abs(m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_complex val = gsl_vector_complex_get(g, cidx);
          complex double gnm = GSL_REAL(val) + I * GSL_IMAG(val);
          complex double Y_nm, dY_nm;

          if (m >= 0)
            {
              Y_nm = Ynm[pidx];
              dY_nm = dYnm[pidx];
            }
          else
            {
              Y_nm = conj(Ynm[pidx]);
              dY_nm = conj(dYnm[pidx]);
            }

          X += rterm * gnm * dY_nm;
          Y -= rterm * I * m * invsint * gnm * Y_nm;
          Z -= (n + 1.0) * rterm * gnm * Y_nm;
        }
    }

  /* external potential field */
  ratio = r / w->R;
  for (n = 1; n <= w->nmax_ext; ++n)
    {
      int ni = (int) GSL_MIN(n, w->mmax_ext);
      double rterm = pow(ratio, n - 1.0);

      for (m = -ni; m <= ni; ++m)
        {
          size_t mabs = (size_t) abs(m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          gsl_complex val = gsl_vector_complex_get(g, cidx);
          complex double knm = GSL_REAL(val) + I * GSL_IMAG(val);
          complex double Y_nm, dY_nm;

          if (m >= 0)
            {
              Y_nm = Ynm[pidx];
              dY_nm = dYnm[pidx];
            }
          else
            {
              Y_nm = conj(Ynm[pidx]);
              dY_nm = conj(dYnm[pidx]);
            }

          X += rterm * knm * dY_nm;
          Y -= rterm * I * m * invsint * knm * Y_nm;
          Z += n * rterm * knm * Y_nm;
        }
    }

  /* shell poloidal field */
  ratio = w->R / r;
  for (j = 0; j <= w->shell_J; ++j)
    {
      for (n = 1; n <= w->nmax_sh; ++n)
        {
          int ni = (int) GSL_MIN(n, w->mmax_sh);
          double An, Bn, ABsum, dABsum;

          poltor_shell_An(n, j, r, &An, w);
          poltor_shell_Bn(n, j, r, &Bn, w);
          ABsum = An + Bn;
          dABsum = -(double)n * An + (n + 1.0) * Bn;

          for (m = -ni; m <= ni; ++m)
            {
              size_t mabs = (size_t) abs(m);
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              size_t cidx = poltor_jnmidx(j, n, m, w);
              gsl_complex val = gsl_vector_complex_get(g, cidx);
              complex double gnm = GSL_REAL(val) + I * GSL_IMAG(val);
              complex double Y_nm, dY_nm;

              if (m >= 0)
                {
                  Y_nm = Ynm[pidx];
                  dY_nm = dYnm[pidx];
                }
              else
                {
                  Y_nm = conj(Ynm[pidx]);
                  dY_nm = conj(dYnm[pidx]);
                }

              X += ratio * gnm * dABsum * dY_nm;
              Y -= ratio * I * m * invsint * gnm * dABsum * Y_nm;
              Z += ratio * n * (n + 1.0) * gnm * ABsum * Y_nm;
            }
        }
    }

  /* shell toroidal field */
  for (n = 1; n <= w->nmax_tor; ++n)
    {
      int ni = (int) GSL_MIN(n, w->mmax_tor);

      for (m = -ni; m <= ni; ++m)
        {
          size_t mabs = (size_t) abs(m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          size_t cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
          gsl_complex val = gsl_vector_complex_get(g, cidx);
          complex double phinm = GSL_REAL(val) + I * GSL_IMAG(val);
          complex double Y_nm, dY_nm;

          if (m >= 0)
            {
              Y_nm = Ynm[pidx];
              dY_nm = dYnm[pidx];
            }
          else
            {
              Y_nm = conj(Ynm[pidx]);
              dY_nm = conj(dYnm[pidx]);
            }

          X += I * m * invsint * phinm * Y_nm;
          Y += phinm * dY_nm;
        }
    }

  B[0] = creal(X);
  B[1] = creal(Y);
  B[2] = creal(Z);

  return s;
} /* poltor_synth_calc() */

/* replace with synthetic data for testing */
int
poltor_synth(poltor_workspace *w)
{
  int s = 0;
  size_t i, j;
  gsl_vector_complex *g = gsl_vector_complex_calloc(w->p);
  magdata_list *list = w->data;

  s = poltor_synth_init(g, w);

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      for (j = 0; j < mptr->n; ++j)
        {
          double B[3];

          s += poltor_synth_calc(mptr->r[j], mptr->theta[j], mptr->phi[j], g, B, w);

          mptr->Bx_nec[j] = B[0];
          mptr->By_nec[j] = B[1];
          mptr->Bz_nec[j] = B[2];
          mptr->F[j] = gsl_hypot3(B[0], B[1], B[2]);

          mptr->Bx_model[j] = 0.0;
          mptr->By_model[j] = 0.0;
          mptr->Bz_model[j] = 0.0;

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DF_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW | MAGDATA_FLG_DF_EW))
            {
              double B_ns[3];

              s += poltor_synth_calc(mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j], g, B_ns, w);

              mptr->Bx_nec_ns[j] = B_ns[0];
              mptr->By_nec_ns[j] = B_ns[1];
              mptr->Bz_nec_ns[j] = B_ns[2];
              mptr->F_ns[j] = gsl_hypot3(B_ns[0], B_ns[1], B_ns[2]);

              mptr->Bx_model_ns[j] = 0.0;
              mptr->By_model_ns[j] = 0.0;
              mptr->Bz_model_ns[j] = 0.0;
            }
        }
    }

  gsl_vector_complex_free(g);

  return 0;
} /* poltor_synth() */

int
poltor_synth_print(poltor_workspace *w)
{
  int s = 0;
  const size_t nmax = 2; /* maximum nmax for printing */
  const char *synth_spectrum_file = "poltor_synth.s";
  gsl_vector_complex *g_synth = gsl_vector_complex_alloc(w->p);
  size_t n;

  poltor_synth_init(g_synth, w);

  /* print internal poloidal coefficients */
  fprintf(stderr, "Internal poloidal coefficients:\n");
  for (n = 1; n <= GSL_MIN(nmax, w->nmax_int); ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_int);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_complex g_nm = poltor_get(cidx, w);
          gsl_complex synth_nm = gsl_vector_complex_get(g_synth, cidx);

          fprintf(stderr, "g(%2zu,%2d) = (%7.4f,%7.4f) [truth: (%7.4f,%7.4f)] relerr = (%12.6e,%12.6e)\n",
                  n,
                  m,
                  GSL_REAL(g_nm),
                  GSL_IMAG(g_nm),
                  GSL_REAL(synth_nm),
                  GSL_IMAG(synth_nm),
                  fabs(GSL_REAL(g_nm) - GSL_REAL(synth_nm)) / GSL_REAL(synth_nm),
                  fabs(GSL_IMAG(g_nm) - GSL_IMAG(synth_nm)) / GSL_IMAG(synth_nm));
        }
    }

  /* print external poloidal coefficients */
  fprintf(stderr, "External poloidal coefficients:\n");
  for (n = 1; n <= GSL_MIN(nmax, w->nmax_ext); ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_ext);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          gsl_complex g_nm = poltor_get(cidx, w);
          gsl_complex synth_nm = gsl_vector_complex_get(g_synth, cidx);

          fprintf(stderr, "k(%2zu,%2d) = (%7.4f,%7.4f) [truth: (%7.4f,%7.4f)] relerr = (%12.6e,%12.6e)\n",
                  n,
                  m,
                  GSL_REAL(g_nm),
                  GSL_IMAG(g_nm),
                  GSL_REAL(synth_nm),
                  GSL_IMAG(synth_nm),
                  fabs(GSL_REAL(g_nm) - GSL_REAL(synth_nm)) / GSL_REAL(synth_nm),
                  fabs(GSL_IMAG(g_nm) - GSL_IMAG(synth_nm)) / GSL_IMAG(synth_nm));
        }
    }

  /* print shell poloidal coefficients */
  fprintf(stderr, "Shell poloidal coefficients:\n");
  {
    size_t j;

    for (j = 0; j <= w->shell_J; ++j)
      {
        for (n = 1; n <= GSL_MIN(nmax, w->nmax_sh); ++n)
          {
            int M = (int) GSL_MIN(n, w->mmax_sh);
            int m;

            for (m = -M; m <= M; ++m)
              {
                size_t cidx = poltor_jnmidx(j, n, m, w);
                gsl_complex g_nm = poltor_get(cidx, w);
                gsl_complex synth_nm = gsl_vector_complex_get(g_synth, cidx);

                /* only output j=0 terms */
                if (j == 0)
                  {
                    fprintf(stderr, "q(%2zu,%2d) = (%7.4f,%7.4f) [truth: (%7.4f,%7.4f)] relerr = (%12.6e,%12.6e)\n",
                            n,
                            m,
                            GSL_REAL(g_nm),
                            GSL_IMAG(g_nm),
                            GSL_REAL(synth_nm),
                            GSL_IMAG(synth_nm),
                            fabs(GSL_REAL(g_nm) - GSL_REAL(synth_nm)) / GSL_REAL(synth_nm),
                            fabs(GSL_IMAG(g_nm) - GSL_IMAG(synth_nm)) / GSL_IMAG(synth_nm));
                  }
              }
          }
      }
  }

  /* print toroidal coefficients */
  fprintf(stderr, "Shell toroidal coefficients:\n");
  for (n = 1; n <= GSL_MIN(nmax, w->nmax_tor); ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_tor);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
          gsl_complex g_nm = poltor_get(cidx, w);
          gsl_complex synth_nm = gsl_vector_complex_get(g_synth, cidx);

          fprintf(stderr, "phi(%2zu,%2d) = (%7.4f,%7.4f) [truth: (%7.4f,%7.4f)] relerr = (%12.6e,%12.6e)\n",
                  n,
                  m,
                  GSL_REAL(g_nm),
                  GSL_IMAG(g_nm),
                  GSL_REAL(synth_nm),
                  GSL_IMAG(synth_nm),
                  fabs(GSL_REAL(g_nm) - GSL_REAL(synth_nm)) / GSL_REAL(synth_nm),
                  fabs(GSL_IMAG(g_nm) - GSL_IMAG(synth_nm)) / GSL_IMAG(synth_nm));
        }
    }

  fprintf(stderr, "poltor_synth_print: printing synthetic spectrum to %s...", synth_spectrum_file);
  poltor_print_spectrum(synth_spectrum_file, g_synth, w);
  fprintf(stderr, "done\n");

  fprintf(stderr, "|| g || = %.12e\n", gsl_blas_dznrm2(w->c));
  fprintf(stderr, "|| g_synth || = %.12e\n", gsl_blas_dznrm2(g_synth));

  gsl_vector_complex_sub(g_synth, w->c);
  fprintf(stderr, "|| g - g_synth || = %.12e\n", gsl_blas_dznrm2(g_synth));

  return s;
}
