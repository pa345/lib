/*
 * mfield_synth.c
 *
 * Routines for handling synthetic test case
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_euler.h>

#include "mfield.h"
#include "mfield_fluxcal.h"
#include "mfield_synth.h"

static int mfield_synth_calc(const double t, const double r, const double theta, const double phi,
                             const gsl_vector * g, double B[3], mfield_workspace *w);
static int mfield_synth_calc_dBdt(const double t, const double r, const double theta, const double phi, const gsl_vector * g,
                                  double dBdt[3], mfield_workspace *w);

/* fill in internal coefficient vector of synthetic gauss coefficients */
int
mfield_synth_g(gsl_vector * g, mfield_workspace * w)
{
  msynth_workspace *core_p = msynth_shc_read(MSYNTH_CHAOS_FILE);
  msynth_workspace *crust_p = msynth_mf7_read(MSYNTH_MF7_FILE);

  gsl_vector_set_zero(g);

  /* initialize internal field coefficients using CHAOS */
  mfield_fill_g(g, core_p, crust_p, w);

  msynth_free(core_p);
  msynth_free(crust_p);

  return 0;
}

/* replace with synthetic data for testing */
int
mfield_synth_replace(mfield_workspace *w)
{
  const mfield_parameters *params = &(w->params);
  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng *rng_p = gsl_rng_alloc(T);
  gsl_vector *g = gsl_vector_alloc(w->p_int);
  mfield_workspace **mfield_array;
  size_t i, j;

  /* Euler angles */
  const double alpha = -13.1 * M_PI / 180.0 * 1.0;
  const double beta = -5.2 * M_PI / 180.0 * 1.0;
  const double gamma = 3.4 * M_PI / 180.0 * 1.0;

  /* fluxgate calibration parameters */
  const double cal_data[] = { 1.01, 0.98, 1.03, 40.0, -11.2, 13.2, 0.01 * M_PI / 180.0, -0.02 * M_PI / 180.0, 0.03 * M_PI / 180.0 };
  gsl_vector_const_view cal_params = gsl_vector_const_view_array(cal_data, FLUXCAL_P);

  /* initialize synthetic gauss coefficients */
  mfield_synth_g(g, w);

  mfield_array = malloc(w->max_threads * sizeof(mfield_workspace *));
  for (i = 0; i < w->max_threads; ++i)
    {
      mfield_array[i] = mfield_copy(w);
    }

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      int fit_euler = params->fit_euler && (mptr->global_flags & MAGDATA_GLOBFLG_EULER);
      int fit_fluxcal = params->fit_fluxcal && (mptr->global_flags & MAGDATA_GLOBFLG_FLUXCAL);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          mfield_workspace *mfield_p = mfield_array[thread_id];
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          double B[3];
          size_t k;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

#if MFIELD_SYNTH_HIGH_LAT_ONLY
          /* only replace field values with synthetic values at high latitudes */
          /*if (fabs(mptr->qdlat[j]) <= params->qdlat_fit_cutoff)*/
          if (fabs(mptr->qdlat[j]) <= 20.0)
            {
              mptr->Bx_nec[j] -= mptr->Bx_model[j];
              mptr->By_nec[j] -= mptr->By_model[j];
              mptr->Bz_nec[j] -= mptr->Bz_model[j];
              mptr->Bx_model[j] = 0.0;
              mptr->By_model[j] = 0.0;
              mptr->Bz_model[j] = 0.0;

              if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                    MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
                {
                  mptr->Bx_nec_ns[j] -= mptr->Bx_model_ns[j];
                  mptr->By_nec_ns[j] -= mptr->By_model_ns[j];
                  mptr->Bz_nec_ns[j] -= mptr->Bz_model_ns[j];
                  mptr->Bx_model_ns[j] = 0.0;
                  mptr->By_model_ns[j] = 0.0;
                  mptr->Bz_model_ns[j] = 0.0;
                }

              continue;
            }
#endif

          /* synthesize magnetic field vector */
          mfield_synth_calc(mptr->t[j], r, theta, phi, g, B, mfield_p);

          if (params->synth_noise)
            {
              /* add some noise to measurements with 0.1 nT sigma in all vector components */
              for (k = 0; k < 3; ++k)
                B[k] += gsl_ran_gaussian(rng_p, 0.1);
            }

          mptr->Bx_nec[j] = B[0];
          mptr->By_nec[j] = B[1];
          mptr->Bz_nec[j] = B[2];
          mptr->F[j] = gsl_hypot3(B[0], B[1], B[2]);

          /* no crustal/external field for synthetic data */
          mptr->Bx_model[j] = 0.0;
          mptr->By_model[j] = 0.0;
          mptr->Bz_model[j] = 0.0;

          /* rotate NEC vector to VFM frame */
          if (fit_euler)
            {
              double *q = &(mptr->q[4*j]);
              double B_vfm[3];

              euler_nec2vfm(mptr->euler_flags, alpha, beta, gamma, q, B, B_vfm);

              if (fit_fluxcal)
                mfield_fluxcal_invapply_datum(&cal_params.vector, B_vfm, B_vfm);

              mptr->Bx_vfm[j] = B_vfm[0];
              mptr->By_vfm[j] = B_vfm[1];
              mptr->Bz_vfm[j] = B_vfm[2];
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DXDT | MAGDATA_FLG_DYDT | MAGDATA_FLG_DZDT))
            {
              double dBdt[3];

              /* synthesize dB/dt vector */
              mfield_synth_calc_dBdt(mptr->t[j], r, theta, phi, g, dBdt, mfield_p);
              mptr->dXdt_nec[j] = dBdt[0];
              mptr->dYdt_nec[j] = dBdt[1];
              mptr->dZdt_nec[j] = dBdt[2];
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DF_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW | MAGDATA_FLG_DF_EW))
            {
              mfield_synth_calc(mptr->t_ns[j], mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j], g, B, mfield_p);

              mptr->Bx_nec_ns[j] = B[0];
              mptr->By_nec_ns[j] = B[1];
              mptr->Bz_nec_ns[j] = B[2];
              mptr->F_ns[j] = gsl_hypot3(B[0], B[1], B[2]);

              mptr->Bx_model_ns[j] = 0.0;
              mptr->By_model_ns[j] = 0.0;
              mptr->Bz_model_ns[j] = 0.0;

              /* rotate NEC vector to VFM frame */
              if (fit_euler)
                {
                  double *q = &(mptr->q_ns[4*j]);
                  double B_vfm[3];

                  euler_nec2vfm(mptr->euler_flags, alpha, beta, gamma, q, B, B_vfm);

                  mptr->Bx_vfm_ns[j] = B_vfm[0];
                  mptr->By_vfm_ns[j] = B_vfm[1];
                  mptr->Bz_vfm_ns[j] = B_vfm[2];
                }
            }

        }
    }

  gsl_rng_free(rng_p);
  gsl_vector_free(g);

  for (i = 0; i < w->max_threads; ++i)
    mfield_free(mfield_array[i]);

  free(mfield_array);

  return 0;
}

/* compute B(r,theta,phi) using Gauss coefficients 'g' */
static int
mfield_synth_calc(const double t, const double r, const double theta, const double phi, const gsl_vector * g,
                  double B[3], mfield_workspace *w)
{
  int s = 0;
  gsl_bspline2_workspace *gauss_spline_p = w->gauss_spline_workspace_p[0];
  const size_t ncontrol = gsl_bspline2_ncontrol(gauss_spline_p);
  const size_t nmax_core = w->nmax_core;
  const size_t nmax = w->nmax;
  const double ratio = w->R / r;
  const double sint = sin(theta);
  const double cost = cos(theta);
  const double tyr = epoch2year(t);
  double rterm = ratio * ratio;
  double *Plm = w->Plm;
  double *dPlm = w->dPlm;
  size_t n;

  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cost, Plm, dPlm);

  for (n = 0; n <= nmax; ++n)
    {
      w->cosmphi[n] = cos(n * phi);
      w->sinmphi[n] = sin(n * phi);
    }

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      int m;

      /* (R/r)^{n+2} */
      rterm *= ratio;

      for (m = 0; m <= ni; ++m)
        {
          double c = w->cosmphi[m];
          double s = w->sinmphi[m];
          size_t pidx = gsl_sf_legendre_array_index(n, m);
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm, hnm = 0.0;

          if (n <= nmax_core)
            {
              gsl_vector_const_view v1 = gsl_vector_const_subvector_with_stride(g, cidx, w->nnm_core, ncontrol);
              gsl_bspline2_eval(tyr, &v1.vector, &gnm, gauss_spline_p);
            }
          else
            {
              gnm = gsl_vector_get(g, cidx + (w->p_core - w->nnm_core));
            }

          if (m > 0)
            {
              cidx = mfield_coeff_nmidx(n, -m);

              if (n <= nmax_core)
                {
                  gsl_vector_const_view v2 = gsl_vector_const_subvector_with_stride(g, cidx, w->nnm_core, ncontrol);
                  gsl_bspline2_eval(tyr, &v2.vector, &hnm, gauss_spline_p);
                }
              else
                {
                  hnm = gsl_vector_get(g, cidx + (w->p_core - w->nnm_core));
                }
            }

          B[0] += rterm * (gnm * c + hnm * s) * dPlm[pidx];
          B[1] += rterm / sint * m * (gnm * s - hnm * c) * Plm[pidx];
          B[2] -= (n + 1.0) * rterm * (gnm * c + hnm * s) * Plm[pidx];
        }
    }

  return s;
}

/* compute d/dt B(r,theta,phi) using Gauss coefficients 'g' */
static int
mfield_synth_calc_dBdt(const double t, const double r, const double theta, const double phi, const gsl_vector * g,
                       double dBdt[3], mfield_workspace *w)
{
  int s = 0;
  gsl_bspline2_workspace *gauss_spline_p = w->gauss_spline_workspace_p[0];
  const size_t ncontrol = gsl_bspline2_ncontrol(gauss_spline_p);
  const size_t nmax_core = w->nmax_core;
  const size_t nmax = w->nmax;
  const double ratio = w->R / r;
  const double sint = sin(theta);
  const double cost = cos(theta);
  const double tyr = epoch2year(t);
  double rterm = ratio * ratio;
  double *Plm = w->Plm;
  double *dPlm = w->dPlm;
  size_t n;

  dBdt[0] = 0.0;
  dBdt[1] = 0.0;
  dBdt[2] = 0.0;

  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cost, Plm, dPlm);

  for (n = 0; n <= nmax; ++n)
    {
      w->cosmphi[n] = cos(n * phi);
      w->sinmphi[n] = sin(n * phi);
    }

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      int m;

      /* (R/r)^{n+2} */
      rterm *= ratio;

      for (m = 0; m <= ni; ++m)
        {
          double c = w->cosmphi[m];
          double s = w->sinmphi[m];
          size_t pidx = gsl_sf_legendre_array_index(n, m);
          size_t cidx = mfield_coeff_nmidx(n, m);
          double dgnm, dhnm = 0.0;

          if (n <= nmax_core)
            {
              gsl_vector_const_view v1 = gsl_vector_const_subvector_with_stride(g, cidx, w->nnm_core, ncontrol);
              gsl_bspline2_eval_deriv(tyr, &v1.vector, 1, &dgnm, gauss_spline_p);
            }

          if (m > 0)
            {
              cidx = mfield_coeff_nmidx(n, -m);

              if (n <= nmax_core)
                {
                  gsl_vector_const_view v2 = gsl_vector_const_subvector_with_stride(g, cidx, w->nnm_core, ncontrol);
                  gsl_bspline2_eval_deriv(tyr, &v2.vector, 1, &dhnm, gauss_spline_p);
                }
            }

          dBdt[0] += rterm * (dgnm * c + dhnm * s) * dPlm[pidx];
          dBdt[1] += rterm / sint * m * (dgnm * s - dhnm * c) * Plm[pidx];
          dBdt[2] -= (n + 1.0) * rterm * (dgnm * c + dhnm * s) * Plm[pidx];
        }
    }

  return s;
}
