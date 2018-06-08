/*
 * mfield_residual.c
 *
 * Routines for computing and printing model residuals
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

#include "euler.h"
#include "magdata.h"
#include "mfield.h"
#include "mfield_residual.h"

static int mfield_residual_print_stat(const char *component_str, const gsl_rstat_workspace *rstat_p);
static int mfield_residual_print_satellite(const char *prefix, const size_t iter, const size_t nsource,
                                           size_t * index, magdata * mptr, mfield_workspace * w);
static int mfield_residual_print_observatory(const char *prefix, const size_t iter,
                                             size_t * index, magdata * mptr, mfield_workspace * w);

int
mfield_residual_print(const char *prefix, const size_t iter, mfield_workspace *w)
{
  int s = 0;
  size_t i;
  mfield_data_workspace *data_p = w->data_workspace_p;
  size_t idx = 0;

  fprintf(stderr, "\n");

  for (i = 0; i < data_p->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, data_p);

      if (mptr->global_flags & MAGDATA_GLOBFLG_SATELLITE)
        mfield_residual_print_satellite(prefix, iter, i, &idx, mptr, w);
      else if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
        mfield_residual_print_observatory(prefix, iter, &idx, mptr, w);
    }

  assert(idx == w->nres);

  return s;
}

static int
mfield_residual_print_stat(const char *component_str, const gsl_rstat_workspace *rstat_p)
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

/*
mfield_residual_print_satellite()
  Print residuals for satellite data

Inputs: prefix  - directory prefix for output files
        iter    - robust iteration number
        nsource - satellite number
        index   - (input/output)
        mptr    - magdata
        w       - mfield workspace
*/

static int
mfield_residual_print_satellite(const char *prefix, const size_t iter, const size_t nsource,
                                size_t * index, magdata * mptr, mfield_workspace * w)
{
  int s = 0;
  const char *fmtstr = "%ld %.8f %.4f %.4f %.4f %.4f %.3e %.3e %.3e %.4f %.4f %.4f %.4f\n";
  const char *fmtstr_F = "%ld %.8f %.4f %.4f %.4f %.4f %.3e %.3e %.3e %.4f %.4f %.4f\n";
  const char *fmtstr_grad = "%ld %.8f %.4f %.4f %.4f %.4f %.3e %.3e %.3e %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n";
  const double qdlat_cutoff = w->params.qdlat_fit_cutoff; /* cutoff latitude for high/low statistics */
  const size_t n = 12; /* number of components to write to disk */
  FILE *fp[12];
  size_t idx = *index;
  size_t j, k;
  char buf[2048];
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

  sprintf(buf, "%s/res%zu_X_iter%zu.dat", prefix, nsource, iter);
  fp[0] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_Y_iter%zu.dat", prefix, nsource, iter);
  fp[1] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_Z_iter%zu.dat", prefix, nsource, iter);
  fp[2] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_F_iter%zu.dat", prefix, nsource, iter);
  fp[3] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_DX_NS_iter%zu.dat", prefix, nsource, iter);
  fp[4] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_DY_NS_iter%zu.dat", prefix, nsource, iter);
  fp[5] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_DZ_NS_iter%zu.dat", prefix, nsource, iter);
  fp[6] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_DF_NS_iter%zu.dat", prefix, nsource, iter);
  fp[7] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_DX_EW_iter%zu.dat", prefix, nsource, iter);
  fp[8] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_DY_EW_iter%zu.dat", prefix, nsource, iter);
  fp[9] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_DZ_EW_iter%zu.dat", prefix, nsource, iter);
  fp[10] = fopen(buf, "w");

  sprintf(buf, "%s/res%zu_DF_EW_iter%zu.dat", prefix, nsource, iter);
  fp[11] = fopen(buf, "w");

  /* header line */
  fprintf(fp[0], "# X vector residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[1], "# Y vector residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[2], "# Z vector residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[3], "# F scalar residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[4], "# DX gradient (N/S) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[5], "# DY gradient (N/S) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[6], "# DZ gradient (N/S) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[7], "# DZ gradient (N/S) scalar residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[8], "# DX gradient (E/W) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[9], "# DY gradient (E/W) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[10], "# DZ gradient (E/W) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);
  fprintf(fp[11], "# DF gradient (E/W) scalar residuals for MF modeling (satellite %zu, iteration %zu)\n", nsource, iter);

  for (j = 0; j < n; ++j)
    {
      k = 1;
      fprintf(fp[j], "# Field %zu: timestamp (UT seconds since 1970-01-01)\n", k++);
      fprintf(fp[j], "# Field %zu: time (decimal year)\n", k++);
      fprintf(fp[j], "# Field %zu: longitude (degrees)\n", k++);
      fprintf(fp[j], "# Field %zu: geocentric latitude (degrees)\n", k++);
      fprintf(fp[j], "# Field %zu: QD latitude (degrees)\n", k++);
      fprintf(fp[j], "# Field %zu: geocentric radius (km)\n", k++);
      fprintf(fp[j], "# Field %zu: spatial weight factor\n", k++);
      fprintf(fp[j], "# Field %zu: robust weight factor\n", k++);
      fprintf(fp[j], "# Field %zu: total weight factor\n", k++);
    }

  fprintf(fp[0], "# Field %zu: X vector measurement (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y vector measurement (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z vector measurement (nT)\n", k);
  fprintf(fp[3], "# Field %zu: F scalar measurement (nT)\n", k);
  fprintf(fp[4], "# Field %zu: X vector measurement (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y vector measurement (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z vector measurement (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F scalar measurement (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X vector measurement (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y vector measurement (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z vector measurement (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F scalar measurement (nT)\n", k);
  ++k;

  fprintf(fp[0], "# Field %zu: X a priori model (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y a priori model (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z a priori model (nT)\n", k);
  fprintf(fp[3], "# Field %zu: | B_prior + B_fitted | (nT)\n", k);
  fprintf(fp[4], "# Field %zu: X a priori model (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y a priori model (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z a priori model (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F a priori model (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X a priori model (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y a priori model (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z a priori model (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F a priori model (nT)\n", k);
  ++k;

  fprintf(fp[0], "# Field %zu: X fitted model (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y fitted model (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z fitted model (nT)\n", k);
  fprintf(fp[3], "# Field %zu: scalar residual (nT)\n", k);
  fprintf(fp[4], "# Field %zu: X fitted model (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y fitted model (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z fitted model (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F fitted model (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X fitted model (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y fitted model (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z fitted model (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F fitted model (nT)\n", k);
  ++k;

  fprintf(fp[0], "# Field %zu: X residual (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y residual (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z residual (nT)\n", k);
  fprintf(fp[4], "# Field %zu: X vector measurement at N/S gradient point (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y vector measurement at N/S gradient point (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z vector measurement at N/S gradient point (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F scalar measurement at N/S gradient point (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X vector measurement at E/W gradient point (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y vector measurement at E/W gradient point (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z vector measurement at E/W gradient point (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F scalar measurement at E/W gradient point (nT)\n", k);
  ++k;

  fprintf(fp[4], "# Field %zu: X a priori model at N/S gradient point (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y a priori model at N/S gradient point (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z a priori model at N/S gradient point (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F a priori model at N/S gradient point (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X a priori model at E/W gradient point (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y a priori model at E/W gradient point (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z a priori model at E/W gradient point (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F a priori model at E/W gradient point (nT)\n", k);
  ++k;

  fprintf(fp[4], "# Field %zu: X fitted model at N/S gradient point (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y fitted model at N/S gradient point (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z fitted model at N/S gradient point (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F fitted model at N/S gradient point (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X fitted model at E/W gradient point (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y fitted model at E/W gradient point (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z fitted model at E/W gradient point (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F fitted model at E/W gradient point (nT)\n", k);
  ++k;

  fprintf(fp[4], "# Field %zu: DX N/S residual (nT)\n", k);
  fprintf(fp[5], "# Field %zu: DY N/S residual (nT)\n", k);
  fprintf(fp[6], "# Field %zu: DZ N/S residual (nT)\n", k);
  fprintf(fp[8], "# Field %zu: DX E/W residual (nT)\n", k);
  fprintf(fp[9], "# Field %zu: DY E/W residual (nT)\n", k);
  fprintf(fp[10], "# Field %zu: DZ E/W residual (nT)\n", k);
  ++k;

  for (j = 0; j < mptr->n; ++j)
    {
      double t = satdata_epoch2year(mptr->t[j]);
      time_t unix_time = satdata_epoch2timet(mptr->t[j]);
      double r = mptr->r[j];
      double theta = mptr->theta[j];
      double phi = mptr->phi[j];
      double lat = 90.0 - theta * 180.0 / M_PI;
      double qdlat = mptr->qdlat[j];
      double B_nec[4], B_nec_grad[4];     /* observations in NEC */
      double B_model[4], B_grad_model[4]; /* a priori models */
      double B_fit[4], B_grad_fit[4];     /* fitted field model */
      double res[4], res_grad[4];         /* residuals */

      if (MAGDATA_Discarded(mptr->flags[j]))
        continue;

      if (MAGDATA_FitMF(mptr->flags[j]))
        {
          /* evaluate internal field models */
          mfield_eval(mptr->t[j], r, theta, phi, B_fit, w);

          if (MAGDATA_ExistDX_NS(mptr->flags[j]) || MAGDATA_ExistDY_NS(mptr->flags[j]) ||
              MAGDATA_ExistDZ_NS(mptr->flags[j]) || MAGDATA_ExistDF_NS(mptr->flags[j]) ||
              MAGDATA_ExistDX_EW(mptr->flags[j]) || MAGDATA_ExistDY_EW(mptr->flags[j]) ||
              MAGDATA_ExistDZ_EW(mptr->flags[j]) || MAGDATA_ExistDF_EW(mptr->flags[j]))
            {
              mfield_eval(mptr->t_ns[j], mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j], B_grad_fit, w);
            }

          B_model[0] = mptr->Bx_model[j];
          B_model[1] = mptr->By_model[j];
          B_model[2] = mptr->Bz_model[j];
          B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

          B_nec_grad[0] = mptr->Bx_nec_ns[j];
          B_nec_grad[1] = mptr->By_nec_ns[j];
          B_nec_grad[2] = mptr->Bz_nec_ns[j];
          B_nec_grad[3] = mptr->F_ns[j];

          B_grad_model[0] = mptr->Bx_model_ns[j];
          B_grad_model[1] = mptr->By_model_ns[j];
          B_grad_model[2] = mptr->Bz_model_ns[j];
          B_grad_model[3] = gsl_hypot3(B_grad_model[0], B_grad_model[1], B_grad_model[2]);

          B_nec[3] = mptr->F[j];

          if (w->params.fit_euler && mptr->global_flags & MAGDATA_GLOBFLG_EULER)
            {
              size_t euler_idx = mfield_euler_idx(nsource, mptr->t[j], w);
              double alpha = gsl_vector_get(w->c, euler_idx);
              double beta = gsl_vector_get(w->c, euler_idx + 1);
              double gamma = gsl_vector_get(w->c, euler_idx + 2);
              double *q = &(mptr->q[4*j]);

              B_nec[0] = mptr->Bx_vfm[j];
              B_nec[1] = mptr->By_vfm[j];
              B_nec[2] = mptr->Bz_vfm[j];

              /* rotate to NEC with computed Euler angles */
              euler_vfm2nec(mptr->euler_flags, alpha, beta, gamma, q, B_nec, B_nec);
            }
          else
            {
              B_nec[0] = mptr->Bx_nec[j];
              B_nec[1] = mptr->By_nec[j];
              B_nec[2] = mptr->Bz_nec[j];
            }

          /* calculate residuals */
          for (k = 0; k < 3; ++k)
            {
              res[k] = B_nec[k] - B_model[k] - B_fit[k];
              res_grad[k] = B_nec_grad[k] - B_grad_model[k] - B_grad_fit[k];
            }

          res[3] = B_nec[3] - gsl_hypot3(B_model[0] + B_fit[0], B_model[1] + B_fit[1], B_model[2] + B_fit[2]);
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

              fprintf(fp[0], fmtstr, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[0], B_model[0], B_fit[0], res[0]);
              gsl_rstat_add(B_nec[0] - B_model[0] - B_fit[0], rstat_x);
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
              fprintf(fp[1], fmtstr, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[1], B_model[1], B_fit[1], res[1]);
              gsl_rstat_add(B_nec[1] - B_model[1] - B_fit[1], rstat_y);
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
              fprintf(fp[2], fmtstr, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[2], B_model[2], B_fit[2], res[2]);

              gsl_rstat_add(B_nec[2] - B_model[2] - B_fit[2], rstat_z);

              if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                gsl_rstat_add(B_nec[2] - B_model[2] - B_fit[2], rstat_lowz);
              else
                gsl_rstat_add(B_nec[2] - B_model[2] - B_fit[2], rstat_highz);
            }

          ++idx;
        }

      if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[3], fmtstr_F, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[3], B_nec[3] - res[3], res[3]);
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
              fprintf(fp[4], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[0], B_model[0], B_fit[0], B_nec_grad[0], B_grad_model[0], B_grad_fit[0], res[0] - res_grad[0]);
              gsl_rstat_add(res[0] - res_grad[0], rstat_dx_ns);
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
              fprintf(fp[5], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[1], B_model[1], B_fit[1], B_nec_grad[1], B_grad_model[1], B_grad_fit[1], res[1] - res_grad[1]);
              gsl_rstat_add(res[1] - res_grad[1], rstat_dy_ns);
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
              fprintf(fp[6], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[2], B_model[2], B_fit[2], B_nec_grad[2], B_grad_model[2], B_grad_fit[2], res[2] - res_grad[2]);

              if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                gsl_rstat_add(res[2] - res_grad[2], rstat_low_dz_ns);
              else
                gsl_rstat_add(res[2] - res_grad[2], rstat_high_dz_ns);
            }

          ++idx;
        }

      if (MAGDATA_ExistDF_NS(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);
          fprintf(fp[7], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[3], B_model[3], B_fit[3], B_nec_grad[3], B_grad_model[3], B_grad_fit[3]);
          gsl_rstat_add(B_nec[3] - B_model[3] - B_fit[3] - (B_nec_grad[3] - B_grad_model[3] - B_grad_fit[3]), rstat_df_ns);

          ++idx;
        }

      if (MAGDATA_ExistDX_EW(mptr->flags[j]))
        {
          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              double ws = gsl_vector_get(w->wts_spatial, idx);
              double wr = gsl_vector_get(w->wts_robust, idx);
              double wf = gsl_vector_get(w->wts_final, idx);
              fprintf(fp[8], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[0], B_model[0], B_fit[0], B_nec_grad[0], B_grad_model[0], B_grad_fit[0], res[0] - res_grad[0]);
              gsl_rstat_add(res[0] - res_grad[0], rstat_dx_ew);
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
              fprintf(fp[9], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[1], B_model[1], B_fit[1], B_nec_grad[1], B_grad_model[1], B_grad_fit[1], res[1] - res_grad[1]);
              gsl_rstat_add(res[1] - res_grad[1], rstat_dy_ew);
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
              fprintf(fp[10], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[2], B_model[2], B_fit[2], B_nec_grad[2], B_grad_model[2], B_grad_fit[2], res[2] - res_grad[2]);

              if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                gsl_rstat_add(res[2] - res_grad[2], rstat_low_dz_ew);
              else
                gsl_rstat_add(res[2] - res_grad[2], rstat_high_dz_ew);
            }

          ++idx;
        }

      if (MAGDATA_ExistDF_EW(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);
          fprintf(fp[11], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, ws, wr, wf, B_nec[3], B_model[3], B_fit[3], B_nec_grad[3], B_grad_model[3], B_grad_fit[3]);
          gsl_rstat_add(B_nec[3] - B_model[3] - B_fit[3] - (B_nec_grad[3] - B_grad_model[3] - B_grad_fit[3]), rstat_df_ew);

          ++idx;
        }
    }

  fprintf(stderr, "=== FIT STATISTICS SATELLITE %zu ===\n", nsource);

  /* print header */
  mfield_residual_print_stat(NULL, NULL);

  mfield_residual_print_stat("X", rstat_x);
  mfield_residual_print_stat("Y", rstat_y);
  mfield_residual_print_stat("Z", rstat_z);
  mfield_residual_print_stat("F", rstat_f);
  mfield_residual_print_stat("low Z", rstat_lowz);
  mfield_residual_print_stat("high Z", rstat_highz);
  mfield_residual_print_stat("low F", rstat_lowf);
  mfield_residual_print_stat("high F", rstat_highf);

  mfield_residual_print_stat("N/S DX", rstat_dx_ns);
  mfield_residual_print_stat("N/S DY", rstat_dy_ns);
  mfield_residual_print_stat("low N/S DZ", rstat_low_dz_ns);
  mfield_residual_print_stat("high N/S DZ", rstat_high_dz_ns);
  mfield_residual_print_stat("N/S DF", rstat_df_ns);

  mfield_residual_print_stat("E/W DX", rstat_dx_ew);
  mfield_residual_print_stat("E/W DY", rstat_dy_ew);
  mfield_residual_print_stat("low E/W DZ", rstat_low_dz_ew);
  mfield_residual_print_stat("high E/W DZ", rstat_high_dz_ew);
  mfield_residual_print_stat("E/W DF", rstat_df_ew);

  *index = idx;

  for (j = 0; j < n; ++j)
    fclose(fp[j]);

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

/*
mfield_residual_print_observatory()
  Print residuals for observatory data

Inputs: prefix  - directory prefix for output files
        iter    - robust iteration number
        index   - (input/output)
        mptr    - magdata
        w       - mfield workspace
*/

static int
mfield_residual_print_observatory(const char *prefix, const size_t iter,
                                  size_t * index, magdata * mptr, mfield_workspace * w)
{
  int s = 0;
  const char *fmtstr = "%ld %8.4f %6.3f %6.3f %6.3f %10.4f %10.4f %10.4f\n";
  const size_t n = 3; /* number of components to write to disk */
  const double r = mptr->r[0];
  const double theta = mptr->theta[0];
  const double phi = mptr->phi[0];
  FILE *fp[3];
  size_t idx = *index;
  size_t j, k;
  char buf[2048];
  gsl_rstat_workspace *rstat_dxdt = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_dydt = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_dzdt = gsl_rstat_alloc();

  sprintf(buf, "%s/obs/res_%s_DXDT_iter%zu.dat", prefix, mptr->name, iter);
  fp[0] = fopen(buf, "w");

  sprintf(buf, "%s/obs/res_%s_DYDT_iter%zu.dat", prefix, mptr->name, iter);
  fp[1] = fopen(buf, "w");

  sprintf(buf, "%s/obs/res_%s_DZDT_iter%zu.dat", prefix, mptr->name, iter);
  fp[2] = fopen(buf, "w");

  for (j = 0; j < n; ++j)
    {
      if (fp[j] == NULL)
        {
          fprintf(stderr, "mfield_residual_print_observatory: fp[%zu] is NULL\n", j);
          return -1;
        }
    }

  /* print header */
  fprintf(fp[0], "# %s observatory\n# dX/dt vector data for MF modeling\n", mptr->name);
  fprintf(fp[1], "# %s observatory\n# dY/dt vector data for MF modeling\n", mptr->name);
  fprintf(fp[2], "# %s observatory\n# dZ/dt vector data for MF modeling\n", mptr->name);

  for (j = 0; j < n; ++j)
    {
      fprintf(fp[j], "# Radius:    %.4f [km]\n", r);
      fprintf(fp[j], "# Longitude: %.4f [deg]\n", phi * 180.0 / M_PI);
      fprintf(fp[j], "# Latitude:  %.4f [deg]\n", 90.0 - theta * 180.0 / M_PI);

      k = 1;
      fprintf(fp[j], "# Field %zu: timestamp (UT seconds since 1970-01-01)\n", k++);
      fprintf(fp[j], "# Field %zu: QD latitude (degrees)\n", k++);
      fprintf(fp[j], "# Field %zu: spatial weight factor\n", k++);
      fprintf(fp[j], "# Field %zu: robust weight factor\n", k++);
      fprintf(fp[j], "# Field %zu: total weight factor\n", k++);
    }

  fprintf(fp[0], "# Field %zu: dX/dt vector measurement (nT/year)\n", k);
  fprintf(fp[1], "# Field %zu: dY/dt vector measurement (nT/year)\n", k);
  fprintf(fp[2], "# Field %zu: dZ/dt vector measurement (nT/year)\n", k);
  ++k;

  fprintf(fp[0], "# Field %zu: dX/dt fitted model (nT/year)\n", k);
  fprintf(fp[1], "# Field %zu: dY/dt fitted model (nT/year)\n", k);
  fprintf(fp[2], "# Field %zu: dZ/dt fitted model (nT/year)\n", k);
  ++k;

  fprintf(fp[0], "# Field %zu: dX/dt residual (nT/year)\n", k);
  fprintf(fp[1], "# Field %zu: dY/dt residual (nT/year)\n", k);
  fprintf(fp[2], "# Field %zu: dZ/dt residual (nT/year)\n", k);
  ++k;

  for (j = 0; j < mptr->n; ++j)
    {
      time_t unix_time = satdata_epoch2timet(mptr->t[j]);
      double dBdt_fit[4], res_dBdt[3];

      if (MAGDATA_Discarded(mptr->flags[j]))
        continue;

      if (!MAGDATA_FitMF(mptr->flags[j]))
        continue;

      /* evaluate internal field models */
      mfield_eval_dBdt(mptr->t[j], r, theta, phi, dBdt_fit, w);

      res_dBdt[0] = mptr->dXdt_nec[j] - dBdt_fit[0];
      res_dBdt[1] = mptr->dYdt_nec[j] - dBdt_fit[1];
      res_dBdt[2] = mptr->dZdt_nec[j] - dBdt_fit[2];

      if (MAGDATA_ExistDXDT(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[0], fmtstr, unix_time, mptr->qdlat[j], ws, wr, wf, mptr->dXdt_nec[j], dBdt_fit[0], res_dBdt[0]);

          ++idx;
        }

      if (MAGDATA_ExistDYDT(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[1], fmtstr, unix_time, mptr->qdlat[j], ws, wr, wf, mptr->dYdt_nec[j], dBdt_fit[1], res_dBdt[1]);

          ++idx;
        }

      if (MAGDATA_ExistDZDT(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[2], fmtstr, unix_time, mptr->qdlat[j], ws, wr, wf, mptr->dZdt_nec[j], dBdt_fit[2], res_dBdt[2]);

          ++idx;
        }
    }

  *index = idx;

  for (j = 0; j < n; ++j)
    fclose(fp[j]);

  gsl_rstat_free(rstat_dxdt);
  gsl_rstat_free(rstat_dydt);
  gsl_rstat_free(rstat_dzdt);

  return s;
}
