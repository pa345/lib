/*
 * invert_residual.c
 *
 * Routines for computing and printing model residuals
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_magdata.h>

#include "invert.h"
#include "invert_residual.h"
#include "invert_multifit.c"

static int invert_residual_print_stat(const char *component_str, const gsl_rstat_workspace *rstat_p);
static int invert_residual_print_satellite(const char *prefix, const size_t iter, const size_t nsource,
                                           const gsl_vector * f, size_t * index, magdata * mptr, invert_workspace * w);
static int invert_residual_print_observatory(const char *prefix, const size_t iter, const size_t mptr_index,
                                             size_t * index, gsl_rstat_workspace ** rstat_X,
                                             gsl_rstat_workspace ** rstat_Y, gsl_rstat_workspace ** rstat_Z,
                                             invert_workspace * w);
static int invert_residual_print_observatory_SV(const char *prefix, const size_t iter,
                                                size_t * index, magdata * mptr,
                                                gsl_rstat_workspace ** rstat_DXDT,
                                                gsl_rstat_workspace ** rstat_DYDT,
                                                gsl_rstat_workspace ** rstat_DZDT,
                                                invert_workspace * w);

int
invert_residual_print(const char *prefix, const size_t iter, invert_workspace *w)
{
  int s = 0;
  size_t i;
  const size_t n = w->nres_tot;   /* number of residuals */
  invert_data_workspace *data_p = w->data_workspace_p;
  gsl_vector * f = gsl_vector_alloc(n);
  gsl_rstat_workspace **rstat_X = malloc(3 * sizeof(gsl_rstat_workspace));
  gsl_rstat_workspace **rstat_Y = malloc(3 * sizeof(gsl_rstat_workspace));
  gsl_rstat_workspace **rstat_Z = malloc(3 * sizeof(gsl_rstat_workspace));
  gsl_rstat_workspace **rstat_DXDT = malloc(3 * sizeof(gsl_rstat_workspace));
  gsl_rstat_workspace **rstat_DYDT = malloc(3 * sizeof(gsl_rstat_workspace));
  gsl_rstat_workspace **rstat_DZDT = malloc(3 * sizeof(gsl_rstat_workspace));
  size_t idx = 0;
  int print_obs_stat = 0;
  struct timeval tv0, tv1;

  for (i = 0; i < 3; ++i)
    {
      rstat_X[i] = gsl_rstat_alloc();
      rstat_Y[i] = gsl_rstat_alloc();
      rstat_Z[i] = gsl_rstat_alloc();
      rstat_DXDT[i] = gsl_rstat_alloc();
      rstat_DYDT[i] = gsl_rstat_alloc();
      rstat_DZDT[i] = gsl_rstat_alloc();
    }

  fprintf(stderr, "\n");

  fprintf(stderr, "invert_residual_print: computing residual vector...");
  gettimeofday(&tv0, NULL);
  invert_calc_f(w->c, w, f);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, || f || = %.12e)\n", time_diff(tv0, tv1), gsl_blas_dnrm2(f));

  for (i = 0; i < data_p->nsources; ++i)
    {
      magdata *mptr = invert_data_ptr(i, data_p);

      if (mptr->global_flags & MAGDATA_GLOBFLG_SATELLITE)
        invert_residual_print_satellite(prefix, iter, i, f, &idx, mptr, w);
      else if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
        invert_residual_print_observatory(prefix, iter, i, &idx, rstat_X, rstat_Y, rstat_Z, w);
#if 0/*XXX*/
      else if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY_SV)
        invert_residual_print_observatory_SV(prefix, iter, &idx, mptr, rstat_DXDT, rstat_DYDT, rstat_DZDT, w);
#endif

      if (mptr->global_flags & (MAGDATA_GLOBFLG_OBSERVATORY | MAGDATA_GLOBFLG_OBSERVATORY_SV))
        print_obs_stat = 1;
    }

  assert(idx == w->nres);

  if (print_obs_stat)
    {
      fprintf(stderr, "=== FIT STATISTICS OBSERVATORY (NORTH POLE) ===\n");

      /* print header */
      invert_residual_print_stat(NULL, NULL);

      invert_residual_print_stat("X", rstat_X[0]);
      invert_residual_print_stat("Y", rstat_Y[0]);
      invert_residual_print_stat("Z", rstat_Z[0]);
      invert_residual_print_stat("dX/dt", rstat_DXDT[0]);
      invert_residual_print_stat("dY/dt", rstat_DYDT[0]);
      invert_residual_print_stat("dZ/dt", rstat_DZDT[0]);

      fprintf(stderr, "=== FIT STATISTICS OBSERVATORY (SOUTH POLE) ===\n");

      /* print header */
      invert_residual_print_stat(NULL, NULL);

      invert_residual_print_stat("X", rstat_X[1]);
      invert_residual_print_stat("Y", rstat_Y[1]);
      invert_residual_print_stat("Z", rstat_Z[1]);
      invert_residual_print_stat("dX/dt", rstat_DXDT[1]);
      invert_residual_print_stat("dY/dt", rstat_DYDT[1]);
      invert_residual_print_stat("dZ/dt", rstat_DZDT[1]);

      fprintf(stderr, "=== FIT STATISTICS OBSERVATORY (MIDDLE LATITUDES) ===\n");

      /* print header */
      invert_residual_print_stat(NULL, NULL);

      invert_residual_print_stat("X", rstat_X[2]);
      invert_residual_print_stat("Y", rstat_Y[2]);
      invert_residual_print_stat("Z", rstat_Z[2]);
      invert_residual_print_stat("dX/dt", rstat_DXDT[2]);
      invert_residual_print_stat("dY/dt", rstat_DYDT[2]);
      invert_residual_print_stat("dZ/dt", rstat_DZDT[2]);
    }

  for (i = 0; i < 3; ++i)
    {
      gsl_rstat_free(rstat_X[i]);
      gsl_rstat_free(rstat_Y[i]);
      gsl_rstat_free(rstat_Z[i]);
      gsl_rstat_free(rstat_DXDT[i]);
      gsl_rstat_free(rstat_DYDT[i]);
      gsl_rstat_free(rstat_DZDT[i]);
    }

  free(rstat_X);
  free(rstat_Y);
  free(rstat_Z);
  free(rstat_DXDT);
  free(rstat_DYDT);
  free(rstat_DZDT);
  gsl_vector_free(f);

  return s;
}

static int
invert_residual_print_stat(const char *component_str, const gsl_rstat_workspace *rstat_p)
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
invert_residual_print_satellite()
  Print residuals for satellite data

Inputs: prefix  - directory prefix for output files
        iter    - robust iteration number
        nsource - satellite number
        f       - residual vector
        index   - (input/output)
        mptr    - magdata
        w       - invert workspace
*/

static int
invert_residual_print_satellite(const char *prefix, const size_t iter, const size_t nsource,
                                const gsl_vector * f, size_t * index, magdata * mptr, invert_workspace * w)
{
  int s = 0;
  const char *fmtstr = "%ld %.4f %.4f %.4f %.4f %.4f %.3e %.3e %.3e %.4f %.4f %.4f %.4f\n";
  const char *fmtstr_F = "%ld %.4f %.4f %.4f %.4f %.4f %.3e %.3e %.3e %.4f %.4f %.4f\n";
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
      fprintf(fp[j], "# Field %zu: local time (hours)\n", k++);
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
      time_t unix_time = satdata_epoch2timet(mptr->t[j]);
      double lt = get_localtime(unix_time, mptr->phi[j]);
      double r = mptr->r[j];
      double theta = mptr->theta[j];
      double phi = mptr->phi[j] * 180.0 / M_PI;
      double lat = 90.0 - theta * 180.0 / M_PI;
      double qdlat = mptr->qdlat[j];
      double B_nec[4];     /* observations in NEC */
      double B_model[4];   /* a priori models */

      B_model[0] = mptr->Bx_model[j];
      B_model[1] = mptr->By_model[j];
      B_model[2] = mptr->Bz_model[j];
      B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

      B_nec[0] = mptr->Bx_nec[j];
      B_nec[1] = mptr->By_nec[j];
      B_nec[2] = mptr->Bz_nec[j];
      B_nec[3] = mptr->F[j];

      /*XXXres[3] = B_nec[3] - gsl_hypot3(B_model[0] + B_fit[0], B_model[1] + B_fit[1], B_model[2] + B_fit[2]);*/

      if ((j > 0) && (mptr->flags[j] & MAGDATA_FLG_TRACK_START))
        {
          for (k = 0; k < n; ++k)
            fprintf(fp[k], "\n\n");
        }

      if (MAGDATA_ExistX(mptr->flags[j]))
        {
          double X_res = gsl_vector_get(f, idx);
          double X_fit = B_nec[0] - B_model[0] - X_res;
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[0], fmtstr, unix_time, lt, phi, lat, qdlat, r, ws, wr, wf, B_nec[0], B_model[0], X_fit, X_res);
          gsl_rstat_add(X_res, rstat_x);

          ++idx;
        }

      if (MAGDATA_ExistY(mptr->flags[j]))
        {
          double Y_res = gsl_vector_get(f, idx);
          double Y_fit = B_nec[1] - B_model[1] - Y_res;
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[1], fmtstr, unix_time, lt, phi, lat, qdlat, r, ws, wr, wf, B_nec[1], B_model[1], Y_fit, Y_res);
          gsl_rstat_add(Y_res, rstat_y);

          ++idx;
        }

      if (MAGDATA_ExistZ(mptr->flags[j]))
        {
          double Z_res = gsl_vector_get(f, idx);
          double Z_fit = B_nec[2] - B_model[2] - Z_res;
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[2], fmtstr, unix_time, lt, phi, lat, qdlat, r, ws, wr, wf, B_nec[2], B_model[2], Z_fit, Z_res);

          gsl_rstat_add(Z_res, rstat_z);

          if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
            gsl_rstat_add(Z_res, rstat_lowz);
          else
            gsl_rstat_add(Z_res, rstat_highz);

          ++idx;
        }

      if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
        {
          double F_res = gsl_vector_get(f, idx);
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[3], fmtstr_F, unix_time, lt, phi, lat, qdlat, r, ws, wr, wf, B_nec[3], B_nec[3] - F_res, F_res);
          gsl_rstat_add(F_res, rstat_f);

          if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
            gsl_rstat_add(F_res, rstat_lowf);
          else
            gsl_rstat_add(F_res, rstat_highf);

          ++idx;
        }
    }

  fprintf(stderr, "=== FIT STATISTICS SATELLITE %zu ===\n", nsource);

  /* print header */
  invert_residual_print_stat(NULL, NULL);

  invert_residual_print_stat("X", rstat_x);
  invert_residual_print_stat("Y", rstat_y);
  invert_residual_print_stat("Z", rstat_z);
  invert_residual_print_stat("F", rstat_f);
  invert_residual_print_stat("low Z", rstat_lowz);
  invert_residual_print_stat("high Z", rstat_highz);
  invert_residual_print_stat("low F", rstat_lowf);
  invert_residual_print_stat("high F", rstat_highf);

  invert_residual_print_stat("N/S DX", rstat_dx_ns);
  invert_residual_print_stat("N/S DY", rstat_dy_ns);
  invert_residual_print_stat("low N/S DZ", rstat_low_dz_ns);
  invert_residual_print_stat("high N/S DZ", rstat_high_dz_ns);
  invert_residual_print_stat("N/S DF", rstat_df_ns);

  invert_residual_print_stat("E/W DX", rstat_dx_ew);
  invert_residual_print_stat("E/W DY", rstat_dy_ew);
  invert_residual_print_stat("low E/W DZ", rstat_low_dz_ew);
  invert_residual_print_stat("high E/W DZ", rstat_high_dz_ew);
  invert_residual_print_stat("E/W DF", rstat_df_ew);

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
invert_residual_print_observatory()
  Print residuals for observatory data

Inputs: prefix     - directory prefix for output files
        iter       - robust iteration number
        mptr_index - magdata index
        index      - (input/output)
        rstat_X    - array of 3 rstat workspaces for X component
                     index 0: north pole stats
                     index 1: south pole stats
                     index 2: mid latitude stats
        rstat_Y    - array of 3 rstat workspaces for Y component
                     index 0: north pole stats
                     index 1: south pole stats
                     index 2: mid latitude stats
        rstat_Z    - array of 3 rstat workspaces for Z component
                     index 0: north pole stats
                     index 1: south pole stats
                     index 2: mid latitude stats
        w          - invert workspace
*/

static int
invert_residual_print_observatory(const char *prefix, const size_t iter, const size_t mptr_index,
                                  size_t * index, gsl_rstat_workspace ** rstat_X,
                                  gsl_rstat_workspace ** rstat_Y, gsl_rstat_workspace ** rstat_Z,
                                  invert_workspace * w)
{
  int s = 0;
  const double qdlat_cutoff = w->params.qdlat_fit_cutoff; /* cutoff latitude for high/low statistics */
  const magdata *mptr = invert_data_ptr(mptr_index, w->data_workspace_p);
  const char *fmtstr = "%ld %8.4f %6.3f %6.3f %6.3f %10.4f %10.4f %10.4f %10.4f\n";
  const size_t n = 3; /* number of components to write to disk */
  const double r = mptr->r[0];
  const double theta = mptr->theta[0];
  const double phi = mptr->phi[0];
  FILE *fp[3];
  size_t idx = *index;
  size_t j, k;
  char buf[2048];

  sprintf(buf, "%s/obs/res_%s_X_iter%zu.dat", prefix, mptr->name, iter);
  fp[0] = fopen(buf, "w");

  sprintf(buf, "%s/obs/res_%s_Y_iter%zu.dat", prefix, mptr->name, iter);
  fp[1] = fopen(buf, "w");

  sprintf(buf, "%s/obs/res_%s_Z_iter%zu.dat", prefix, mptr->name, iter);
  fp[2] = fopen(buf, "w");

  for (j = 0; j < n; ++j)
    {
      if (fp[j] == NULL)
        {
          fprintf(stderr, "invert_residual_print_observatory: fp[%zu] is NULL\n", j);
          return -1;
        }
    }

  /* print header */
  fprintf(fp[0], "# %s observatory\n# X vector data for MF modeling\n", mptr->name);
  fprintf(fp[1], "# %s observatory\n# Y vector data for MF modeling\n", mptr->name);
  fprintf(fp[2], "# %s observatory\n# Z vector data for MF modeling\n", mptr->name);

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

  fprintf(fp[0], "# Field %zu: X vector measurement (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y vector measurement (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z vector measurement (nT)\n", k);
  ++k;

  fprintf(fp[0], "# Field %zu: X fitted model (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y fitted model (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z fitted model (nT)\n", k);
  ++k;

  fprintf(fp[0], "# Field %zu: X a priori model (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y a priori model (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z a priori model (nT)\n", k);
  ++k;

  fprintf(fp[0], "# Field %zu: X residual (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y residual (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z residual (nT)\n", k);
  ++k;

  for (j = 0; j < mptr->n; ++j)
    {
      time_t unix_time = satdata_epoch2timet(mptr->t[j]);
      double B_nec[3];
      double B_model[3], B_fit[4], res_B[3];
      size_t l, rstat_idx;

      if (MAGDATA_Discarded(mptr->flags[j]))
        continue;

      if (!MAGDATA_FitMF(mptr->flags[j]))
        continue;

      if (mptr->qdlat[j] > qdlat_cutoff)
        rstat_idx = 0; /* north pole latitudes */
      else if (mptr->qdlat[j] < -qdlat_cutoff)
        rstat_idx = 1; /* south pole latitudes */
      else
        rstat_idx = 2; /* mid latitudes */

      B_nec[0] = mptr->Bx_nec[j];
      B_nec[1] = mptr->By_nec[j];
      B_nec[2] = mptr->Bz_nec[j];

      B_model[0] = mptr->Bx_model[j];
      B_model[1] = mptr->By_model[j];
      B_model[2] = mptr->Bz_model[j];

      /* evaluate internal field models */
      invert_eval(mptr->t[j], r, theta, phi, B_fit, w);

      for (l = 0; l < 3; ++l)
        res_B[l] = B_nec[l] - B_fit[l] - B_model[l];

      if (MAGDATA_ExistX(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[0], fmtstr, unix_time, mptr->qdlat[j], ws, wr, wf, mptr->Bx_nec[j], B_fit[0], B_model[0], res_B[0]);
          gsl_rstat_add(res_B[0], rstat_X[rstat_idx]);

          ++idx;
        }

      if (MAGDATA_ExistY(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[1], fmtstr, unix_time, mptr->qdlat[j], ws, wr, wf, mptr->By_nec[j], B_fit[1], B_model[1], res_B[1]);
          gsl_rstat_add(res_B[1], rstat_Y[rstat_idx]);

          ++idx;
        }

      if (MAGDATA_ExistZ(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[2], fmtstr, unix_time, mptr->qdlat[j], ws, wr, wf, mptr->Bz_nec[j], B_fit[2], B_model[2], res_B[2]);
          gsl_rstat_add(res_B[2], rstat_Z[rstat_idx]);

          ++idx;
        }
    }

  *index = idx;

  for (j = 0; j < n; ++j)
    fclose(fp[j]);

  return s;
}

/*
invert_residual_print_observatory_SV()
  Print residuals for observatory data

Inputs: prefix  - directory prefix for output files
        iter    - robust iteration number
        index   - (input/output)
        mptr    - magdata
        rstat_DXDT - array of 3 rstat workspaces for dX/dt component
                     index 0: north pole stats
                     index 1: south pole stats
                     index 2: mid latitude stats
        rstat_DYDT - array of 3 rstat workspaces for dY/dt component
                     index 0: north pole stats
                     index 1: south pole stats
                     index 2: mid latitude stats
        rstat_DZDT - array of 3 rstat workspaces for dZ/dt component
                     index 0: north pole stats
                     index 1: south pole stats
                     index 2: mid latitude stats
        w       - invert workspace
*/

#if 0 /*XXX*/
static int
invert_residual_print_observatory_SV(const char *prefix, const size_t iter,
                                     size_t * index, magdata * mptr,
                                     gsl_rstat_workspace ** rstat_DXDT,
                                     gsl_rstat_workspace ** rstat_DYDT,
                                     gsl_rstat_workspace ** rstat_DZDT,
                                     invert_workspace * w)
{
  int s = 0;
  const double qdlat_cutoff = w->params.qdlat_fit_cutoff; /* cutoff latitude for high/low statistics */
  const char *fmtstr = "%ld %8.4f %6.3f %6.3f %6.3f %10.4f %10.4f %10.4f\n";
  const size_t n = 3; /* number of components to write to disk */
  const double r = mptr->r[0];
  const double theta = mptr->theta[0];
  const double phi = mptr->phi[0];
  FILE *fp[3];
  size_t idx = *index;
  size_t j, k;
  char buf[2048];

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
          fprintf(stderr, "invert_residual_print_observatory: fp[%zu] is NULL\n", j);
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
      size_t rstat_idx;

      if (MAGDATA_Discarded(mptr->flags[j]))
        continue;

      if (!MAGDATA_FitMF(mptr->flags[j]))
        continue;

      if (mptr->qdlat[j] > qdlat_cutoff)
        rstat_idx = 0; /* north pole latitudes */
      else if (mptr->qdlat[j] < -qdlat_cutoff)
        rstat_idx = 1; /* south pole latitudes */
      else
        rstat_idx = 2; /* mid latitudes */

      /* evaluate internal field models */
      invert_eval_dBdt(mptr->t[j], r, theta, phi, dBdt_fit, w);

      res_dBdt[0] = mptr->dXdt_nec[j] - dBdt_fit[0];
      res_dBdt[1] = mptr->dYdt_nec[j] - dBdt_fit[1];
      res_dBdt[2] = mptr->dZdt_nec[j] - dBdt_fit[2];

      if (MAGDATA_ExistDXDT(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[0], fmtstr, unix_time, mptr->qdlat[j], ws, wr, wf, mptr->dXdt_nec[j], dBdt_fit[0], res_dBdt[0]);
          gsl_rstat_add(res_dBdt[0], rstat_DXDT[rstat_idx]);

          ++idx;
        }

      if (MAGDATA_ExistDYDT(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[1], fmtstr, unix_time, mptr->qdlat[j], ws, wr, wf, mptr->dYdt_nec[j], dBdt_fit[1], res_dBdt[1]);
          gsl_rstat_add(res_dBdt[1], rstat_DYDT[rstat_idx]);

          ++idx;
        }

      if (MAGDATA_ExistDZDT(mptr->flags[j]))
        {
          double ws = gsl_vector_get(w->wts_spatial, idx);
          double wr = gsl_vector_get(w->wts_robust, idx);
          double wf = gsl_vector_get(w->wts_final, idx);

          fprintf(fp[2], fmtstr, unix_time, mptr->qdlat[j], ws, wr, wf, mptr->dZdt_nec[j], dBdt_fit[2], res_dBdt[2]);
          gsl_rstat_add(res_dBdt[2], rstat_DZDT[rstat_idx]);

          ++idx;
        }
    }

  *index = idx;

  for (j = 0; j < n; ++j)
    fclose(fp[j]);

  return s;
}
#endif
