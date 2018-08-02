/*
 * mfield_main.c
 *
 * Usage:
 * ./mfield [flags]
 *
 * Flags:
 *   -c coef_output_file
 *   -n max_iterations
 *   -e epoch_decimal_year
 *   -p euler_period_days
 *   -r residual_file
 *   -l Lcurve_data_file
 *
 * After each iteration, the file 'res.#.dat' is written
 * where # is the iteration number. This file contains the
 * residuals of a sample of the DMSP dataset.
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <libconfig.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>

#include <satdata/satdata.h>

#include <common/common.h>
#include <common/oct.h>
#include <msynth/msynth.h>

#include "euler.h"
#include "magdata.h"
#include "mfield.h"
#include "mfield_error.h"
#include "mfield_residual.h"
#include "mfield_synth.h"
#include "track.h"

#define MAX_BUFFER           2048

#if 0

/*
parse_input()
  Read in index file for a given satellite, perform data
selection and downsampling
*/

satdata_mag *
parse_input(const size_t sat_idx)
{
  satdata_mag *data = NULL;
  const char *idxfile = index_files[sat_idx];
  size_t i;

  fprintf(stderr, "parse_input: reading %s...", idxfile);

  if (sat_idx >= IDX_SWA && sat_idx <= IDX_SWC)
    data = satdata_swarm_read_idx(idxfile, 0);
  else if (sat_idx == IDX_CHAMP)
    data = satdata_champ_read_idx(idxfile, 0);
  else if (sat_idx == IDX_OERSTED)
    data = satdata_oersted_read_idx(idxfile, 0);
  else if (sat_idx >= IDX_F15 && sat_idx <= IDX_F18)
    data = satdata_dmsp_read_idx(idxfile, 0);

  if (!data)
    return NULL;

  fprintf(stderr, "done (%zu points read)\n", data->n);

  if (sat_idx >= IDX_SWA && sat_idx <= IDX_SWC)
    {
      size_t nrms;
      track_workspace *track_p = track_alloc();
      double thresh[] = { 20.0, 25.0, 15.0, 15.0 };

      satdata_swarm_filter_instrument(1, data);

      /* filter by track rms */

      track_init(data, NULL, track_p);

      nrms = track_flag_rms("swarm_rms.dat", thresh, data, track_p);
      fprintf(stderr, "parse_input: flagged (%zu/%zu) (%.1f%%) points due to high rms\n",
              nrms, data->n, (double) nrms / (double) data->n * 100.0);

      track_free(track_p);

      satdata_filter_wmm(1, data);
    }

  fprintf(stderr, "parse_input: downsampling data by factor %d...", DOWNSAMPLE);

  for (i = 0; i < data->n; ++i)
    {
      if (i % DOWNSAMPLE != 0)
        data->flags[i] |= SATDATA_FLG_OUTLIER;
    }

  fprintf(stderr, "done\n");

  /* flag local time */
  if (!(sat_idx >= IDX_F15 && sat_idx <= IDX_F18))
    {
      size_t nlt;
      const double lt_min = 5.0;
      double lt_max = 22.0;
      const double euler_lt_min = 6.0;
      const double euler_lt_max = 18.0;

      /* in the first half of 2013, Oersted is in a ~10am/10pm orbit */
      if (sat_idx == IDX_OERSTED)
        lt_max = 20.0;

      fprintf(stderr, "parse_input: flagging points inside LT window [%g,%g], euler [%g,%g]...",
              lt_min, lt_max, euler_lt_min, euler_lt_max);

      nlt = flag_local_time(lt_min, lt_max, euler_lt_min, euler_lt_max, data);

      fprintf(stderr, "done (%zu/%zu data flagged)\n", nlt, data->n);
    }

  {
    size_t nflagged = satdata_nflagged(data);
    fprintf(stderr, "parse_input: total flagged points: %zu/%zu (%.1f%%) (%zu remaining)\n",
            nflagged, data->n, (double)nflagged / (double)data->n * 100.0,
            data->n - nflagged);
  }

  return data;
} /* parse_input() */

#endif /* 0 */

/*
initial_guess()
  Construct initial guess for main field coefficients. These
are based on the relevant IGRF coefficients, extrapolated forward
to the desired epoch using the SV coefficients. Initial SA coefficients
are set to 0.
*/

int
initial_guess(gsl_vector *c, mfield_workspace *w)
{
  gsl_vector_set_zero(c);

#if 0 /* crustal field modeling? */
  return 0; /*XXX*/
#endif

  {
    msynth_workspace *msynth_p = msynth_igrf_read(MSYNTH_IGRF_FILE);
    const size_t nmax = GSL_MIN(w->nmax_max, msynth_p->nmax);
    const double t = w->epoch;                       /* desired epoch */
    const double t0 = msynth_get_epoch(t, msynth_p); /* IGRF epoch */
    const double dt = t - t0;
    size_t n;
    int m;

    for (n = 1; n <= nmax; ++n)
      {
        int ni = (int) n;

        for (m = -ni; m <= ni; ++m)
          {
            size_t midx = msynth_nmidx(n, m, msynth_p);
            size_t cidx = mfield_coeff_nmidx(n, m);
            double gnm = msynth_get_mf(t, midx, msynth_p);
            double dgnm = msynth_get_sv(t, midx, msynth_p);

            /*
             * use SV prediction to update main field coefficients for new
             * epoch
             */
            mfield_set_mf(c, cidx, gnm + dt * dgnm, w);
            mfield_set_sv(c, cidx, dgnm, w);
            mfield_set_sa(c, cidx, 0.0, w);
          }
      }

    msynth_free(msynth_p);
  }

  return 0;
} /* initial_guess() */

int
print_spectrum(const char *filename, mfield_workspace *w)
{
  const double c = 3485.0;               /* Earth core radius */
  const double ratio = w->params.R / c;  /* a / c */
  size_t n;
  FILE *fp = fopen(filename, "w");

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
  fprintf(fp, "# Field %zu: MF power R_n at Earth surface\n", n++);
  fprintf(fp, "# Field %zu: SV power R_n at Earth surface\n", n++);
  fprintf(fp, "# Field %zu: SA power R_n at Earth surface\n", n++);
  fprintf(fp, "# Field %zu: MF power R_n at CMB\n", n++);
  fprintf(fp, "# Field %zu: SV power R_n at CMB\n", n++);
  fprintf(fp, "# Field %zu: SA power R_n at CMB\n", n++);

  fprintf(stderr, "print_spectrum: writing spectrum to %s...", filename);
  for (n = 1; n <= w->nmax_max; ++n)
    {
      double gn = mfield_spectrum(n, w);
      double dgn = mfield_spectrum_sv(n, w);
      double ddgn = mfield_spectrum_sa(n, w);
      double rterm = pow(ratio, 2.0*n + 4.0);

      fprintf(fp, "%zu %.12e %.12e %.12e %.12e %.12e %.12e\n",
              n,
              gn,
              dgn,
              ddgn,
              rterm * gn,
              rterm * dgn,
              rterm * ddgn);
    }
  fprintf(stderr, "done\n");

  fclose(fp);

  return 0;
} /* print_spectrum() */

int
mfield_print_residual_stat(const char *component_str, const gsl_rstat_workspace *rstat_p)
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
mfield_print_residual(const char *prefix, const size_t iter, mfield_workspace *w)
{
  int s = 0;
  char buf[2048];
  FILE *fp[12];
  const size_t n = 12; /* number of components to write to disk */
  const char *fmtstr = "%ld %.8f %.4f %.4f %.4f %.4f %.3e %.3e %.3e %.4f %.4f %.4f %.4f\n";
  const char *fmtstr_F = "%ld %.8f %.4f %.4f %.4f %.4f %.3e %.3e %.3e %.4f %.4f %.4f\n";
  const char *fmtstr_grad = "%ld %.8f %.4f %.4f %.4f %.4f %.3e %.3e %.3e %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n";
  const double qdlat_cutoff = w->params.qdlat_fit_cutoff; /* cutoff latitude for high/low statistics */
  size_t i;
  mfield_data_workspace *data_p = w->data_workspace_p;
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

  for (i = 0; i < data_p->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, data_p);
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
      fprintf(fp[0], "# X vector residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[1], "# Y vector residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[2], "# Z vector residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[3], "# F scalar residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[4], "# DX gradient (N/S) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[5], "# DY gradient (N/S) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[6], "# DZ gradient (N/S) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[7], "# DZ gradient (N/S) scalar residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[8], "# DX gradient (E/W) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[9], "# DY gradient (E/W) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[10], "# DZ gradient (E/W) vector residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);
      fprintf(fp[11], "# DF gradient (E/W) scalar residuals for MF modeling (satellite %zu, iteration %zu)\n", i, iter);

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
                  size_t euler_idx = mfield_euler_idx(i, mptr->t[j], w);
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

      fprintf(stderr, "=== FIT STATISTICS SATELLITE %zu ===\n", i);

      /* print header */
      mfield_print_residual_stat(NULL, NULL);

      mfield_print_residual_stat("X", rstat_x);
      mfield_print_residual_stat("Y", rstat_y);
      mfield_print_residual_stat("Z", rstat_z);
      mfield_print_residual_stat("F", rstat_f);
      mfield_print_residual_stat("low Z", rstat_lowz);
      mfield_print_residual_stat("high Z", rstat_highz);
      mfield_print_residual_stat("low F", rstat_lowf);
      mfield_print_residual_stat("high F", rstat_highf);

      mfield_print_residual_stat("N/S DX", rstat_dx_ns);
      mfield_print_residual_stat("N/S DY", rstat_dy_ns);
      mfield_print_residual_stat("low N/S DZ", rstat_low_dz_ns);
      mfield_print_residual_stat("high N/S DZ", rstat_high_dz_ns);
      mfield_print_residual_stat("N/S DF", rstat_df_ns);

      mfield_print_residual_stat("E/W DX", rstat_dx_ew);
      mfield_print_residual_stat("E/W DY", rstat_dy_ew);
      mfield_print_residual_stat("low E/W DZ", rstat_low_dz_ew);
      mfield_print_residual_stat("high E/W DZ", rstat_high_dz_ew);
      mfield_print_residual_stat("E/W DF", rstat_df_ew);
    }

  assert(idx == w->nres);

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

static int
check_parameters(const mfield_parameters * mfield_params,
                 const mfield_data_parameters * data_params)
{
  int s = 0;

  if (data_params->qdlat_fit_cutoff < 0.0)
    {
      fprintf(stderr, "check_parameters: qdlat_fit_cutoff must be > 0\n");
      ++s;
    }

  if (mfield_params->epoch < 0.0)
    {
      fprintf(stderr, "check_parameters: epoch must be > 0\n");
      ++s;
    }

  if (mfield_params->qdlat_fit_cutoff < 0.0)
    {
      fprintf(stderr, "check_parameters: qdlat_fit_cutoff must be > 0\n");
      ++s;
    }

  if (mfield_params->R < 0.0)
    {
      fprintf(stderr, "check_parameters: R must be > 0\n");
      ++s;
    }

  if (mfield_params->fit_euler)
    {
      if (!data_params->fit_X)
        {
          fprintf(stderr, "check_parameters: fitting Euler angles but not fitting X component\n");
          ++s;
        }

      if (!data_params->fit_Y)
        {
          fprintf(stderr, "check_parameters: fitting Euler angles but not fitting Y component\n");
          ++s;
        }

      if (!data_params->fit_Z)
        {
          fprintf(stderr, "check_parameters: fitting Euler angles but not fitting Z component\n");
          ++s;
        }
    }

  return s;
}

static int
parse_config_file(const char *filename, mfield_parameters *mfield_params,
                  mfield_data_parameters *data_params)
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

  if (config_lookup_int(&cfg, "nmax_mf", &ival))
    mfield_params->nmax_mf = (size_t) ival;
  if (config_lookup_int(&cfg, "nmax_sv", &ival))
    mfield_params->nmax_sv = (size_t) ival;
  if (config_lookup_int(&cfg, "nmax_sa", &ival))
    mfield_params->nmax_sa = (size_t) ival;

  if (config_lookup_float(&cfg, "epoch", &fval))
    mfield_params->epoch = fval;
  if (config_lookup_float(&cfg, "R", &fval))
    mfield_params->R = fval;
  if (config_lookup_float(&cfg, "euler_period", &fval))
    mfield_params->euler_period = fval;

  if (config_lookup_int(&cfg, "max_iter", &ival))
    mfield_params->max_iter = (size_t) ival;
  if (config_lookup_float(&cfg, "qdlat_fit_cutoff", &fval))
    {
      data_params->qdlat_fit_cutoff = fval;
      mfield_params->qdlat_fit_cutoff = fval;
    }
  if (config_lookup_int(&cfg, "fit_mf", &ival))
    mfield_params->fit_mf = ival;
  if (config_lookup_int(&cfg, "fit_sv", &ival))
    mfield_params->fit_sv = ival;
  if (config_lookup_int(&cfg, "fit_sa", &ival))
    mfield_params->fit_sa = ival;
  if (config_lookup_int(&cfg, "fit_euler", &ival))
    mfield_params->fit_euler = ival;
  if (config_lookup_int(&cfg, "fit_ext", &ival))
    mfield_params->fit_ext = ival;
  if (config_lookup_int(&cfg, "fit_cbias", &ival))
    mfield_params->fit_cbias = ival;

  if (config_lookup_int(&cfg, "scale_time", &ival))
    mfield_params->scale_time = ival;
  if (config_lookup_int(&cfg, "use_weights", &ival))
    mfield_params->use_weights = ival;
  if (config_lookup_int(&cfg, "regularize", &ival))
    mfield_params->regularize = ival;

  if (config_lookup_float(&cfg, "lambda_mf", &fval))
    mfield_params->lambda_mf = fval;
  if (config_lookup_float(&cfg, "lambda_sv", &fval))
    mfield_params->lambda_sv = fval;
  if (config_lookup_float(&cfg, "lambda_sa", &fval))
    mfield_params->lambda_sa = fval;

  if (config_lookup_float(&cfg, "weight_X", &fval))
    mfield_params->weight_X = fval;
  if (config_lookup_float(&cfg, "weight_Y", &fval))
    mfield_params->weight_Y = fval;
  if (config_lookup_float(&cfg, "weight_Z", &fval))
    mfield_params->weight_Z = fval;
  if (config_lookup_float(&cfg, "weight_F", &fval))
    mfield_params->weight_F = fval;
  if (config_lookup_float(&cfg, "weight_DXDT", &fval))
    mfield_params->weight_DXDT = fval;
  if (config_lookup_float(&cfg, "weight_DYDT", &fval))
    mfield_params->weight_DYDT = fval;
  if (config_lookup_float(&cfg, "weight_DZDT", &fval))
    mfield_params->weight_DZDT = fval;
  if (config_lookup_float(&cfg, "weight_DX", &fval))
    mfield_params->weight_DX = fval;
  if (config_lookup_float(&cfg, "weight_DY", &fval))
    mfield_params->weight_DY = fval;
  if (config_lookup_float(&cfg, "weight_DZ", &fval))
    mfield_params->weight_DZ = fval;

  if (config_lookup_int(&cfg, "fit_X", &ival))
    data_params->fit_X = ival;
  if (config_lookup_int(&cfg, "fit_Y", &ival))
    data_params->fit_Y = ival;
  if (config_lookup_int(&cfg, "fit_Z", &ival))
    data_params->fit_Z = ival;
  if (config_lookup_int(&cfg, "fit_F", &ival))
    data_params->fit_F = ival;

  if (config_lookup_int(&cfg, "fit_DXDT", &ival))
    data_params->fit_DXDT = ival;
  if (config_lookup_int(&cfg, "fit_DYDT", &ival))
    data_params->fit_DYDT = ival;
  if (config_lookup_int(&cfg, "fit_DZDT", &ival))
    data_params->fit_DZDT = ival;

  if (config_lookup_int(&cfg, "fit_DX_NS", &ival))
    data_params->fit_DX_NS = ival;
  if (config_lookup_int(&cfg, "fit_DY_NS", &ival))
    data_params->fit_DY_NS = ival;
  if (config_lookup_int(&cfg, "fit_DZ_NS", &ival))
    data_params->fit_DZ_NS = ival;
  if (config_lookup_int(&cfg, "fit_DF_NS", &ival))
    data_params->fit_DF_NS = ival;

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

  if (config_lookup_int(&cfg, "synth_data", &ival))
    mfield_params->synth_data = ival;
  if (config_lookup_int(&cfg, "synth_noise", &ival))
    mfield_params->synth_noise = ival;
  if (config_lookup_int(&cfg, "synth_nmin", &ival))
    mfield_params->synth_nmin = (size_t) ival;

  config_destroy(&cfg);

  return 0;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options] sat1.dat sat2.dat ...\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --maxit | -n num_iterations     - number of robust iterations\n");
  fprintf(stderr, "\t --output_file | -o file         - coefficient output file (ASCII)\n");
  fprintf(stderr, "\t --epoch | -e epoch              - model epoch in decimal years\n");
  fprintf(stderr, "\t --euler | -p period             - Euler bin size in days\n");
  fprintf(stderr, "\t --print_residuals | -r          - write residuals at each iteration\n");
  fprintf(stderr, "\t --lcurve_file | -l file         - L-curve data file\n");
  fprintf(stderr, "\t --tmin | -b min_time            - minimum data period time in decimal years\n");
  fprintf(stderr, "\t --tmax | -c max_time            - maximum data period time in decimal years\n");
  fprintf(stderr, "\t --print_data | -d               - print data used for MF modeling to output directory\n");
  fprintf(stderr, "\t --print_map | -m                - print spatial data map files to output directory\n");
  fprintf(stderr, "\t --config_file | -C file         - configuration file\n");
  fprintf(stderr, "\t --lambda_sv | -v lambda_sv      - secular variation damping parameter\n");
  fprintf(stderr, "\t --lambda_sa | -a lambda_sa      - secular acceleration damping parameter\n");
} /* print_help() */

int
main(int argc, char *argv[])
{
  int status;
  const char *error_file = "error.txt";
  char *outfile = NULL;
  char *Lfile = NULL;
  char *datamap_prefix = "output";
  char *data_prefix = "output";
  char *residual_prefix = "output";
  char *config_file = "MF.cfg"; /* default config file */
  mfield_workspace *mfield_workspace_p;
  mfield_parameters mfield_params;
  mfield_data_workspace *mfield_data_p;
  mfield_data_parameters data_params;
  gsl_vector *coeffs; /* model coefficients */
  size_t iter = 0;
  size_t maxit = 0;
  double epoch = -1.0;        /* model epoch */
  double euler_period = -1.0; /* set to 0 for single set of angles */
  double tmin = -1.0;         /* minimum time for data in years */
  double tmax = -1.0;         /* maximum time for data in years */
  int nsource = 0;            /* number of data sources (satellites/observatories) */
  int print_data = 0;         /* print data for MF modeling */
  int print_map = 0;          /* print data maps */
  int print_residuals = 0;    /* print residuals at each iteration */
  double lambda_sv = -1.0;    /* SV damping parameter */
  double lambda_sa = -1.0;    /* SA damping parameter */
  double sigma = -1.0;        /* sigma for artificial noise */
  double bias = 0.0;          /* bias for artificial noise */
  struct timeval tv0, tv1;
  char buf[MAX_BUFFER];

  data_params.qdlat_fit_cutoff = -1.0;

  mfield_init_params(&mfield_params);

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "print_residuals", no_argument, NULL, 'r' },
          { "output_file", required_argument, NULL, 'o' },
          { "epoch", required_argument, NULL, 'e' },
          { "lcurve_file", required_argument, NULL, 'l' },
          { "maxit", required_argument, NULL, 'n' },
          { "tmin", required_argument, NULL, 'b' },
          { "tmax", required_argument, NULL, 'c' },
          { "euler", required_argument, NULL, 'p' },
          { "print_data", required_argument, NULL, 'd' },
          { "print_map", required_argument, NULL, 'm' },
          { "config_file", required_argument, NULL, 'C' },
          { "lambda_sa", required_argument, NULL, 'a' },
          { "lambda_sv", required_argument, NULL, 'v' },
          { "sigma", required_argument, NULL, 'S' },
          { "bias", required_argument, NULL, 'B' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:B:c:C:de:l:mn:o:p:rv:S:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            lambda_sa = atof(optarg);
            break;

          case 'v':
            lambda_sv = atof(optarg);
            break;

          case 'b':
            tmin = atof(optarg);
            break;

          case 'c':
            tmax = atof(optarg);
            break;

          case 'C':
            config_file = optarg;
            break;

          case 'd':
            print_data = 1;
            break;

          case 'm':
            print_map = 1;
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'n':
            maxit = (size_t) atoi(optarg);
            break;

          case 'e':
            epoch = atof(optarg);
            break;

          case 'p':
            euler_period = atof(optarg);
            break;

          case 'r':
            print_residuals = 1;
            break;

          case 'l':
            Lfile = optarg;
            break;

          case 'B':
            bias = atof(optarg);
            break;

          case 'S':
            sigma = atof(optarg);
            break;

          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  nsource = argc - optind;
  if (nsource == 0)
    {
      print_help(argv);
      exit(1);
    }

  /* parse configuration file */
  fprintf(stderr, "main: parsing configuration file %s...", config_file);
  parse_config_file(config_file, &mfield_params, &data_params);
  fprintf(stderr, "done\n");

  /* check if any command-line arguments should override config file values */
  if (epoch > 0.0)
    mfield_params.epoch = epoch;
  if (euler_period >= 0.0)
    mfield_params.euler_period = euler_period;
  if (maxit > 0)
    mfield_params.max_iter = maxit;
  if (lambda_sv >= 0.0)
    mfield_params.lambda_sv = lambda_sv;
  if (lambda_sa >= 0.0)
    mfield_params.lambda_sa = lambda_sa;

  status = check_parameters(&mfield_params, &data_params);
  if (status)
    exit(1);

  /* make sure config flags are consistent with what components are fitted */
  if (!mfield_params.fit_mf)
    {
      data_params.fit_X = data_params.fit_Y = data_params.fit_Z = data_params.fit_F = 0;
      data_params.fit_Z_highlat = data_params.fit_F_highlat = 0;
      data_params.fit_DZ_NS_highlat = data_params.fit_DF_NS_highlat = 0;
      data_params.fit_DZ_EW_highlat = data_params.fit_DF_EW_highlat = 0;
      data_params.fit_DX_NS = data_params.fit_DY_NS = data_params.fit_DZ_NS = data_params.fit_DF_NS = 0;
      data_params.fit_DX_EW = data_params.fit_DY_EW = data_params.fit_DZ_EW = data_params.fit_DF_EW = 0;
    }

  if (!mfield_params.fit_sv)
    {
      data_params.fit_DXDT = data_params.fit_DYDT = data_params.fit_DZDT = 0;
    }

  fprintf(stderr, "main: epoch = %.2f\n", mfield_params.epoch);
  fprintf(stderr, "main: radius = %g [km]\n", mfield_params.R);

  if (mfield_params.fit_mf)
    {
      fprintf(stderr, "main: MF nmax = %zu\n", mfield_params.nmax_mf);
      fprintf(stderr, "main: MF damping = %g\n", mfield_params.lambda_mf);
    }

  if (mfield_params.fit_sv)
    {
      fprintf(stderr, "main: SV nmax = %zu\n", mfield_params.nmax_sv);
      fprintf(stderr, "main: SV damping = %g\n", mfield_params.lambda_sv);
    }

  if (mfield_params.fit_sv && mfield_params.fit_sa)
    {
      fprintf(stderr, "main: SA nmax = %zu\n", mfield_params.nmax_sa);
      fprintf(stderr, "main: SA damping = %g\n", mfield_params.lambda_sa);
    }

  fprintf(stderr, "main: euler period = %g [days]\n", mfield_params.euler_period);
  fprintf(stderr, "main: tmin = %g\n", tmin);
  fprintf(stderr, "main: tmax = %g\n", tmax);
  fprintf(stderr, "main: number of robust iterations = %zu\n", mfield_params.max_iter);
  fprintf(stderr, "main: number of data sources = %d\n", nsource);
  fprintf(stderr, "main: number of threads = %d\n", omp_get_max_threads());
  fprintf(stderr, "main: print_residuals = %d\n", print_residuals);
  if (outfile)
    fprintf(stderr, "main: output coefficient file = %s\n", outfile);
  if (Lfile)
    fprintf(stderr, "main: L-curve output file = %s\n", Lfile);

  /* allocate data workspace */
  data_params.epoch = mfield_params.epoch;
  mfield_data_p = mfield_data_alloc(nsource, &data_params);

  {
    int satnum = 0;

    while (optind < argc)
      {
        magdata **mdata = &(mfield_data_p->mdata[satnum]);

        assert(satnum++ < nsource);

        fprintf(stderr, "main: reading %s...", argv[optind]);
        gettimeofday(&tv0, NULL);
        *mdata = magdata_read(argv[optind], NULL);
        gettimeofday(&tv1, NULL);

        if (!(*mdata))
          exit(1);

        fprintf(stderr, "done (%zu data total, %g seconds)\n",
                (*mdata)->n, time_diff(tv0, tv1));

        magdata_init(*mdata);
        magdata_calc(*mdata);

        ++optind;
      }
  }

  {
    size_t nflag;

    /* flag any datapoints outside of [tmin,tmax] */
    fprintf(stderr, "main: flagging points outside of time [%g,%g]...", tmin, tmax);
    nflag = mfield_data_filter_time(tmin, tmax, mfield_data_p);
    fprintf(stderr, "done (%zu data flagged)\n", nflag);

    if (!mfield_params.fit_euler)
      {
        fprintf(stderr, "main: flagging Euler-only data points...");
        nflag = mfield_data_filter_euler(mfield_data_p);
        fprintf(stderr, "done (%zu data flagged)\n", nflag);
      }

    fprintf(stderr, "main: flagging non-fitted components...");
    nflag = mfield_data_filter_comp(mfield_data_p);
    fprintf(stderr, "done (%zu data flagged)\n", nflag);

    fprintf(stderr, "main: flagging sparse observatory data...");
    nflag = mfield_data_filter_observatory(mfield_data_p);
    fprintf(stderr, "done (%zu observatories flagged)\n", nflag);
  }

  if (bias > 0.0 || sigma > 0.0)
    {
      fprintf(stderr, "main: adding noise (sigma = %.1f [nT], bias = %.1f [nT])...", sigma, bias);
      mfield_data_add_noise(sigma, bias, mfield_data_p);
      fprintf(stderr, "done\n");
    }

  fprintf(stderr, "main: data epoch = %.2f\n", mfield_data_epoch(mfield_data_p));
  fprintf(stderr, "main: data tmin  = %.2f\n", satdata_epoch2year(mfield_data_p->t0_data));
  fprintf(stderr, "main: data tmax  = %.2f\n", satdata_epoch2year(mfield_data_p->t1_data));

  if (print_map)
    {
      /* print spatial coverage maps for each satellite */
      mfield_data_map(datamap_prefix, mfield_data_p);
    }

  mfield_params.nsat = nsource;
  mfield_params.mfield_data_p = mfield_data_p;

  /* allocate mfield workspace */
  mfield_workspace_p = mfield_alloc(&mfield_params);
  if (mfield_workspace_p == NULL)
    {
      fprintf(stderr, "main: mfield_alloc failed\n");
      exit(1);
    }

  fprintf(stderr, "main: number of main field parameters:           %zu\n", mfield_workspace_p->nnm_mf);
  fprintf(stderr, "main: number of secular variation parameters:    %zu\n", mfield_workspace_p->nnm_sv);
  fprintf(stderr, "main: number of secular acceleration parameters: %zu\n", mfield_workspace_p->nnm_sa);
  fprintf(stderr, "main: number of Euler angle parameters:          %zu\n", mfield_workspace_p->neuler);
  fprintf(stderr, "main: number of external field parameters:       %zu\n", mfield_workspace_p->next);
  fprintf(stderr, "main: number of crustal bias parameters:         %zu\n", mfield_workspace_p->nbias);

  if (mfield_params.synth_data)
    {
      fprintf(stderr, "main: replacing with synthetic data...");
      gettimeofday(&tv0, NULL);
      mfield_synth_replace(mfield_workspace_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
    }

  /* initialize model parameters */
  mfield_init(mfield_workspace_p);

  /* print out dataset if requested - do this after mfield_init() so
   * spatial weights are computed */
  if (print_data)
    {
      /* print data used for MF modeling for each satellite */
      fprintf(stderr, "main: printing data for MF modeling to %s...", data_prefix);
      mfield_data_print(data_prefix, mfield_workspace_p->wts_spatial, mfield_data_p);
      fprintf(stderr, "done\n");
    }

  /* construct initial guess vector from IGRF */
  coeffs = gsl_vector_alloc(mfield_workspace_p->p);
  fprintf(stderr, "main: constructing initial coefficient vector...");
  initial_guess(coeffs, mfield_workspace_p);
  fprintf(stderr, "done\n");

  gettimeofday(&tv0, NULL);

  while (iter++ < mfield_params.max_iter)
    {
      fprintf(stderr, "main: ROBUST ITERATION %zu/%zu\n", iter, mfield_params.max_iter);

      mfield_calc_nonlinear(coeffs, mfield_workspace_p);

      /* output coefficients for this iteration */
      sprintf(buf, "coef.txt.iter%zu", iter);
      fprintf(stderr, "main: writing coefficient file %s...", buf);
      mfield_write_ascii(buf, mfield_workspace_p->epoch, 0, mfield_workspace_p);
      fprintf(stderr, "done\n");

      /* output spectrum for this iteration */
      sprintf(buf, "mfield.s.iter%zu", iter);
      print_spectrum(buf, mfield_workspace_p);

      if (print_residuals)
        {
          fprintf(stderr, "main: printing residuals to %s...", residual_prefix);
#if 0
          mfield_print_residual(residual_prefix, iter, mfield_workspace_p);
#else
          mfield_residual_print(residual_prefix, iter, mfield_workspace_p);
#endif
          fprintf(stderr, "done\n");
        }

      /* reset workspace for a new iteration */
      mfield_reset(mfield_workspace_p);
    }

  gettimeofday(&tv1, NULL);

  fprintf(stderr, "main: total time for inversion: %.2f seconds\n", time_diff(tv0, tv1));

  /* calculate covariance matrix */
  fprintf(stderr, "main: calculating covariance matrix...");
  gettimeofday(&tv0, NULL);
  status = mfield_covariance(mfield_workspace_p->covar, mfield_workspace_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (status = %d, %g seconds)\n", status, time_diff(tv0, tv1));

  /* calculate errors in coefficients */
  fprintf(stderr, "main: printing coefficient uncertainties to %s...", error_file);
  gettimeofday(&tv0, NULL);
  mfield_print_uncertainties(error_file, mfield_workspace_p->covar, mfield_workspace_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  sprintf(buf, "mfield_coeffs.txt");
  fprintf(stderr, "main: writing coefficients to %s...", buf);
  mfield_write_ascii(buf, mfield_workspace_p->epoch, 1, mfield_workspace_p);
  fprintf(stderr, "done\n");

  {
    gsl_vector *evals = gsl_vector_alloc(mfield_workspace_p->p);
    FILE *fp;

    fprintf(stderr, "main: calculating eigenvalues of J^T J...");
    gettimeofday(&tv0, NULL);
    mfield_calc_evals(evals, mfield_workspace_p);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

    sprintf(buf, "mfield_evals.txt");
    fprintf(stderr, "main: writing eigenvalues to %s...", buf);
    fp = fopen(buf, "w");
    gsl_vector_fprintf(fp, evals, "%.12e");
    fclose(fp);
    fprintf(stderr, "done\n");

    gsl_vector_free(evals);
  }

  /* L-curve data */
  if (Lfile)
    {
      gsl_vector *f = gsl_multifit_nlinear_residual(mfield_workspace_p->multifit_nlinear_p);
      double xnorm = gsl_blas_dnrm2(mfield_workspace_p->c);
      double fnorm = gsl_blas_dnrm2(f);
      FILE *fp = fopen(Lfile, "a");

      if (!fp)
        {
          fprintf(stderr, "main: unable to open %s: %s\n",
                  Lfile, strerror(errno));
        }
      else
        {
          fprintf(fp, "%.12e %.12e %f %f\n",
                  log(fnorm),
                  log(xnorm),
                  mfield_params.lambda_sv,
                  mfield_params.lambda_sa);

          fclose(fp);
        }
    }

  if (outfile)
    {
      fprintf(stderr, "main: writing ASCII coefficients to %s...", outfile);
      /*mfield_write(outfile, mfield_workspace_p);*/
      mfield_write_ascii(outfile, mfield_workspace_p->epoch, 0, mfield_workspace_p);
      fprintf(stderr, "done\n");
    }

  print_spectrum("mfield.s", mfield_workspace_p);

  /* print coefficients */
  {
    size_t n;
    int m;

#if MFIELD_FIT_EXTFIELD
    char *ext_file = "coeffs.ext";
    FILE *fp = fopen(ext_file, "w");

    fprintf(stderr, "main: printing external coefficients to %s...", ext_file);
    
    for (n = 0; n < mfield_workspace_p->next; ++n)
      {
        size_t idx = mfield_workspace_p->ext_offset + n;
        double k = gsl_vector_get(coeffs, idx);

        fprintf(fp, "%d %g\n", mfield_workspace_p->ext_fdayi[n], k);
      }

    fprintf(stderr, "done\n");

    fclose(fp);
#endif

    /* print Euler angles */
    if (mfield_params.fit_euler)
      {
        for (n = 0; n < mfield_workspace_p->nsat; ++n)
          {
            magdata *mptr = mfield_data_ptr(n, mfield_workspace_p->data_workspace_p);

            if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
              {
                double t0 = mfield_workspace_p->data_workspace_p->t0[n];
                size_t euler_idx = mfield_euler_idx(n, t0, mfield_workspace_p);
                double alpha = gsl_vector_get(coeffs, euler_idx);
                double beta = gsl_vector_get(coeffs, euler_idx + 1);
                double gamma = gsl_vector_get(coeffs, euler_idx + 2);
                char filename[2048];

                fprintf(stderr, "main: satellite %zu: alpha = %f beta = %f gamma = %f [deg]\n",
                        n,
                        wrap180(alpha * 180.0 / M_PI),
                        wrap180(beta * 180.0 / M_PI),
                        wrap180(gamma * 180.0 / M_PI));

                sprintf(filename, "euler.%zu", n);
                fprintf(stderr, "main: satellite %zu: printing Euler angles to %s...", n, filename);
                mfield_euler_print(filename, n, mfield_workspace_p);
                fprintf(stderr, "done\n");
              }
          }
      }

    fprintf(stderr, "main: printing internal coefficients up to degree 3\n");
    for (n = 1; n <= GSL_MIN(3, mfield_workspace_p->nmax_max); ++n)
      {
        int ni = (int) n;
        for (m = -ni; m <= ni; ++m)
          {
            int mabs = abs(m);
            size_t cidx = mfield_coeff_nmidx(n, m);
            char c = (m < 0) ? 'h' : 'g';
            double gnm = mfield_get_mf(coeffs, cidx, mfield_workspace_p);
            double dgnm = mfield_get_sv(coeffs, cidx, mfield_workspace_p);
            double ddgnm = mfield_get_sa(coeffs, cidx, mfield_workspace_p);

            fprintf(stderr, "%c(%d,%d) = %12g (%12g,%12g)\n", c, ni, mabs,
                    gnm, dgnm, ddgnm);
          }
      }
  }

  /* print residual norm between synthetic and computed coefficients */
  if (mfield_params.synth_data)
    {
      gsl_vector *g_synth = gsl_vector_alloc(mfield_workspace_p->p_int);
      gsl_vector_view g = gsl_vector_subvector(coeffs, 0, mfield_workspace_p->p_int);
      double norm;

      /* synthetic internal coefficients */
      mfield_synth_g(g_synth, mfield_workspace_p);

      /* subtract model internal coefficients */
      gsl_vector_sub(g_synth, &g.vector);

      norm = gsl_blas_dnrm2(g_synth);

      fprintf(stderr, "main: || g_synth - g || = %.12e [nT]\n", norm);

      gsl_vector_free(g_synth);
    }

  mfield_free(mfield_workspace_p);
  mfield_data_free(mfield_data_p);
  gsl_vector_free(coeffs);

  return 0;
} /* main() */
