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
#include <track/track.h>

#include "euler.h"
#include "magdata.h"
#include "mfield.h"
#include "mfield_euler.h"
#include "mfield_fluxcal.h"
#include "mfield_error.h"
#include "mfield_residual.h"
#include "mfield_synth.h"

/*XXX*/
#include "task2.c"

#define MAX_BUFFER           2048

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
  size_t i;
  msynth_workspace *msynth_p = msynth_igrf_read(MSYNTH_IGRF_FILE);
  /*msynth_workspace *msynth_p = msynth_shc_read(MSYNTH_CHAOS_FILE);*/

  gsl_vector_set_zero(c);

  mfield_fill_g(c, msynth_p, NULL, w);

  mfield_write_ascii("initial.txt", 2015.5, c, w);
  mfield_write_shc("initial.shc", c, w);

  /* initialize fluxgate calibration scale factors to 1.0, leaving offsets and angles at 0 */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      int fit_fluxcal = w->params.fit_fluxcal && (mptr->global_flags & MAGDATA_GLOBFLG_FLUXCAL);

      if (fit_fluxcal)
        {
          size_t ncontrol = gsl_bspline2_ncontrol(w->fluxcal_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)]);
          gsl_vector_view tmp = gsl_vector_subvector(c, w->fluxcal_offset + w->offset_fluxcal[i], FLUXCAL_P * ncontrol);
          gsl_matrix_view control_pts = gsl_matrix_view_vector(&tmp.vector, FLUXCAL_P, ncontrol);
          gsl_matrix_view m = gsl_matrix_submatrix(&control_pts.matrix, FLUXCAL_IDX_SX, 0, 3, ncontrol);
          gsl_matrix_set_all(&m.matrix, 1.0);
        }
    }

  msynth_free(msynth_p);

  return 0;
}

int
print_spectrum(const char *filename, mfield_workspace *w)
{
  const double c = 3485.0;               /* Earth core radius */
  const double ratio = w->params.R / c;  /* a / c */
  const double epoch = w->params.epoch;
  size_t n;
  FILE *fp = fopen(filename, "w");

  n = 1;
  fprintf(fp, "# Spectrum epoch: %g\n", epoch);
  fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
  fprintf(fp, "# Field %zu: MF power R_n at Earth surface\n", n++);
  fprintf(fp, "# Field %zu: SV power R_n at Earth surface\n", n++);
  fprintf(fp, "# Field %zu: SA power R_n at Earth surface\n", n++);
  fprintf(fp, "# Field %zu: MF power R_n at CMB\n", n++);
  fprintf(fp, "# Field %zu: SV power R_n at CMB\n", n++);
  fprintf(fp, "# Field %zu: SA power R_n at CMB\n", n++);

  fprintf(stderr, "print_spectrum: writing spectrum to %s...", filename);
  for (n = 1; n <= w->nmax; ++n)
    {
      double gn = mfield_spectrum(epoch, n, 0, w);
      double dgn = mfield_spectrum(epoch, n, 1, w);
      double ddgn = mfield_spectrum(epoch, n, 2, w);
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
}

int
print_gnm(const char *prefix, const size_t iter, mfield_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax_core;
  const double t0 = epoch2year(w->data_workspace_p->t0_data);
  const double t1 = epoch2year(w->data_workspace_p->t1_data);
  const double dt = 0.05; /* step size in years */
  gsl_bspline2_workspace *gauss_spline_p = w->gauss_spline_workspace_p[0];
  const size_t ncontrol = gsl_bspline2_ncontrol(gauss_spline_p);
  const size_t order = gsl_bspline2_order(gauss_spline_p);
  const size_t nderiv = GSL_MIN(3, order - 1);
  char filename[2048];
  size_t n;

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      int m;

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          gsl_vector_const_view v = gsl_vector_const_subvector_with_stride(w->c, cidx, w->nnm_core, ncontrol);
          char gh = (m < 0) ? 'h' : 'g';
          double t;
          FILE *fp;
          size_t i;

          sprintf(filename, "%s/%c%zu_%d_iter%zu.txt", prefix, gh, n, abs(m), iter);
          fp = fopen(filename, "w");
          if (!fp)
            {
              fprintf(stderr, "print_gnm: unable to open %s: %s\n",
                      filename, strerror(errno));
              ++s;
              continue;
            }

          i = 1;
          fprintf(fp, "# %c(%zu,%d) time series\n", gh, n, abs(m));
          fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
          fprintf(fp, "# Field %zu: %c(%zu,%d)(t) (nT)\n", i++, gh, n, abs(m));

          if (nderiv > 0)
            fprintf(fp, "# Field %zu: d/dt %c(%zu,%d)(t) (nT/year)\n", i++, gh, n, abs(m));
          if (nderiv > 1)
            fprintf(fp, "# Field %zu: d^2/dt^2 %c(%zu,%d)(t) (nT/year^2)\n", i++, gh, n, abs(m));
          if (nderiv > 2)
            fprintf(fp, "# Field %zu: d^3/dt^3 %c(%zu,%d)(t) (nT/year^3)\n", i++, gh, n, abs(m));

          for (t = t0; t <= t1; t += dt)
            {
              double gnm;

              fprintf(fp, "%10.4f ", t);

              for (i = 0; i <= nderiv; ++i)
                {
                  gsl_bspline2_eval_deriv(t, &v.vector, i, &gnm, gauss_spline_p);
                  fprintf(fp, "%20.12e ", gnm);
                }

              fprintf(fp, "\n");
            }

          fclose(fp);
        }
    }

  return s;
}

int
print_euler(const char *prefix, const size_t iter, mfield_workspace *w)
{
  int s = 0;
  char filename[2048];
  size_t n;

  if (!w->params.fit_euler)
    return 0;

  for (n = 0; n < w->nsat; ++n)
    {
      magdata *mptr = mfield_data_ptr(n, w->data_workspace_p);

      if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
        {
          sprintf(filename, "%s/euler.%zu.iter%zu.txt", prefix, n, iter);
          mfield_euler_print(filename, n, w);
        }
    }

  return s;
}

int
print_fluxcal(const char *prefix, const size_t iter, mfield_workspace *w)
{
  int s = 0;
  char filename[2048];
  size_t n;

  if (!w->params.fit_fluxcal)
    return 0;

  for (n = 0; n < w->nsat; ++n)
    {
      magdata *mptr = mfield_data_ptr(n, w->data_workspace_p);

      if (mptr->global_flags & MAGDATA_GLOBFLG_FLUXCAL)
        {
          sprintf(filename, "%s/fluxcal.%zu.iter%zu.txt", prefix, n, iter);
          mfield_fluxcal_print(filename, n, w);
        }
    }

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

  if (config_lookup_int(&cfg, "nmax_core", &ival))
    mfield_params->nmax_core = (size_t) ival;
  if (config_lookup_int(&cfg, "nmax", &ival))
    mfield_params->nmax = (size_t) ival;
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
  if (config_lookup_float(&cfg, "gauss_period", &fval))
    mfield_params->gauss_period = fval;
  if (config_lookup_float(&cfg, "euler_period", &fval))
    mfield_params->euler_period = fval;
  if (config_lookup_float(&cfg, "fluxcal_period", &fval))
    mfield_params->fluxcal_period = fval;

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
  if (config_lookup_int(&cfg, "fit_fluxcal", &ival))
    mfield_params->fit_fluxcal = ival;
  if (config_lookup_int(&cfg, "fit_cbias", &ival))
    mfield_params->fit_cbias = ival;

  if (config_lookup_int(&cfg, "gauss_spline_order", &ival))
    mfield_params->gauss_spline_order = ival;
  if (config_lookup_int(&cfg, "euler_spline_order", &ival))
    mfield_params->euler_spline_order = ival;
  if (config_lookup_int(&cfg, "fluxcal_spline_order", &ival))
    mfield_params->fluxcal_spline_order = ival;

  if (config_lookup_int(&cfg, "scale_time", &ival))
    mfield_params->scale_time = ival;
  if (config_lookup_int(&cfg, "use_weights", &ival))
    mfield_params->use_weights = ival;
  if (config_lookup_int(&cfg, "regularize", &ival))
    mfield_params->regularize = ival;

  if (config_lookup_float(&cfg, "lambda_0", &fval))
    mfield_params->lambda_0 = fval;
  if (config_lookup_float(&cfg, "lambda_1", &fval))
    mfield_params->lambda_1 = fval;
  if (config_lookup_float(&cfg, "lambda_2", &fval))
    mfield_params->lambda_2 = fval;
  if (config_lookup_float(&cfg, "lambda_3", &fval))
    mfield_params->lambda_3 = fval;
  if (config_lookup_float(&cfg, "lambda_s", &fval))
    mfield_params->lambda_s = fval;
  if (config_lookup_float(&cfg, "lambda_o", &fval))
    mfield_params->lambda_o = fval;
  if (config_lookup_float(&cfg, "lambda_u", &fval))
    mfield_params->lambda_u = fval;

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

  if (config_lookup_int(&cfg, "fit_seplat", &ival))
    data_params->fit_seplat = ival;

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
  fprintf(stderr, "\t --lambda_0 | -J lambda_0        - main field damping parameter\n");
  fprintf(stderr, "\t --lambda_1 | -K lambda_1        - 1st time derivative of main field damping parameter\n");
  fprintf(stderr, "\t --lambda_2 | -L lambda_2        - 2nd time derivative of main field damping parameter\n");
  fprintf(stderr, "\t --lambda_3 | -M lambda_3        - 3rd time derivative of main field damping parameter\n");
  fprintf(stderr, "\t --lambda_s | -N lambda_s        - scale factor damping parameter\n");
  fprintf(stderr, "\t --lambda_o | -O lambda_o        - offset factor damping parameter\n");
  fprintf(stderr, "\t --lambda_u | -P lambda_u        - non-orthogonality factor damping parameter\n");
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
  double lambda_0 = -1.0;     /* MF damping parameter */
  double lambda_1 = -1.0;     /* 1st time derivative of MF damping parameter */
  double lambda_2 = -1.0;     /* 2nd time derivative of MF damping parameter */
  double lambda_3 = -1.0;     /* 3rd time derivative of MF damping parameter */
  double lambda_s = -1.0;     /* scale factor damping parameter */
  double lambda_o = -1.0;     /* offset factor damping parameter */
  double lambda_u = -1.0;     /* non-orthogonality factor damping parameter */
  double sigma = -1.0;        /* sigma for artificial noise */
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
          { "lambda_0", required_argument, NULL, 'J' },
          { "lambda_1", required_argument, NULL, 'K' },
          { "lambda_2", required_argument, NULL, 'L' },
          { "lambda_3", required_argument, NULL, 'M' },
          { "lambda_s", required_argument, NULL, 'N' },
          { "lambda_o", required_argument, NULL, 'O' },
          { "lambda_u", required_argument, NULL, 'P' },
          { "sigma", required_argument, NULL, 'S' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "b:c:C:de:l:mJ:K:L:M:N:O:P:n:o:p:rS:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'J':
            lambda_0 = atof(optarg);
            break;

          case 'K':
            lambda_1 = atof(optarg);
            break;

          case 'L':
            lambda_2 = atof(optarg);
            break;

          case 'M':
            lambda_3 = atof(optarg);
            break;

          case 'N':
            lambda_s = atof(optarg);
            break;

          case 'O':
            lambda_o = atof(optarg);
            break;

          case 'P':
            lambda_u = atof(optarg);
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

  mfield_params.nmax = GSL_MAX(mfield_params.nmax_core, mfield_params.nmax);

  /* check if any command-line arguments should override config file values */
  if (epoch > 0.0)
    mfield_params.epoch = epoch;
  if (euler_period >= 0.0)
    mfield_params.euler_period = euler_period;
  if (maxit > 0)
    mfield_params.max_iter = maxit;
  if (lambda_0 >= 0.0)
    mfield_params.lambda_0 = lambda_0;
  if (lambda_1 >= 0.0)
    mfield_params.lambda_1 = lambda_1;
  if (lambda_2 >= 0.0)
    mfield_params.lambda_2 = lambda_2;
  if (lambda_3 >= 0.0)
    mfield_params.lambda_3 = lambda_3;
  if (lambda_s >= 0.0)
    mfield_params.lambda_s = lambda_s;
  if (lambda_o >= 0.0)
    mfield_params.lambda_o = lambda_o;
  if (lambda_u >= 0.0)
    mfield_params.lambda_u = lambda_u;

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

  fprintf(stderr, "main: epoch           = %.2f\n", mfield_params.epoch);
  fprintf(stderr, "main: radius          = %g [km]\n", mfield_params.R);
  fprintf(stderr, "main: core field nmax = %zu\n", mfield_params.nmax_core);
  fprintf(stderr, "main: total nmax      = %zu\n", mfield_params.nmax);

  if (mfield_params.fit_mf)
    {
      fprintf(stderr, "main: MF nmax = %zu\n", mfield_params.nmax_mf);
      fprintf(stderr, "main: MF damping                       = %g\n", mfield_params.lambda_0);
      fprintf(stderr, "main: MF damping (1st time derivative) = %g\n", mfield_params.lambda_1);
      fprintf(stderr, "main: MF damping (2nd time derivative) = %g\n", mfield_params.lambda_2);
      fprintf(stderr, "main: MF damping (3rd time derivative) = %g\n", mfield_params.lambda_3);
    }

  if (mfield_params.fit_sv)
    fprintf(stderr, "main: SV nmax = %zu\n", mfield_params.nmax_sv);

  if (mfield_params.fit_sv && mfield_params.fit_sa)
    fprintf(stderr, "main: SA nmax = %zu\n", mfield_params.nmax_sa);

  if (mfield_params.fit_fluxcal)
    {
      fprintf(stderr, "main: scale factor damping      = %g\n", mfield_params.lambda_s);
      fprintf(stderr, "main: offset damping            = %g\n", mfield_params.lambda_o);
      fprintf(stderr, "main: non-orthogonality damping = %g\n", mfield_params.lambda_u);
    }

  fprintf(stderr, "main: Gauss period   = %g [days]\n", mfield_params.gauss_period);
  fprintf(stderr, "main: euler period   = %g [days]\n", mfield_params.euler_period);
  fprintf(stderr, "main: fluxcal period = %g [days]\n", mfield_params.fluxcal_period);
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

#if 0
  /*XXX*/
  {
    fprintf(stderr, "main: TASK2 adding noise (sigma = %.1f [nT])...", sigma);
    task2_add_noise(sigma, mfield_data_p);
    fprintf(stderr, "done\n");
  }
#endif

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

  fprintf(stderr, "main: number of core field spline parameters:    %zu\n", mfield_workspace_p->p_core);
  fprintf(stderr, "main: number of crustal field parameters:        %zu\n", mfield_workspace_p->p_crust);
  fprintf(stderr, "main: number of Euler angle parameters:          %zu\n", mfield_workspace_p->p_euler);
  fprintf(stderr, "main: number of fluxgate calibration parameters: %zu\n", mfield_workspace_p->p_fluxcal);
  fprintf(stderr, "main: number of external field parameters:       %zu\n", mfield_workspace_p->p_ext);
  fprintf(stderr, "main: number of crustal bias parameters:         %zu\n", mfield_workspace_p->p_bias);
  fprintf(stderr, "main: number of total parameters:                %zu\n", mfield_workspace_p->p);

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

  if (print_residuals)
    {
      /* print residuals with initial guess vector */
      fprintf(stderr, "main: printing initial residuals to %s...", residual_prefix);
      gsl_vector_memcpy(mfield_workspace_p->c, coeffs);
      mfield_residual_print(residual_prefix, 0, mfield_workspace_p);
      fprintf(stderr, "done\n");
    }

  gettimeofday(&tv0, NULL);

  status = GSL_CONTINUE;
  while (iter++ < mfield_params.max_iter)
    {
      fprintf(stderr, "main: ROBUST ITERATION %zu/%zu\n", iter, mfield_params.max_iter);

      status = mfield_calc_nonlinear(coeffs, mfield_workspace_p);

      /* output coefficients for this iteration */
      sprintf(buf, "coef.txt.iter%zu", iter);
      fprintf(stderr, "main: writing coefficient file %s...", buf);
      mfield_write_ascii(buf, mfield_workspace_p->epoch, coeffs, mfield_workspace_p);
      fprintf(stderr, "done\n");

      /* output coefficients for this iteration */
      sprintf(buf, "coef.shc.iter%zu", iter);
      fprintf(stderr, "main: writing coefficient file %s...", buf);
      mfield_write_shc(buf, coeffs, mfield_workspace_p);
      fprintf(stderr, "done\n");

      /* output spectrum for this iteration */
      sprintf(buf, "mfield.s.iter%zu", iter);
      print_spectrum(buf, mfield_workspace_p);

      /* output Gauss coefficient time series for this iteration */
      print_gnm("output/gnm", iter, mfield_workspace_p);

      /* output Euler angle time series for this iteration */
      print_euler("output", iter, mfield_workspace_p);

      /* output fluxgate calibration time series for this iteration */
      print_fluxcal("output", iter, mfield_workspace_p);

      if (print_residuals)
        {
          fprintf(stderr, "main: printing residuals to %s...", residual_prefix);
          mfield_residual_print(residual_prefix, iter, mfield_workspace_p);
          fprintf(stderr, "done\n");
        }

      /* reset workspace for a new iteration */
      mfield_reset(mfield_workspace_p);
    }

  gettimeofday(&tv1, NULL);

  fprintf(stderr, "main: total time for inversion: %.2f seconds\n", time_diff(tv0, tv1));

#if 1
  /* calculate covariance matrix */
  if (mfield_workspace_p->old_fdf == 0)
    {
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

      /* calculate correlation matrix */
      fprintf(stderr, "main: calculating correlation matrix...");
      gettimeofday(&tv0, NULL);
      status = mfield_correlation(mfield_workspace_p->covar, mfield_workspace_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      printsym_octave(mfield_workspace_p->covar, "corr");
    }
#endif

  sprintf(buf, "mfield_coeffs.txt");
  fprintf(stderr, "main: writing coefficients to %s...", buf);
  mfield_write_ascii(buf, mfield_workspace_p->epoch, coeffs, mfield_workspace_p);
  fprintf(stderr, "done\n");

  sprintf(buf, "mfield_coeffs.shc");
  fprintf(stderr, "main: writing SHC coefficients to %s...", buf);
  mfield_write_shc(buf, coeffs, mfield_workspace_p);
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
          fprintf(fp, "%.12e %.12e\n",
                  log(fnorm),
                  log(xnorm));

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
    
    for (n = 0; n < mfield_workspace_p->p_ext; ++n)
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
                double t0 = mptr->t[0];
                double t1 = mptr->t[mptr->n - 1];
                gsl_bspline2_workspace *euler_spline_p = mfield_workspace_p->euler_spline_workspace_p[CIDX2(n, mfield_workspace_p->nsat, 0, mfield_workspace_p->max_threads)];
                size_t ncontrol = gsl_bspline2_ncontrol(euler_spline_p);
                size_t euler_idx = mfield_workspace_p->euler_offset + mfield_workspace_p->offset_euler[n];
                double euler_data[EULER_P];
                gsl_vector_view euler_params = gsl_vector_view_array(euler_data, EULER_P);
                char filename[2048];

                /* Euler angle control points for this dataset */
                gsl_vector_const_view tmp = gsl_vector_const_subvector(coeffs, euler_idx, EULER_P * ncontrol);
                gsl_matrix_const_view control_pts = gsl_matrix_const_view_vector(&tmp.vector, EULER_P, ncontrol);

                gsl_bspline2_vector_eval(0.5*(t0+t1), &control_pts.matrix, &euler_params.vector, euler_spline_p);

                fprintf(stderr, "main: satellite %zu: alpha = %f beta = %f gamma = %f [deg]\n",
                        n,
                        wrap180(euler_data[0] * 180.0 / M_PI),
                        wrap180(euler_data[1] * 180.0 / M_PI),
                        wrap180(euler_data[2] * 180.0 / M_PI));

                sprintf(filename, "euler.%zu", n);
                fprintf(stderr, "main: satellite %zu: printing Euler angles to %s...", n, filename);
                mfield_euler_print(filename, n, mfield_workspace_p);
                fprintf(stderr, "done\n");
              }
          }
      }

    /* print fluxgate calibration parameters */
    if (mfield_params.fit_fluxcal)
      {
        for (n = 0; n < mfield_workspace_p->nsat; ++n)
          {
            magdata *mptr = mfield_data_ptr(n, mfield_workspace_p->data_workspace_p);

            if (mptr->global_flags & MAGDATA_GLOBFLG_FLUXCAL)
              {
                double t0 = mfield_workspace_p->data_workspace_p->t0[n];
                double t1 = mfield_workspace_p->data_workspace_p->t1[n];
                gsl_bspline2_workspace *fluxcal_spline_p = mfield_workspace_p->fluxcal_spline_workspace_p[CIDX2(n, mfield_workspace_p->nsat, 0, mfield_workspace_p->max_threads)];
                size_t ncontrol = gsl_bspline2_ncontrol(fluxcal_spline_p);
                size_t fluxcal_idx = mfield_workspace_p->fluxcal_offset + mfield_workspace_p->offset_fluxcal[n];
                double cal_data[FLUXCAL_P];
                gsl_vector_view cal_params = gsl_vector_view_array(cal_data, FLUXCAL_P);
                char filename[2048];

                /* fluxgate calibration control points for this dataset */
                gsl_vector_const_view tmp = gsl_vector_const_subvector(coeffs, fluxcal_idx, FLUXCAL_P * ncontrol);
                gsl_matrix_const_view control_pts = gsl_matrix_const_view_vector(&tmp.vector, FLUXCAL_P, ncontrol);

                gsl_bspline2_vector_eval(0.5*(t0+t1), &control_pts.matrix, &cal_params.vector, fluxcal_spline_p);

                fprintf(stderr, "main: satellite %zu: S = %e %e %e\n",
                        n,
                        cal_data[FLUXCAL_IDX_SX],
                        cal_data[FLUXCAL_IDX_SY],
                        cal_data[FLUXCAL_IDX_SZ]);
                fprintf(stderr, "                   O = %e %e %e [nT]\n",
                        cal_data[FLUXCAL_IDX_OX],
                        cal_data[FLUXCAL_IDX_OY],
                        cal_data[FLUXCAL_IDX_OZ]);
                fprintf(stderr, "                   U = %e %e %e [deg]\n",
                        cal_data[FLUXCAL_IDX_U1] * 180.0 / M_PI,
                        cal_data[FLUXCAL_IDX_U2] * 180.0 / M_PI,
                        cal_data[FLUXCAL_IDX_U3] * 180.0 / M_PI);

                sprintf(filename, "fluxcal.%zu", n);
                fprintf(stderr, "main: satellite %zu: printing fluxgate calibration splines to %s...", n, filename);
                mfield_fluxcal_print(filename, n, mfield_workspace_p);
                fprintf(stderr, "done\n");
              }
          }
      }

    {
      gsl_bspline2_workspace *gauss_spline_p = mfield_workspace_p->gauss_spline_workspace_p[0];
      size_t ncontrol = gsl_bspline2_ncontrol(gauss_spline_p);
      double t = mfield_params.epoch;

      fprintf(stderr, "main: printing internal coefficients up to degree 3\n");
      for (n = 1; n <= GSL_MIN(3, mfield_workspace_p->nmax); ++n)
        {
          int ni = (int) n;
          for (m = -ni; m <= ni; ++m)
            {
              int mabs = abs(m);
              size_t cidx = mfield_coeff_nmidx(n, m);
              gsl_vector_view v = gsl_vector_subvector_with_stride(coeffs, cidx, mfield_workspace_p->nnm_core, ncontrol);
              char c = (m < 0) ? 'h' : 'g';
              double gnm, dgnm, ddgnm;

              gsl_bspline2_eval(t, &v.vector, &gnm, gauss_spline_p);
              gsl_bspline2_eval_deriv(t, &v.vector, 1, &dgnm, gauss_spline_p);
              gsl_bspline2_eval_deriv(t, &v.vector, 2, &ddgnm, gauss_spline_p);

              fprintf(stderr, "%c(%d,%d) = %12g (%12g,%12g)\n", c, ni, mabs, gnm, dgnm, ddgnm);
            }
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
