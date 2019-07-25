/*
 * invert_main.c
 *
 * Usage:
 * ./invert [flags]
 *
 * Flags:
 *   -c coef_output_file
 *   -n max_iterations
 *   -e epoch_decimal_year
 *   -r residual_file
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

#include "magdata.h"
#include "invert.h"

#include "invert_multifit.c"

#define MAX_BUFFER           2048

/*
initial_guess()
  Construct initial guess for main field coefficients. These
are based on the relevant IGRF coefficients, extrapolated forward
to the desired epoch using the SV coefficients. Initial SA coefficients
are set to 0.
*/

int
initial_guess(gsl_vector *c, invert_workspace *w)
{
  gsl_vector_set_zero(c);

  return 0;
}

static int
check_parameters(const invert_parameters * invert_params,
                 const invert_data_parameters * data_params)
{
  int s = 0;

  if (data_params->qdlat_fit_cutoff < 0.0)
    {
      fprintf(stderr, "check_parameters: qdlat_fit_cutoff must be > 0\n");
      ++s;
    }

  if (invert_params->epoch < 0.0)
    {
      fprintf(stderr, "check_parameters: epoch must be > 0\n");
      ++s;
    }

  if (invert_params->qdlat_fit_cutoff < 0.0)
    {
      fprintf(stderr, "check_parameters: qdlat_fit_cutoff must be > 0\n");
      ++s;
    }

  if (invert_params->R < 0.0)
    {
      fprintf(stderr, "check_parameters: R must be > 0\n");
      ++s;
    }

  return s;
}

static int
parse_config_file(const char *filename, invert_parameters *invert_params,
                  invert_data_parameters *data_params)
{
  int s;
  config_t cfg;
  double fval;
  int ival;
  const char * sval;

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

  if (config_lookup_float(&cfg, "epoch", &fval))
    invert_params->epoch = fval;
  if (config_lookup_float(&cfg, "R", &fval))
    invert_params->R = fval;
  if (config_lookup_int(&cfg, "num_spatial_modes", &ival))
    invert_params->nspatmodes = (size_t) ival;

  if (config_lookup_string(&cfg, "tmode_file", &sval))
    strcpy(invert_params->tmode_file, sval);

  if (config_lookup_int(&cfg, "max_iter", &ival))
    invert_params->max_iter = (size_t) ival;
  if (config_lookup_float(&cfg, "qdlat_fit_cutoff", &fval))
    {
      data_params->qdlat_fit_cutoff = fval;
      invert_params->qdlat_fit_cutoff = fval;
    }
  if (config_lookup_int(&cfg, "use_weights", &ival))
    invert_params->use_weights = ival;
  if (config_lookup_int(&cfg, "regularize", &ival))
    invert_params->regularize = ival;

  if (config_lookup_float(&cfg, "weight_X", &fval))
    invert_params->weight_X = fval;
  if (config_lookup_float(&cfg, "weight_Y", &fval))
    invert_params->weight_Y = fval;
  if (config_lookup_float(&cfg, "weight_Z", &fval))
    invert_params->weight_Z = fval;
  if (config_lookup_float(&cfg, "weight_F", &fval))
    invert_params->weight_F = fval;
  if (config_lookup_float(&cfg, "weight_DXDT", &fval))
    invert_params->weight_DXDT = fval;
  if (config_lookup_float(&cfg, "weight_DYDT", &fval))
    invert_params->weight_DYDT = fval;
  if (config_lookup_float(&cfg, "weight_DZDT", &fval))
    invert_params->weight_DZDT = fval;
  if (config_lookup_float(&cfg, "weight_DX", &fval))
    invert_params->weight_DX = fval;
  if (config_lookup_float(&cfg, "weight_DY", &fval))
    invert_params->weight_DY = fval;
  if (config_lookup_float(&cfg, "weight_DZ", &fval))
    invert_params->weight_DZ = fval;

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
    invert_params->synth_data = ival;
  if (config_lookup_int(&cfg, "synth_noise", &ival))
    invert_params->synth_noise = ival;
  if (config_lookup_int(&cfg, "synth_nmin", &ival))
    invert_params->synth_nmin = (size_t) ival;

  config_destroy(&cfg);

  return 0;
}

void
print_J_grid(const char * filename, const double t, invert_workspace * w)
{
  FILE *fp = fopen(filename, "w");
  const double r = R_EARTH_KM + 110.0;
  double lon, lat;

  for (lon = -180.0; lon <= 180.0; lon += 5.0)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -89.9; lat <= 89.9; lat += 5.0)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double J[3];

          invert_nonlinear_model_J(0, w->c, t, r, theta, phi, 0, J, w);

          fprintf(fp, "%f %f %.6e %.6e %.6e\n",
                  lon,
                  lat,
                  J[0],
                  J[1],
                  J[2]);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options] sat1.dat sat2.dat ...\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --maxit | -n num_iterations     - number of robust iterations\n");
  fprintf(stderr, "\t --output_file | -o file         - coefficient output file (ASCII)\n");
  fprintf(stderr, "\t --epoch | -e epoch              - model epoch in decimal years\n");
  fprintf(stderr, "\t --print_residuals | -r          - write residuals at each iteration\n");
  fprintf(stderr, "\t --tmin | -b min_time            - minimum data period time in decimal years\n");
  fprintf(stderr, "\t --tmax | -c max_time            - maximum data period time in decimal years\n");
  fprintf(stderr, "\t --print_data | -d               - print data used for MF modeling to output directory\n");
  fprintf(stderr, "\t --print_map | -m                - print spatial data map files to output directory\n");
  fprintf(stderr, "\t --config_file | -C file         - configuration file\n");
} /* print_help() */

int
main(int argc, char *argv[])
{
  int status;
  const char *error_file = "error.txt";
  char *outfile = NULL;
  char *datamap_prefix = "output";
  char *data_prefix = "output";
  char *residual_prefix = "output";
  char *config_file = "INVERT.cfg"; /* default config file */
  invert_workspace *invert_workspace_p;
  invert_parameters invert_params;
  invert_data_workspace *invert_data_p;
  invert_data_parameters data_params;
  gsl_vector *coeffs; /* model coefficients */
  size_t iter = 0;
  size_t maxit = 0;
  double epoch = -1.0;        /* model epoch */
  double tmin = -1.0;         /* minimum time for data in years */
  double tmax = -1.0;         /* maximum time for data in years */
  int nsource = 0;            /* number of data sources (satellites/observatories) */
  int print_data = 0;         /* print data for MF modeling */
  int print_map = 0;          /* print data maps */
  int print_residuals = 0;    /* print residuals at each iteration */
  struct timeval tv0, tv1;
  char buf[MAX_BUFFER];

  data_params.qdlat_fit_cutoff = -1.0;

  invert_init_params(&invert_params);

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "print_residuals", no_argument, NULL, 'r' },
          { "output_file", required_argument, NULL, 'o' },
          { "epoch", required_argument, NULL, 'e' },
          { "maxit", required_argument, NULL, 'n' },
          { "tmin", required_argument, NULL, 'b' },
          { "tmax", required_argument, NULL, 'c' },
          { "print_data", required_argument, NULL, 'd' },
          { "print_map", required_argument, NULL, 'm' },
          { "config_file", required_argument, NULL, 'C' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "b:c:C:de:l:mn:o:r", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
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

          case 'r':
            print_residuals = 1;
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
  parse_config_file(config_file, &invert_params, &data_params);
  fprintf(stderr, "done\n");

  /* check if any command-line arguments should override config file values */
  if (epoch > 0.0)
    invert_params.epoch = epoch;
  if (maxit > 0)
    invert_params.max_iter = maxit;

  status = check_parameters(&invert_params, &data_params);
  if (status)
    exit(1);

  fprintf(stderr, "main: epoch             = %.2f\n", invert_params.epoch);
  fprintf(stderr, "main: radius            = %g [km]\n", invert_params.R);
  fprintf(stderr, "main: num_spatial_modes = %zu\n", invert_params.nspatmodes);

  fprintf(stderr, "main: tmin = %g\n", tmin);
  fprintf(stderr, "main: tmax = %g\n", tmax);
  fprintf(stderr, "main: number of robust iterations = %zu\n", invert_params.max_iter);
  fprintf(stderr, "main: number of data sources = %d\n", nsource);
  fprintf(stderr, "main: number of threads = %d\n", omp_get_max_threads());
  fprintf(stderr, "main: print_residuals = %d\n", print_residuals);
  if (outfile)
    fprintf(stderr, "main: output coefficient file = %s\n", outfile);

  /* allocate data workspace */
  data_params.epoch = invert_params.epoch;
  invert_data_p = invert_data_alloc(nsource, &data_params);

  {
    int satnum = 0;

    while (optind < argc)
      {
        magdata **mdata = &(invert_data_p->mdata[satnum]);

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
    nflag = invert_data_filter_time(tmin, tmax, invert_data_p);
    fprintf(stderr, "done (%zu data flagged)\n", nflag);

    fprintf(stderr, "main: flagging non-fitted components...");
    nflag = invert_data_filter_comp(invert_data_p);
    fprintf(stderr, "done (%zu data flagged)\n", nflag);

    fprintf(stderr, "main: flagging sparse observatory data...");
    nflag = invert_data_filter_observatory(invert_data_p);
    fprintf(stderr, "done (%zu observatories flagged)\n", nflag);
  }

  fprintf(stderr, "main: data epoch = %.2f\n", invert_data_epoch(invert_data_p));
  fprintf(stderr, "main: data tmin  = %.2f\n", satdata_epoch2year(invert_data_p->t0_data));
  fprintf(stderr, "main: data tmax  = %.2f\n", satdata_epoch2year(invert_data_p->t1_data));

  if (print_map)
    {
      /* print spatial coverage maps for each satellite */
      invert_data_map(datamap_prefix, invert_data_p);
    }

  invert_params.nsat = nsource;
  invert_params.invert_data_p = invert_data_p;

  /* allocate invert workspace */
  invert_workspace_p = invert_alloc(&invert_params);
  if (invert_workspace_p == NULL)
    {
      fprintf(stderr, "main: invert_alloc failed\n");
      exit(1);
    }

  fprintf(stderr, "main: number of total parameters:                %zu\n", invert_workspace_p->p);

#if 0 /*XXX*/
  if (invert_params.synth_data)
    {
      fprintf(stderr, "main: replacing with synthetic data...");
      gettimeofday(&tv0, NULL);
      invert_synth_replace(invert_workspace_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
    }
#endif

  /* initialize model parameters */
  invert_init(invert_workspace_p);

  /* print out dataset if requested - do this after invert_init() so
   * spatial weights are computed */
  if (print_data)
    {
      /* print data used for MF modeling for each satellite */
      fprintf(stderr, "main: printing data for MF modeling to %s...", data_prefix);
      invert_data_print(data_prefix, invert_workspace_p->wts_spatial, invert_data_p);
      fprintf(stderr, "done\n");
    }

  /* construct initial guess vector from IGRF */
  coeffs = gsl_vector_alloc(invert_workspace_p->p);
  fprintf(stderr, "main: constructing initial coefficient vector...");
  initial_guess(coeffs, invert_workspace_p);
  fprintf(stderr, "done\n");

  if (print_residuals)
    {
      /* print residuals with initial guess vector */
      fprintf(stderr, "main: printing initial residuals to %s...", residual_prefix);
      gsl_vector_memcpy(invert_workspace_p->c, coeffs);
      invert_residual_print(residual_prefix, 0, invert_workspace_p);
      fprintf(stderr, "done\n");
    }

  gettimeofday(&tv0, NULL);

  status = GSL_CONTINUE;
  while (iter++ < invert_params.max_iter)
    {
      fprintf(stderr, "main: ROBUST ITERATION %zu/%zu\n", iter, invert_params.max_iter);

      status = invert_calc_nonlinear(coeffs, invert_workspace_p);

      if (print_residuals)
        {
          fprintf(stderr, "main: printing residuals to %s...", residual_prefix);
          invert_residual_print(residual_prefix, iter, invert_workspace_p);
          fprintf(stderr, "done\n");
        }

      /* reset workspace for a new iteration */
      invert_reset(invert_workspace_p);
    }

  gettimeofday(&tv1, NULL);

  fprintf(stderr, "main: total time for inversion: %.2f seconds\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: printing J grid...");
  print_J_grid("grid.txt", invert_workspace_p->t0_data, invert_workspace_p);
  fprintf(stderr, "done\n");

#if 0
  /* calculate covariance matrix */
  if (invert_workspace_p->old_fdf == 0)
    {
      fprintf(stderr, "main: calculating covariance matrix...");
      gettimeofday(&tv0, NULL);
      status = invert_covariance(invert_workspace_p->covar, invert_workspace_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (status = %d, %g seconds)\n", status, time_diff(tv0, tv1));

      /* calculate errors in coefficients */
      fprintf(stderr, "main: printing coefficient uncertainties to %s...", error_file);
      gettimeofday(&tv0, NULL);
      invert_print_uncertainties(error_file, invert_workspace_p->covar, invert_workspace_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      /* calculate correlation matrix */
      fprintf(stderr, "main: calculating correlation matrix...");
      gettimeofday(&tv0, NULL);
      status = invert_correlation(invert_workspace_p->covar, invert_workspace_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      printsym_octave(invert_workspace_p->covar, "corr");
    }
#endif

  sprintf(buf, "invert_coeffs.txt");
  fprintf(stderr, "main: writing coefficients to %s...", buf);
  invert_write_ascii(buf, invert_workspace_p->epoch, coeffs, invert_workspace_p);
  fprintf(stderr, "done\n");

  if (outfile)
    {
      fprintf(stderr, "main: writing ASCII coefficients to %s...", outfile);
      /*invert_write(outfile, invert_workspace_p);*/
      invert_write_ascii(outfile, invert_workspace_p->epoch, 0, invert_workspace_p);
      fprintf(stderr, "done\n");
    }

#if 0/*XXX*/
  /* print residual norm between synthetic and computed coefficients */
  if (invert_params.synth_data)
    {
      gsl_vector *g_synth = gsl_vector_alloc(invert_workspace_p->p_int);
      gsl_vector_view g = gsl_vector_subvector(coeffs, 0, invert_workspace_p->p_int);
      double norm;

      /* synthetic internal coefficients */
      invert_synth_g(g_synth, invert_workspace_p);

      /* subtract model internal coefficients */
      gsl_vector_sub(g_synth, &g.vector);

      norm = gsl_blas_dnrm2(g_synth);

      fprintf(stderr, "main: || g_synth - g || = %.12e [nT]\n", norm);

      gsl_vector_free(g_synth);
    }
#endif

  invert_free(invert_workspace_p);
  invert_data_free(invert_data_p);
  gsl_vector_free(coeffs);

  return 0;
} /* main() */
