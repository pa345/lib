/*
 * invert_print.c
 *
 * Usage:
 * ./invert_print [flags]
 *
 * Flags:
 *   -i coef_input_file
 *   -C config_file
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/stat.h>
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

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_track.h>
#include <mainlib/ml_magdata.h>

#include "invert.h"
#include "invert_residual.h"
#include "invert_multifit.c"
#include "invert_synth.h"

#define MAX_BUFFER           2048

static int
check_parameters(const invert_parameters * invert_params)
{
  int s = 0;

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
parse_config_file(const char *filename, invert_parameters *invert_params)
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

  if (config_lookup_float(&cfg, "R", &fval))
    invert_params->R = fval;
  if (config_lookup_int(&cfg, "nfreq", &ival))
    invert_params->nfreq = (size_t) ival;

  if (config_lookup_string(&cfg, "tmode_file", &sval))
    strcpy(invert_params->tmode_file, sval);

  if (config_lookup_int(&cfg, "max_iter", &ival))
    invert_params->max_iter = (size_t) ival;
  if (config_lookup_float(&cfg, "qdlat_fit_cutoff", &fval))
    invert_params->qdlat_fit_cutoff = fval;
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

  config_destroy(&cfg);

  return 0;
}

void
print_spectrum(const char * filename_spatial, const char * filename_grid,
               const gsl_vector * coeffs, invert_workspace * w)
{
  FILE *fp_spatial = fopen(filename_spatial, "w");
  FILE *fp_grid = fopen(filename_grid, "w");
  invert_tmode_workspace * tmode_p = w->tmode_workspace_p;
  invert_smode_workspace * smode_p = w->smode_workspace_p;
  size_t nt, ns, nf;

  fprintf(stderr, "print_spectrum: printing spectrum to %s and %s...", filename_spatial, filename_grid);

  nf = 1;
  fprintf(fp_grid, "# Field %zu: spatial mode\n", nf++);
  fprintf(fp_grid, "# Field %zu: temporal mode\n", nf++);
  fprintf(fp_grid, "# Field %zu: spectral power (nT^2)\n", nf++);

  nf = 1;
  fprintf(fp_spatial, "# Field %zu: spatial mode\n", nf++);
  fprintf(fp_spatial, "# Field %zu: spectral power summed over all temporal modes (nT^2)\n", nf++);

  for (nf = 0; nf < w->nfreq; ++nf)
    {
      fprintf(fp_grid, "# Index: %zu Frequency band: %g [cpd]\n", nf, smode_p->freqs[nf]);
      fprintf(fp_spatial, "# Index: %zu Frequency band: %g [cpd]\n", nf, smode_p->freqs[nf]);

      for (ns = 0; ns < smode_p->nmodes[nf]; ++ns)
        {
          double power_spatial = 0.0;

          for (nt = 0; nt < tmode_p->nmodes[nf]; ++nt)
            {
              size_t idx = invert_coeff_idx(nf, nt, ns, w);
              double c_re = gsl_vector_get(coeffs, idx);
              double c_im = gsl_vector_get(coeffs, idx + w->p_complex);
              double csq = c_re*c_re + c_im*c_im;

              fprintf(fp_grid, "%zu %zu %.12e\n", ns, nt, csq);
              power_spatial += csq;
            }

          if (ns != smode_p->nmodes[nf] - 1)
            fprintf(fp_grid, "\n");

          fprintf(fp_spatial, "%zu %.12e\n", ns, power_spatial);
        }

      fprintf(fp_grid, "\n\n");
      fprintf(fp_spatial, "\n\n");
    }

  fclose(fp_grid);
  fclose(fp_spatial);

  fprintf(stderr, "done\n");
}

/*
print_J_grid()
  Print lat/lon grid of J

Inputs: filename - output file
        t        - timestamp (CDF_EPOCH)
        r        - radius (km)
        w        - workspace
*/

void
print_J_grid(const char * filename, const double t, const double r, invert_workspace * w)
{
  FILE *fp = fopen(filename, "w");
  const double dlon = 2.0;
  const double dlat = 2.0;
  const double lat_min = -89.5;
  const double lat_max = 89.5;
  const size_t nlon = (size_t) (360.0 / dlon) + 1;
  const size_t nlat = (size_t) ((lat_max - lat_min) / dlat) + 1;
  gsl_matrix *Jr = gsl_matrix_alloc(nlon, nlat);
  gsl_matrix *Jt = gsl_matrix_alloc(nlon, nlat);
  gsl_matrix *Jp = gsl_matrix_alloc(nlon, nlat);
  size_t i;

  i = 1;
  fprintf(fp, "# Radius: %.4f [km] (altitude: %.4f [km])\n", r, r - R_EARTH_KM);
  fprintf(fp, "# Timestamp: %ld\n", epoch2timet(t));
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: J_r (uA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: J_t (uA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: J_p (uA/m^2)\n", i++);

#pragma omp parallel for private(i)
  for (i = 0; i < nlon; ++i)
    {
      int thread_id = omp_get_thread_num();
      double lon = -180.0 + i * dlon;
      double phi = lon * M_PI / 180.0;
      size_t j;

      for (j = 0; j < nlat; ++j)
        {
          double lat = lat_min + j * dlat;
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double J[3];

          invert_nonlinear_model_J(w->c, t, r, theta, phi, thread_id, J, w);

          gsl_matrix_set(Jr, i, j, -J[2] * 1.0e-6);
          gsl_matrix_set(Jt, i, j, -J[0] * 1.0e-6);
          gsl_matrix_set(Jp, i, j, J[1] * 1.0e-6);
        }
    }

  /* write file */
  for (i = 0; i < nlon; ++i)
    {
      double lon = -180.0 + i * dlon;
      size_t j;

      for (j = 0; j < nlat; ++j)
        {
          double lat = lat_min + j * dlat;

          fprintf(fp, "%f %f %.6e %.6e %.6e\n",
                  lon,
                  lat,
                  gsl_matrix_get(Jr, i, j),
                  gsl_matrix_get(Jt, i, j),
                  gsl_matrix_get(Jp, i, j));
        }

      fprintf(fp, "\n");
    }

#if 0
  for (lon = -180.0; lon <= 180.0; lon += 2.0)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -89.9; lat <= 89.9; lat += 2.0)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double J[3];

          invert_nonlinear_model_J(0, w->c, t, r, theta, phi, 0, J, w);

          fprintf(fp, "%f %f %.6e %.6e %.6e\n",
                  lon,
                  lat,
                  -J[2],
                  -J[0],
                  J[1]);
        }

      fprintf(fp, "\n");
    }
#endif

  fclose(fp);

  gsl_matrix_free(Jr);
  gsl_matrix_free(Jt);
  gsl_matrix_free(Jp);
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --config_file | -C file         - configuration file\n");
  fprintf(stderr, "\t --input_coef_file | -i file     - input coefficient file (binary)\n");
}

int
main(int argc, char *argv[])
{
  int status;
  char *config_file = "INVERT.cfg"; /* default config file */
  char *input_coef_file = NULL;
  invert_workspace *invert_workspace_p;
  invert_parameters invert_params;
  gsl_vector *coeffs; /* model coefficients */
  struct timeval tv0, tv1;

  invert_init_params(&invert_params);

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "config_file", required_argument, NULL, 'C' },
          { "input_coef_file", required_argument, NULL, 'i' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "C:i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'C':
            config_file = optarg;
            break;

          case 'i':
            input_coef_file = optarg;
            break;

          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  if (input_coef_file == NULL)
    {
      print_help(argv);
      exit(1);
    }

  /* parse configuration file */
  fprintf(stderr, "main: parsing configuration file %s...", config_file);
  parse_config_file(config_file, &invert_params);
  fprintf(stderr, "done\n");

  status = check_parameters(&invert_params);
  if (status)
    exit(1);

  fprintf(stderr, "main: radius            = %g [km]\n", invert_params.R);
  fprintf(stderr, "main: nfreq             = %zu\n", invert_params.nfreq);

  fprintf(stderr, "main: number of robust iterations = %zu\n", invert_params.max_iter);
  fprintf(stderr, "main: number of threads = %d\n", omp_get_max_threads());
  fprintf(stderr, "main: input coef file  = %s\n", input_coef_file);

  invert_params.nsat = 0;
  invert_params.invert_data_p = NULL;

  /* allocate invert workspace */
  invert_workspace_p = invert_alloc(&invert_params);
  if (invert_workspace_p == NULL)
    {
      fprintf(stderr, "main: invert_alloc failed\n");
      exit(1);
    }

  fprintf(stderr, "main: number of total parameters:                %zu\n", invert_workspace_p->p);

  coeffs = gsl_vector_alloc(invert_workspace_p->p);

  fprintf(stderr, "main: reading input coefficients...");
  invert_read(input_coef_file, coeffs);
  fprintf(stderr, "done\n");

  print_spectrum("spectrum.txt", "spectrum_grid.txt", coeffs, invert_workspace_p);

#if 0
  /*XXX*/
  {
    const double ndays = 30.0;
    /*const double t0 = invert_workspace_p->t0_data;*/
    const double t0 = computeEPOCH(2014, 3, 15, 0, 0, 0, 0);
    const double t1 = t0 + ndays * 8.64e7;
    const double dt = 3.6e6;
    char buf[256];
    double t;
    size_t idx = 1;

    gsl_vector_memcpy(invert_workspace_p->c, coeffs);

    for (t = t0; t < t1; t += dt)
      {
        sprintf(buf, "output/J_grid_%04zu.txt", idx++);
        fprintf(stderr, "main: writing %s...", buf);
        gettimeofday(&tv0, NULL);
        print_J_grid(buf, t, invert_workspace_p);
        gettimeofday(&tv1, NULL);
        fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
      }

    exit(1);
  }
#endif

  fprintf(stderr, "main: printing J grid...");
  gsl_vector_memcpy(invert_workspace_p->c, coeffs);
  print_J_grid("grid0.txt", satdata_timet2epoch(1385479822), R_EARTH_KM + 110.0, invert_workspace_p);
  fprintf(stderr, "done\n");

  invert_free(invert_workspace_p);
  gsl_vector_free(coeffs);

  return 0;
}
