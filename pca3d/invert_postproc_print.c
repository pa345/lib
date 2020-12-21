/*
 * invert_postproc_print.c
 *
 * Usage:
 * ./invert_postproc_print [flags]
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
#include <gsl/gsl_matrix.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_track.h>
#include <mainlib/ml_magdata.h>
#include <mainlib/ml_matio.h>

#include "invert.h"
#include "invert_residual.h"
#include "invert_multifit.c"
#include "invert_synth.h"
#include "io.h"
#include "pca3d.h"

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
  char *config_file = "INVERT.cfg"; /* default config file */
  char *input_coef_file = NULL;
  invert_workspace *invert_workspace_p;
  invert_parameters invert_params;
  gsl_vector *coeffs; /* model coefficients */
  struct timeval tv0, tv1;
  size_t p;

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

  fprintf(stderr, "main: nfreq            = %zu\n", invert_params.nfreq);
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

  p = invert_workspace_p->p;
  fprintf(stderr, "main: number of total parameters:                %zu\n", p);

  fprintf(stderr, "main: reading input coefficients from %s...", input_coef_file);
  coeffs = vecread(input_coef_file);
  fprintf(stderr, "done\n");

  if (coeffs->size != p)
    {
      fprintf(stderr, "main: error: coefficient vector does not match workspace (%zu/%zu)\n",
              coeffs->size, p);
      exit(1);
    }

  invert_free(invert_workspace_p);
  gsl_vector_free(coeffs);

  return 0;
}
