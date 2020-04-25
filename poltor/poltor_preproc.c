/*
 * poltor_preproc.c
 *
 * Pre-process satellite data and save in magdata format. This
 * program should be called with Stage1 data, in order to do track
 * rms test, and along-track gradient calculation
 *
 * Pre-processing steps are:
 * 1. Instrument flags (recommended CHAMP flags except 1 star camera allowed)
 * 2. Track rms test
 * 3. filter for WMM criteria
 * 4. Downsample by factor 20
 * 5. Compute and store along-track gradients
 *
 * The result is an output file in magdata format containing all data point to
 * be used in the modeling. All data points will have a MAGDATA_FLG_FIT_xxx flag
 * set, and other flags will vary.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <omp.h>

#include <flow/flow.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_track.h>
#include <mainlib/ml_euler.h>
#include <mainlib/ml_magdata.h>

#include "poltor.h"

satdata_mag *
read_swarm(const char *filename, const int asmv)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;

  fprintf(stderr, "read_swarm: reading %s...", filename);
  gettimeofday(&tv0, NULL);

  if (asmv)
    data = satdata_swarm_asmv_read_idx(filename, 0);
  else
    data = satdata_swarm_read_idx(filename, 0);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */

  fprintf(stderr, "read_swarm: filtering for instrument flags...");
  nflag = satdata_swarm_filter_instrument(1, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
} /* read_swarm() */

satdata_mag *
read_champ(const char *filename)
{
  satdata_mag *data;
  struct timeval tv0, tv1;
  size_t nflag;

  fprintf(stderr, "read_champ: reading %s...", filename);
  gettimeofday(&tv0, NULL);
  data = satdata_swarm_read_idx(filename, 0);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */
  fprintf(stderr, "read_champ: filtering for instrument flags...");
  nflag = satdata_champ_filter_instrument(1, 0, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
} /* read_champ() */

/*
replace_longitude()
*/

int
replace_longitude(const double ut, satdata_mag *data)
{
  int s = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double lt = get_localtime(unix_time, data->longitude[i] * M_PI / 180.0);

      data->longitude[i] = wrap180(15.0 * (lt - ut));
    }

  return s;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --champ_file      | -c champ_index_file       - CHAMP index file\n");
  fprintf(stderr, "\t --swarm_file      | -s swarm_index_file       - Swarm index file\n");
  fprintf(stderr, "\t --swarm_file2     | -t swarm_index_file2      - Swarm index file 2 (for E/W gradients)\n");
  fprintf(stderr, "\t --swarm_asmv_file | -a swarm_asmv_index_file  - Swarm ASM-V index file\n");
  fprintf(stderr, "\t --downsample      | -d downsample             - downsampling factor\n");
  fprintf(stderr, "\t --gradient_ns     | -g num_samples            - number of samples between N/S gradient points\n");
  fprintf(stderr, "\t --euler_file      | -e euler_file             - Euler angles file\n");
  fprintf(stderr, "\t --euler_file2     | -f euler_file2            - Euler angles file 2 (for E/W gradients)\n");
  fprintf(stderr, "\t --output_file     | -o output_file            - binary output data file (magdata format)\n");
  fprintf(stderr, "\t --config_file     | -C config_file            - configuration file\n");
}

int
main(int argc, char *argv[])
{
  int status;
  char *datamap_prefix = "output/map";
  char *data_prefix = "output/data";
  char *output_file = NULL;
  char *config_file = "PT_preproc.cfg";
  satdata_mag *data = NULL;
  satdata_mag *data2 = NULL;
  magdata *mdata;
  euler_workspace *euler_p = NULL;
  euler_workspace *euler_p2 = NULL;
  struct timeval tv0, tv1;
  track_workspace *track_p = NULL;
  track_workspace *track_p2 = NULL;
  magdata_preprocess_parameters params;
  size_t downsample = 0;     /* downsample factor */
  size_t gradient_ns = 0;    /* number of seconds between N/S gradient points */
  size_t magdata_flags = 0;       /* MAGDATA_GLOBFLG_xxx */
  size_t magdata_flags2 = 0;
  size_t magdata_align_flags = 0; /* ATT_FLG_xxx */
  int flag_vec_rms = 1;

  /* initialize parameters */
  params = magdata_preprocess_default_parameters();

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "swarm_file", required_argument, NULL, 's' },
          { "swarm_file2", required_argument, NULL, 't' },
          { "swarm_asmv_file", required_argument, NULL, 'a' },
          { "champ_file", required_argument, NULL, 'c' },
          { "downsample", required_argument, NULL, 'd' },
          { "output_file", required_argument, NULL, 'o' },
          { "euler_file", required_argument, NULL, 'e' },
          { "euler_file2", required_argument, NULL, 'f' },
          { "config_file", required_argument, NULL, 'C' },
          { "gradient_ns", required_argument, NULL, 'g' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:c:C:d:e:f:g:o:s:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          /* Swarm ASM-V */
          case 'a':
            data = read_swarm(optarg, 0);
            /*magdata_flags = MAGDATA_GLOBFLG_EULER;*/
            flag_vec_rms = 0; /* no NEC data for rms flagging */
            break;

          /* Swarm official */
          case 's':
            data = read_swarm(optarg, 0);
            /*magdata_flags = MAGDATA_GLOBFLG_EULER;*/
            break;

          /* For E/W gradients */
          case 't':
            data2 = read_swarm(optarg, 0);
            /*magdata_flags2 = MAGDATA_GLOBFLG_EULER;*/
            break;

          case 'c':
            data = read_champ(optarg);
            break;

          case 'C':
            config_file = optarg;
            break;

          case 'd':
            downsample = (size_t) atoi(optarg);
            break;

          case 'e':
              fprintf(stderr, "main: reading Euler angles from %s...", optarg);
              euler_p = euler_read(optarg);
              if (!euler_p)
                exit(1);
              fprintf(stderr, "done (%zu sets of angles read)\n", euler_p->n);
              break;

          case 'f':
              fprintf(stderr, "main: reading Euler angles from %s...", optarg);
              euler_p2 = euler_read(optarg);
              if (!euler_p2)
                exit(1);
              fprintf(stderr, "done (%zu sets of angles read)\n", euler_p2->n);
              break;

          case 'g':
            gradient_ns = (size_t) atoi(optarg);
            break;

          case 'o':
            output_file = optarg;
            break;

          default:
            break;
        }
    }

  /* parse configuration file */
  fprintf(stderr, "main: parsing config file %s...", config_file);
  magdata_preprocess_parse(config_file, &params);
  fprintf(stderr, "done\n");

  /* replace config values with command-line arguments */
  if (downsample > 0)
    params.downsample = downsample;
  if (gradient_ns > 0)
    params.gradient_ns = gradient_ns;
  if (flag_vec_rms == 0)
    {
      params.rms_thresh[0] = -1.0;
      params.rms_thresh[1] = -1.0;
      params.rms_thresh[2] = -1.0;
    }

  /* check parameters */
  status = magdata_preprocess_check(&params);
  if (status)
    exit(1);

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: downsample       = %zu\n", params.downsample);
  fprintf(stderr, "main: gradient_ns      = %zu [s]\n", params.gradient_ns);
  fprintf(stderr, "main: rms X threshold  = %.1f [nT]\n", params.rms_thresh[0]);
  fprintf(stderr, "main: rms Y threshold  = %.1f [nT]\n", params.rms_thresh[1]);
  fprintf(stderr, "main: rms Z threshold  = %.1f [nT]\n", params.rms_thresh[2]);
  fprintf(stderr, "main: rms F threshold  = %.1f [nT]\n", params.rms_thresh[3]);
  fprintf(stderr, "main: LT minimum       = %.1f\n", params.min_LT);
  fprintf(stderr, "main: LT maximum       = %.1f\n", params.max_LT);
  fprintf(stderr, "main: Euler LT minimum = %.1f\n", params.euler_min_LT);
  fprintf(stderr, "main: Euler LT maximum = %.1f\n", params.euler_max_LT);

  if (euler_p)
    {
      fprintf(stderr, "main: rotating VFM measurements with new Euler angles...");
      euler_apply(data, euler_p);
      fprintf(stderr, "done\n");
    }

  if (euler_p2 && data2)
    {
      fprintf(stderr, "main: rotating VFM measurements with new Euler angles for satellite 2...");
      euler_apply(data2, euler_p2);
      fprintf(stderr, "done\n");
    }

  fprintf(stderr, "main: === PREPROCESSING SATELLITE 1 ===\n");
  track_p = magdata_preprocess(&params, magdata_flags, data);

  if (data2)
    {
      fprintf(stderr, "main: === PREPROCESSING SATELLITE 2 ===\n");
      track_p2 = magdata_preprocess(&params, magdata_flags2, data2);
    }

  fprintf(stderr, "main: computing vector residuals...");
  gettimeofday(&tv0, NULL);
  mdata = magdata_preprocess_fill(magdata_flags, data, track_p, magdata_flags2, data2, track_p2, &params);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu residuals copied, %g seconds)\n",
          mdata->n, time_diff(tv0, tv1));

  /* free data after copying arrays to free up memory */
  satdata_mag_free(data);

  if (data2)
    satdata_mag_free(data2);

#if 1
  fprintf(stderr, "main: recomputing longitude with LT...");
  magdata_replace_phi_LT(12.0, mdata);
  fprintf(stderr, "done\n");
#endif

  magdata_init(mdata);

  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: writing data to %s...", data_prefix);
  magdata_print(data_prefix, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing data map to %s...", datamap_prefix);
  magdata_map(datamap_prefix, MAGDATA_FLG_FIT_MF, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: satellite rmin = %.1f (%.1f) [km]\n",
          mdata->rmin, mdata->rmin - mdata->R);
  fprintf(stderr, "main: satellite rmax = %.1f (%.1f) [km]\n",
          mdata->rmax, mdata->rmax - mdata->R);

  if (output_file)
    {
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");
    }

  magdata_free(mdata);
  track_free(track_p);

  if (euler_p)
    euler_free(euler_p);

  if (euler_p2)
    euler_free(euler_p2);

  if (track_p2)
    track_free(track_p2);

  return 0;
} /* main() */
