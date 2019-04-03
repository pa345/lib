/*
 * print.c
 *
 * Print contents of Swarm data files
 *
 * ./print <-i swarm_index_file>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <common/common.h>
#include <msynth/msynth.h>
#include <pomme/pomme.h>

typedef struct
{
  double lt_min;
  double lt_max;
  double ut_min;
  double ut_max;
  double qd_min;
  double qd_max;
  double kp_min;
  double kp_max;
  double alt_min;
  double alt_max;
  msynth_workspace *core_p;
} print_parameters;

/*
print_data()

Inputs: down_sample - number of samples to throw out (>= 1)
                      (ie: if this is 5, every 5th sample is kept and
                       the rest discarded)
        params      - selection parameters
        data        - satellite data input
*/

int
print_data(const int down_sample, const print_parameters *params,
           const satdata_mag *data)
{
  int s = 0;
  size_t i, j;
  f107_workspace *f107_p = f107_alloc(F107_IDX_FILE);

  i = 1;
  printf("# Field %zu: time (UT)\n", i++);
  printf("# Field %zu: time (decimal year)\n", i++);
  printf("# Field %zu: UT (hours)\n", i++);
  printf("# Field %zu: local time (hours)\n", i++);
  printf("# Field %zu: season (day of year)\n", i++);
  printf("# Field %zu: EUVAC\n", i++);
  printf("# Field %zu: longitude (degrees)\n", i++);
  printf("# Field %zu: latitude (degrees)\n", i++);
  printf("# Field %zu: altitude (km)\n", i++);
  printf("# Field %zu: QD latitude (degrees)\n", i++);
  printf("# Field %zu: satellite direction\n", i++);
  printf("# Field %zu: electron density (cm^{-3})\n", i++);
  printf("# Field %zu: scalar field (nT)\n", i++);
  printf("# Field %zu: X field (nT)\n", i++);
  printf("# Field %zu: Y field (nT)\n", i++);
  printf("# Field %zu: Z field (nT)\n", i++);
  printf("# Field %zu: scalar residual (nT)\n", i++);
  printf("# Field %zu: X residual (nT)\n", i++);
  printf("# Field %zu: Y residual (nT)\n", i++);
  printf("# Field %zu: Z residual (nT)\n", i++);
  printf("# Field %zu: X main (nT)\n", i++);
  printf("# Field %zu: Y main (nT)\n", i++);
  printf("# Field %zu: Z main (nT)\n", i++);
  printf("# Field %zu: X crust (nT)\n", i++);
  printf("# Field %zu: Y crust (nT)\n", i++);
  printf("# Field %zu: Z crust (nT)\n", i++);
  printf("# Field %zu: X external (nT)\n", i++);
  printf("# Field %zu: Y external (nT)\n", i++);
  printf("# Field %zu: Z external (nT)\n", i++);

  for (i = 0; i < data->n; i += down_sample)
    {
      time_t unix_time;
      double tyr = satdata_epoch2year(data->t[i]);
      double phi = data->longitude[i] * M_PI / 180.0;
      double lt, ut, euvac;
      double qdlat = data->qdlat[i];
      double B_obs[3], B_main[4], B_crust[3], B_ext[3], B_res[4], B_model[4];

      if (data->flags[i])
        continue; /* nan vector components */

      if (qdlat < params->qd_min || qdlat > params->qd_max)
        continue;

      if (data->altitude[i] < params->alt_min || data->altitude[i] > params->alt_max)
        continue;

      unix_time = satdata_epoch2timet(data->t[i]);

      lt = get_localtime(unix_time, phi);
      if (lt < params->lt_min || lt > params->lt_max)
        continue;

      ut = get_ut(unix_time);
      if (ut < params->ut_min || ut > params->ut_max)
        continue;

      f107_get_euvac(unix_time, &euvac, f107_p);

      B_obs[0] = SATDATA_VEC_X(data->B, i);
      B_obs[1] = SATDATA_VEC_Y(data->B, i);
      B_obs[2] = SATDATA_VEC_Z(data->B, i);

      if (params->core_p != NULL)
        {
          double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
          msynth_eval(tyr, data->r[i], theta, phi, B_main, params->core_p);
        }
      else
        {
          B_main[0] = SATDATA_VEC_X(data->B_main, i);
          B_main[1] = SATDATA_VEC_Y(data->B_main, i);
          B_main[2] = SATDATA_VEC_Z(data->B_main, i);
        }

      B_crust[0] = SATDATA_VEC_X(data->B_crust, i);
      B_crust[1] = SATDATA_VEC_Y(data->B_crust, i);
      B_crust[2] = SATDATA_VEC_Z(data->B_crust, i);

      B_ext[0] = SATDATA_VEC_X(data->B_ext, i);
      B_ext[1] = SATDATA_VEC_Y(data->B_ext, i);
      B_ext[2] = SATDATA_VEC_Z(data->B_ext, i);

      for (j = 0; j < 3; ++j)
        {
          B_model[j] = B_main[j] + B_crust[j] + B_ext[j];
          B_res[j] = B_obs[j] - B_model[j];
        }

      /* compute scalar residual */
      B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);
      B_res[3] = data->F[i] - B_model[3];

      printf("%ld %.8f %.2f %.2f %6.2f %5.1f %10.4f %8.4f %10.4f %10.4f %d %10.4e %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
             unix_time,
             tyr,
             get_ut(unix_time),
             lt,
             get_season(unix_time),
             euvac,
             data->longitude[i],
             data->latitude[i],
             data->altitude[i],
             qdlat,
             satdata_mag_satdir(i, data),
             data->ne[i],
             data->F[i],
             SATDATA_VEC_X(data->B, i),
             SATDATA_VEC_Y(data->B, i),
             SATDATA_VEC_Z(data->B, i),
             B_res[3],
             B_res[0],
             B_res[1],
             B_res[2],
             B_main[0],
             B_main[1],
             B_main[2],
             B_crust[0],
             B_crust[1],
             B_crust[2],
             B_ext[0],
             B_ext[1],
             B_ext[2]);
    }

  f107_free(f107_p);

  return s;
}

int
main(int argc, char *argv[])
{
  satdata_mag *data;
  struct timeval tv0, tv1;
  char *infile = NULL;
  int down_sample = 20;
  print_parameters params;

  params.lt_min = -100.0;
  params.lt_max = 100.0;
  params.ut_min = -100.0;
  params.ut_max = 100.0;
  params.qd_min = -100.0;
  params.qd_max = 100.0;
  params.kp_min = 0.0;
  params.kp_max = 20.0;
  params.alt_min = 0.0;
  params.alt_max = 1000.0;
  params.core_p = NULL;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "lt_min", required_argument, NULL, 'a' },
          { "lt_max", required_argument, NULL, 'b' },
          { "qd_min", required_argument, NULL, 'c' },
          { "qd_max", required_argument, NULL, 'e' },
          { "kp_min", required_argument, NULL, 'f' },
          { "kp_max", required_argument, NULL, 'g' },
          { "ut_min", required_argument, NULL, 'h' },
          { "ut_max", required_argument, NULL, 'j' },
          { "alt_min", required_argument, NULL, 'k' },
          { "alt_max", required_argument, NULL, 'l' },
          { "core_file", required_argument, NULL, 'm' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:d:m:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'd':
            down_sample = atoi(optarg);
            break;

          case 'a':
            params.lt_min = atof(optarg);
            break;

          case 'b':
            params.lt_max = atof(optarg);
            break;

          case 'h':
            params.ut_min = atof(optarg);
            break;

          case 'j':
            params.ut_max = atof(optarg);
            break;

          case 'c':
            params.qd_min = atof(optarg);
            break;

          case 'e':
            params.qd_max = atof(optarg);
            break;

          case 'f':
            params.kp_min = atof(optarg);
            break;

          case 'g':
            params.kp_max = atof(optarg);
            break;

          case 'k':
            params.alt_min = atof(optarg);
            break;

          case 'l':
            params.alt_max = atof(optarg);
            break;

          case 'm':
            fprintf(stderr, "main: reading core field coefficients from %s...", optarg);
            params.core_p = msynth_shc_read(optarg);
            msynth_set(1, 15, params.core_p);
            fprintf(stderr, "done\n");
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i swarm_index_file> [-d down_sample] [--lt_min lt_min] [--lt_max lt_max] [--alt_min alt_min] [--alt_max alt_max] [--qd_min qd_min] [--qd_max qd_max] [--kp_min kp_min] [--kp_max kp_max] [--ut_min ut_min] [--ut_max ut_max] [-m core_shc_file]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);
  fprintf(stderr, "downsample factor = %d\n", down_sample);
  fprintf(stderr, "LT min  = %.2f\n", params.lt_min);
  fprintf(stderr, "LT max  = %.2f\n", params.lt_max);
  fprintf(stderr, "UT min  = %.2f\n", params.ut_min);
  fprintf(stderr, "UT max  = %.2f\n", params.ut_max);
  fprintf(stderr, "QD min  = %.2f [deg]\n", params.qd_min);
  fprintf(stderr, "QD max  = %.2f [deg]\n", params.qd_max);
  fprintf(stderr, "kp min  = %.2f\n", params.kp_min);
  fprintf(stderr, "kp max  = %.2f\n", params.kp_max);
  fprintf(stderr, "alt min = %.2f [km]\n", params.alt_min);
  fprintf(stderr, "alt max = %.2f [km]\n", params.alt_max);

  fprintf(stderr, "Reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = satdata_swarm_read_idx(infile, 0);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  {
    size_t nkp;
    fprintf(stderr, "main: filtering kp outside [%g,%g]...",
            params.kp_min, params.kp_max);
    nkp = satdata_filter_kp(params.kp_min, params.kp_max, data);
    fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
            nkp, data->n, (double)nkp / (double)data->n * 100.0);
  }

  print_data(down_sample, &params, data);

  satdata_mag_free(data);

  return 0;
}
