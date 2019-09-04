/*
 * print_efi.c
 *
 * Print contents of Swarm LP data files
 *
 * ./print <-i swarm_efi_file>
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

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_apex.h>
#include <mainlib/ml_common.h>

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
           const satdata_efi *data)
{
  int s = 0;
  size_t i;
  f107_workspace *f107_p = f107_alloc(F107_IDX_FILE);
  double epoch = satdata_epoch2year(data->t[0]);
  apex_workspace *apex_p = apex_alloc();

  i = 1;
  printf("# Field %zu: time (CDF_EPOCH)\n", i++);
  printf("# Field %zu: time (decimal year)\n", i++);
  printf("# Field %zu: UT (hours)\n", i++);
  printf("# Field %zu: local time (hours)\n", i++);
  printf("# Field %zu: season (day of year)\n", i++);
  printf("# Field %zu: EUVAC\n", i++);
  printf("# Field %zu: longitude (degrees)\n", i++);
  printf("# Field %zu: latitude (degrees)\n", i++);
  printf("# Field %zu: radius (km)\n", i++);
  printf("# Field %zu: QD latitude (degrees)\n", i++);
  printf("# Field %zu: E_x (mV/m)\n", i++);
  printf("# Field %zu: E_y (mV/m)\n", i++);
  printf("# Field %zu: E_z (mV/m)\n", i++);
  printf("# Field %zu: v_x (m/s)\n", i++);
  printf("# Field %zu: v_y (m/s)\n", i++);
  printf("# Field %zu: v_z (m/s)\n", i++);
  printf("# Field %zu: ion temperature (K)\n", i++);
  printf("# Field %zu: satellite direction\n", i++);

  for (i = 0; i < data->n; i += down_sample)
    {
      time_t unix_time;
      double phi = data->longitude[i] * M_PI / 180.0;
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double lt, ut, euvac;
      double qdlat, alon, alat;

      if (data->flags[i] > 20)
        continue;

      apex_transform(epoch, theta, phi, data->r[i], &alon, &alat, &qdlat,
                     NULL, NULL, NULL, apex_p);

      if (qdlat < params->qd_min || qdlat > params->qd_max)
        continue;

      if (data->r[i] - R_EARTH_KM < params->alt_min || data->r[i] - R_EARTH_KM > params->alt_max)
        continue;

      unix_time = satdata_epoch2timet(data->t[i]);

      lt = get_localtime(unix_time, phi);
      if (lt < params->lt_min || lt > params->lt_max)
        continue;

      ut = get_ut(unix_time);
      if (ut < params->ut_min || ut > params->ut_max)
        continue;

      f107_get_euvac(unix_time, &euvac, f107_p);

      printf("%f %.8f %6.2f %6.2f %6.2f %5.1f %10.4f %8.4f %10.4f %10.4f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.1f %d\n",
             data->t[i],
             satdata_epoch2year(data->t[i]),
             get_ut(unix_time),
             lt,
             get_season(unix_time),
             euvac,
             data->longitude[i],
             data->latitude[i],
             data->r[i],
             qdlat,
             SATDATA_VEC_X(data->E, i),
             SATDATA_VEC_Y(data->E, i),
             SATDATA_VEC_Z(data->E, i),
             SATDATA_VEC_X(data->v_ion, i),
             SATDATA_VEC_Y(data->v_ion, i),
             SATDATA_VEC_Z(data->v_ion, i),
             data->Ti[i],
             satdata_efi_satdir(i, data));
    }

  f107_free(f107_p);
  apex_free(apex_p);

  return s;
}

/*
print_data2()

Inputs: down_sample - number of samples to throw out (>= 1)
                      (ie: if this is 5, every 5th sample is kept and
                       the rest discarded)
        params      - selection parameters
        data        - satellite data input
*/

int
print_data2(const int down_sample, const print_parameters *params,
           const satdata_efi *data)
{
  int s = 0;
  size_t i;
  f107_workspace *f107_p = f107_alloc(F107_IDX_FILE);
  double epoch = satdata_epoch2year(data->t[0]);
  apex_workspace *apex_p = apex_alloc();
  double qdlat_prev = 0.0;
  const double alpha = 1.0e-2;
  double mean_Ey = SATDATA_VEC_Y(data->E, 0);

  i = 1;
  printf("# Field %zu: time (UT)\n", i++);
  printf("# Field %zu: time (decimal year)\n", i++);
  printf("# Field %zu: UT (hours)\n", i++);
  printf("# Field %zu: local time (hours)\n", i++);
  printf("# Field %zu: season (day of year)\n", i++);
  printf("# Field %zu: EUVAC\n", i++);
  printf("# Field %zu: longitude (degrees)\n", i++);
  printf("# Field %zu: latitude (degrees)\n", i++);
  printf("# Field %zu: radius (km)\n", i++);
  printf("# Field %zu: QD latitude (degrees)\n", i++);
  printf("# Field %zu: E_x (mV/m)\n", i++);
  printf("# Field %zu: E_y (mV/m)\n", i++);
  printf("# Field %zu: E_z (mV/m)\n", i++);
  printf("# Field %zu: v_x (m/s)\n", i++);
  printf("# Field %zu: v_y (m/s)\n", i++);
  printf("# Field %zu: v_z (m/s)\n", i++);
  printf("# Field %zu: moving average E_y (mV/m)\n", i++);
  printf("# Field %zu: ion temperature (K)\n", i++);
  printf("# Field %zu: satellite direction\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time;
      double phi = data->longitude[i] * M_PI / 180.0;
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double lt, euvac;
      double qdlat, alon, alat;

      apex_transform(epoch, theta, phi, data->r[i], &alon, &alat, &qdlat,
                     NULL, NULL, NULL, apex_p);

      if (qdlat_prev * qdlat < 0.0 && fabs(qdlat) < 5.0)
        {
          /* equator crossing */

          if (data->flags[i] > 20)
            continue;

          unix_time = satdata_epoch2timet(data->t[i]);
          lt = get_localtime(unix_time, phi);

          if (lt < 6.0 || lt > 18.0)
            continue;

          f107_get_euvac(unix_time, &euvac, f107_p);

          printf("%ld %.8f %6.2f %6.2f %6.2f %5.1f %10.4f %8.4f %10.4f %10.4f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.1f %d\n",
                 unix_time,
                 satdata_epoch2year(data->t[i]),
                 get_ut(unix_time),
                 lt,
                 get_season(unix_time),
                 euvac,
                 data->longitude[i],
                 data->latitude[i],
                 data->r[i],
                 qdlat,
                 SATDATA_VEC_X(data->E, i),
                 SATDATA_VEC_Y(data->E, i),
                 SATDATA_VEC_Z(data->E, i),
                 SATDATA_VEC_X(data->v_ion, i),
                 SATDATA_VEC_Y(data->v_ion, i),
                 SATDATA_VEC_Z(data->v_ion, i),
                 mean_Ey,
                 data->Ti[i],
                 satdata_efi_satdir(i, data));

        }

      qdlat_prev = qdlat;

      mean_Ey = alpha * SATDATA_VEC_Y(data->E, i) + (1.0 - alpha) * mean_Ey;
    }

  f107_free(f107_p);
  apex_free(apex_p);

  return s;
}

int
main(int argc, char *argv[])
{
  satdata_efi *data;
  struct timeval tv0, tv1;
  char *infile = NULL;
  int down_sample = 1;
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
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:d:", long_options, &option_index);
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

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i swarm_index_file> [-d down_sample] [--lt_min lt_min] [--lt_max lt_max] [--alt_min alt_min] [--alt_max alt_max] [--qd_min qd_min] [--qd_max qd_max] [--kp_min kp_min] [--kp_max kp_max] [--ut_min ut_min] [--ut_max ut_max]\n",
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

  data = satdata_swarm_efi_read_idx(infile, 0);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  print_data2(down_sample, &params, data);

  satdata_efi_free(data);

  return 0;
} /* main() */
