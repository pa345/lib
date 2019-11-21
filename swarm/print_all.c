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
#include <gsl/gsl_statistics.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_pomme.h>
#include <mainlib/ml_euler.h>

/*
print_data()

Inputs: down_sample - number of samples to throw out (>= 1)
                      (ie: if this is 5, every 5th sample is kept and
                       the rest discarded)
        data        - satellite data input
*/

int
print_data(const int down_sample, const satdata_mag *data)
{
  int s = 0;
  size_t i;

  i = 1;
  printf("# Field %zu: time (UT)\n", i++);
  printf("# Field %zu: time (decimal year)\n", i++);
  printf("# Field %zu: local time (hours)\n", i++);
  printf("# Field %zu: longitude (degrees)\n", i++);
  printf("# Field %zu: latitude (degrees)\n", i++);
  printf("# Field %zu: altitude (km)\n", i++);
  printf("# Field %zu: QD latitude (degrees)\n", i++);
  printf("# Field %zu: satellite direction\n", i++);
  printf("# Field %zu: scalar field (nT)\n", i++);
  printf("# Field %zu: NEC X component (nT)\n", i++);
  printf("# Field %zu: NEC Y component (nT)\n", i++);
  printf("# Field %zu: NEC Z component (nT)\n", i++);
  printf("# Field %zu: main field X (nT)\n", i++);
  printf("# Field %zu: main field Y (nT)\n", i++);
  printf("# Field %zu: main field Z (nT)\n", i++);
  printf("# Field %zu: crustal field X (nT)\n", i++);
  printf("# Field %zu: crustal field Y (nT)\n", i++);
  printf("# Field %zu: crustal field Z (nT)\n", i++);
  printf("# Field %zu: external field X (nT)\n", i++);
  printf("# Field %zu: external field Y (nT)\n", i++);
  printf("# Field %zu: external field Z (nT)\n", i++);

  for (i = 0; i < data->n; i += down_sample)
    {
      double year = satdata_epoch2year(data->t[i]);
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double phi = data->longitude[i] * M_PI / 180.0;
      double lt = get_localtime(unix_time, phi);
      double B_model[4], B_nec[3];

      B_model[0] = SATDATA_VEC_X(data->B_main, i) +
                   SATDATA_VEC_X(data->B_crust, i) +
                   SATDATA_VEC_X(data->B_ext, i);
      B_model[1] = SATDATA_VEC_Y(data->B_main, i) +
                   SATDATA_VEC_Y(data->B_crust, i) +
                   SATDATA_VEC_Y(data->B_ext, i);
      B_model[2] = SATDATA_VEC_Z(data->B_main, i) +
                   SATDATA_VEC_Z(data->B_crust, i) +
                   SATDATA_VEC_Z(data->B_ext, i);
      B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

      B_nec[0] = SATDATA_VEC_X(data->B, i);
      B_nec[1] = SATDATA_VEC_Y(data->B, i);
      B_nec[2] = SATDATA_VEC_Z(data->B, i);

      printf("%ld %f %6.2f %10.4f %10.4f %10.4f %10.4f %2d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
             satdata_epoch2timet(data->t[i]),
             year,
             lt,
             data->longitude[i],
             data->latitude[i],
             data->r[i] - R_EARTH_KM,
             data->qdlat[i],
             satdata_mag_satdir(i, data),
             data->F[i],
             B_nec[0],
             B_nec[1],
             B_nec[2],
             SATDATA_VEC_X(data->B_main, i),
             SATDATA_VEC_Y(data->B_main, i),
             SATDATA_VEC_Z(data->B_main, i),
             SATDATA_VEC_X(data->B_crust, i),
             SATDATA_VEC_Y(data->B_crust, i),
             SATDATA_VEC_Z(data->B_crust, i),
             SATDATA_VEC_X(data->B_ext, i),
             SATDATA_VEC_Y(data->B_ext, i),
             SATDATA_VEC_Z(data->B_ext, i));
    }

  return s;
}

int
main(int argc, char *argv[])
{
  satdata_mag *data;
  struct timeval tv0, tv1;
  int c;
  char *infile = NULL;
  int down_sample = 20;
  euler_workspace *euler_p = NULL;

  while ((c = getopt(argc, argv, "i:d:e:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'd':
            down_sample = atoi(optarg);
            break;

          case 'e':
              fprintf(stderr, "main: reading Euler angles from %s...", optarg);
              euler_p = euler_read(optarg);
              if (!euler_p)
                exit(1);
              fprintf(stderr, "done (%zu sets of angles read)\n", euler_p->n);
              break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i swarm_index_file> [-d down_sample] [-e euler_file]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);
  fprintf(stderr, "downsample factor = %d\n", down_sample);

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

  if (euler_p)
    {
      fprintf(stderr, "main: rotating VFM measurements with new Euler angles...");
      euler_apply(data, euler_p);
      fprintf(stderr, "done\n");
    }

  print_data(down_sample, data);

  satdata_mag_free(data);

  if (euler_p)
    euler_free(euler_p);

  return 0;
} /* main() */
