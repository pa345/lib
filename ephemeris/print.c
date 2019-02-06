/*
 * print.c
 * Patrick Alken
 *
 * Usage: ./print [-b bowman_ephemeris_gz_file]
 * 
 * This program reads a Bowman ephemeris file, and prints
 * the data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>

#include <satdata/satdata.h>

#include <common/common.h>
#include <common/eci.h>

#include "eph_data.h"

int
print_eph(const eph_data *data)
{
  size_t i;

  i = 1;
  printf("# Field %zu: time (decimal year)\n", i++);
  printf("# Field %zu: ECI X (m)\n", i++);
  printf("# Field %zu: ECI Y (m)\n", i++);
  printf("# Field %zu: ECI Z (m)\n", i++);
  printf("# Field %zu: ECEF X (m)\n", i++);
  printf("# Field %zu: ECEF Y (m)\n", i++);
  printf("# Field %zu: ECEF Z (m)\n", i++);
  printf("# Field %zu: ECI VX (m/s)\n", i++);
  printf("# Field %zu: ECI VY (m/s)\n", i++);
  printf("# Field %zu: ECI VZ (m/s)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      double r_ECI[3], r_ECEF[3];

      r_ECI[0] = data->X[i];
      r_ECI[1] = data->Y[i];
      r_ECI[2] = data->Z[i];

      eci2ecef(satdata_epoch2timet(data->t[i]), r_ECI, r_ECEF);

      printf("%.12f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
             satdata_epoch2year(data->t[i]),
             data->X[i],
             data->Y[i],
             data->Z[i],
             r_ECEF[0],
             r_ECEF[1],
             r_ECEF[2],
             data->VX[i],
             data->VY[i],
             data->VZ[i]);
    }

  return 0;
}

int
main(int argc, char *argv[])
{
  eph_data *data = NULL;
  int c;
  struct timeval tv0, tv1;

  while ((c = getopt(argc, argv, "b:")) != (-1))
    {
      switch (c)
        {
          case 'b':
            fprintf(stderr, "main: reading Bowman ephemerides from %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = eph_data_read_bowman(optarg);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu read, %g seconds)\n", data->n, time_diff(tv0, tv1));
        }
    }

  if (data == NULL)
    {
      fprintf(stderr, "Usage: %s [-b <bowman_gz_file>]\n", argv[0]);
      exit(1);
    }

  print_eph(data);

  eph_data_free(data);

  return 0;
} /* main() */
