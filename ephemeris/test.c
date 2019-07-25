/*
 * test.c
 * Patrick Alken
 *
 * Usage: ./test [-b bowman_ephemeris_gz_file]
 * 
 * This program reads a Bowman ephemeris file, tests
 * lat/lon matches the ECI positions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>

#include <gsl/gsl_test.h>

#include <satdata/satdata.h>
#include <sofa/sofa.h>
#include <spice/spice.h>

#include <common/common.h>
#include <common/eci.h>

#include "eph.h"
#include "eph_data.h"

/*
 * 0 = original ECI code
 * 1 = SPICE library
 * 2 = SOFA library
 */
#define ECI_LIB   0

int
test_geo(const double eps_lat, const double eps_lon,
         double *max_dlat, double *max_dlon, const eph_data *data)
{
  int s = 0;
  size_t i;
  double maxdlat = -1.0, maxdlon = -1.0;
  sofa_workspace *sofa_p = sofa_alloc();
  spice_workspace *spice_p = spice_alloc();

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double lat, lon;
      double r_ECI[3], r_sph[3];
      double eps;

      r_ECI[0] = data->X[i];
      r_ECI[1] = data->Y[i];
      r_ECI[2] = data->Z[i];

      /* convert ECI position to geocentric spherical coordinates */
#if ECI_LIB == 0
      eci2sph_pos(unix_time, r_ECI, r_sph);
#elif ECI_LIB == 1
      spice_c2sph(unix_time, r_ECI, r_sph);
#else
      sofa_c2sph(unix_time, r_ECI, r_sph, sofa_p);
#endif

      /* extract geocentric latitude and longitude in deg */
      lon = r_sph[2] * 180.0 / M_PI;
      lat = 90.0 - r_sph[1] * 180.0 / M_PI;

      /* Bruce gives lat/lon to 3 decimal places, so these tolerances are chosen accordingly */

      eps = wrap180(lon - data->longitude[i]);
      /*gsl_test_abs(eps, 0.0, eps_lon, "longitude discrepency");*/
      if (fabs(eps) > maxdlon)
        maxdlon = fabs(eps);

      eps = lat - data->latitude[i];
      /*gsl_test_abs(eps, 0.0, eps_lat, "latitude discrepency");*/
      if (fabs(eps) > maxdlat)
        maxdlat = fabs(eps);
    }

  *max_dlon = maxdlon;
  *max_dlat = maxdlat;

  sofa_free(sofa_p);
  spice_free(spice_p);

  return s;
} /* test_geo() */

/* test interpolation at ephemeris timestamps produce correct ECI values */
int
test_interp(const eph_data *data)
{
  int s = 0;
  eph_workspace *w = eph_alloc(data);
  size_t i;
  const double tol = 1.0e-8;

  for (i = 0; i < data->n; ++i)
    {
      double t = data->t[i];
      double pos[3], vel[3];

      eph_interp(t, pos, vel, w);

      gsl_test_rel(pos[0], data->X[i], tol, "X");
      gsl_test_rel(pos[1], data->Y[i], tol, "Y");
      gsl_test_rel(pos[2], data->Z[i], tol, "Z");

      gsl_test_rel(vel[0], data->VX[i], tol, "VX");
      gsl_test_rel(vel[1], data->VY[i], tol, "VY");
      gsl_test_rel(vel[2], data->VZ[i], tol, "VZ");
    }

  eph_free(w);

  return s;
} /* test_interp() */

int
main(int argc, char *argv[])
{
  eph_data *data = NULL;
  int c;
  struct timeval tv0, tv1;
  int test_bowman = 0;

  while ((c = getopt(argc, argv, "b:t:")) != (-1))
    {
      switch (c)
        {
          case 'b':
            fprintf(stderr, "main: reading Bowman ephemerides from %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = eph_data_read_bowman(optarg);
            gettimeofday(&tv1, NULL);
            if (!data)
              break;
            fprintf(stderr, "done (%zu read, %g seconds)\n", data->n, time_diff(tv0, tv1));
            test_bowman = 1;
            break;

          case 't':
            fprintf(stderr, "main: reading TENA ephemerides from %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = eph_data_read_tena(optarg);
            gettimeofday(&tv1, NULL);
            if (!data)
              break;
            fprintf(stderr, "done (%zu read, %g seconds)\n", data->n, time_diff(tv0, tv1));
            break;
        }
    }

  if (data == NULL)
    {
      fprintf(stderr, "Usage: %s [-b <bowman_gz_file>] [-t tena_file]\n", argv[0]);
      exit(1);
    }

  if (test_bowman)
    {
      const double eps_lat = 0.00055; /* latitude tolerance in deg */
      const double eps_lon = 0.00055; /* longitude tolerance in deg */
      double maxdlat, maxdlon;

      fprintf(stderr, "main: testing lat/lon values...");
      test_geo(eps_lat, eps_lon, &maxdlat, &maxdlon, data);
      fprintf(stderr, "done (max dlat = %f [deg], max dlon = %f [deg]\n",
              maxdlat, maxdlon);
    }

  fprintf(stderr, "main: testing Hermite interpolation...");
  test_interp(data);
  fprintf(stderr, "done\n");

  eph_data_free(data);

  exit (gsl_test_summary());
}
