/*
 * print.c
 *
 * Print contents of DMSP data files
 *
 * ./print <-i dmsp_index_file>
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

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <common/common.h>
#include <common/ellipsoid.h>
#include <common/quat.h>

/* given B in NEC, find component along geodetic vertical */
double
calc_mu3(double r, double theta, double phi, double B[3])
{
  double r_ECEF[3];
  double rhat[3], that[3], phat[3];
  double e_mu[3], e_nu[3], e_phi[3];
  double B3;

  sph2ecef(r, theta, phi, r_ECEF);

  /* compute ECEF spherical basis vectors */
  ecef2sph_basis(r_ECEF, rhat, that, phat);

  /* compute ECEF ellipsoid basis vectors */
#if 0
  ellipsoid_basis(r_ECEF, e_mu, e_nu, e_phi);
#else
  ellipsoid_basis_mu(r_ECEF, WGS84_MU, e_mu, e_nu, e_phi);
#endif

  B3 = -B[0] * vec_dot(that, e_mu) +
        B[1] * vec_dot(phat, e_mu) -
        B[2] * vec_dot(rhat, e_mu);

  /* reverse sign to point inward */
  return -B3;
}

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
  printf("# Field %zu: X VFM (nT)\n", i++);
  printf("# Field %zu: Y VFM (nT)\n", i++);
  printf("# Field %zu: Z VFM (nT)\n", i++);
  printf("# Field %zu: scalar field (nT)\n", i++);
  printf("# Field %zu: X model rotated to S/C fixed (nT)\n", i++);
  printf("# Field %zu: Y model rotated to S/C fixed (nT)\n", i++);
  printf("# Field %zu: Z model rotated to S/C fixed (nT)\n", i++);
  printf("# Field %zu: modeled scalar field (nT)\n", i++);

  for (i = 0; i < data->n; i += down_sample)
    {
      double year = satdata_epoch2year(data->t[i]);
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double *q = &(data->q[4 * i]);
      double *B_VFM = &(data->B_VFM[3 * i]);
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double lt = get_localtime(unix_time, phi);
      double B_model[4], B_model_VFM[3];

      satdata_mag_model(i, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      printf("%ld %f %6.2f %10.4f %10.4f %10.4f %10.4f %2d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
             satdata_epoch2timet(data->t[i]),
             year,
             lt,
             data->longitude[i],
             data->latitude[i],
             data->r[i] - data->R,
             data->qdlat[i],
             satdata_mag_satdir(i, data),
             B_VFM[0],
             B_VFM[1],
             B_VFM[2],
             data->F[i],
             B_model_VFM[0],
             B_model_VFM[1],
             B_model_VFM[2],
             B_model[3]);
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

  while ((c = getopt(argc, argv, "i:d:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'd':
            down_sample = atoi(optarg);
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i dmsp_index_file> [-d down_sample]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);
  fprintf(stderr, "downsample factor = %d\n", down_sample);

  fprintf(stderr, "Reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = satdata_dmsp_read_idx(infile, 0);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  print_data(down_sample, data);

  satdata_mag_free(data);

  return 0;
} /* main() */
