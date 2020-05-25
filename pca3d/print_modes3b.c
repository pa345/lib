/*
 * print_modes3b.c
 *
 * Read eigenvector file (left singular vector file) containing modes U
 * (output by stage3b) * and print them in format suitable for plotting
 *
 * ./print_modes [-a altitude(km)]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_bsearch.h>

#include <magfield/magfield.h>
#include <magfield/magfield_eval_complex.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

/* print only magnetic field for Gary */
#define PRINT_B_ONLY          1

/*
fill_magfield()
  Given an eigenmode, fill in the 'qtcoeff', 'pcoeff', and 'qcoeff' arrays
in magfield and compute radial splines to prepare for magfield_eval_complex_J()
and magfield_eval_complex_B() routines

Inputs: u - singular vector (3*N-by-1, where N = nr*nlm_complex)
        w - magfield_eval_complex workspace
*/

int
fill_magfield(const gsl_vector_complex * u, magfield_eval_complex_workspace * w)
{
  const size_t N = w->nr * w->nlm;
  gsl_vector_complex * qtcoeff = gsl_vector_complex_alloc(N);
  gsl_vector_complex * qcoeff = gsl_vector_complex_alloc(N);
  gsl_vector_complex * pcoeff = gsl_vector_complex_alloc(N);
  size_t ir;

  for (ir = 0; ir < w->nr; ++ir)
    {
      size_t l;

      for (l = w->lmin; l <= w->lmax; ++l)
        {
          int M = (int) GSL_MIN(l, w->mmax);
          int m;

          for (m = -M; m <= M; ++m)
            {
              size_t lmidx = magfield_complex_lmidx(l, m, w->mmax);
              size_t uidx = CIDX2(ir, w->nr, lmidx, w->nlm);
              size_t midx = MAG_CMPLX_COEFIDX(ir, lmidx, w);

              gsl_complex qt = gsl_vector_complex_get(u, uidx);
              gsl_complex p = gsl_vector_complex_get(u, uidx + N);
              gsl_complex q = gsl_vector_complex_get(u, uidx + 2*N);

              gsl_vector_complex_set(qtcoeff, midx, qt);
              gsl_vector_complex_set(pcoeff, midx, p);
              gsl_vector_complex_set(qcoeff, midx, q);
            }
        }
    }

  /* compute spline coefficients */
  magfield_eval_complex_init(qtcoeff, qcoeff, pcoeff, w);

  gsl_vector_complex_free(qtcoeff);
  gsl_vector_complex_free(qcoeff);
  gsl_vector_complex_free(pcoeff);

  return 0;
}

int
print_J_modes(const char * prefix, const size_t iband, const double r, const gsl_matrix_complex *U, magfield_params * params)
{
  int status;
  const size_t T = GSL_MIN(U->size2, 20);
  magfield_eval_complex_workspace * w = magfield_eval_complex_alloc(params);
  double lat, lon;
  size_t t;
  char buf[2048];

  sprintf(buf, "%s/band_%zu", prefix, iband);
  fprintf(stderr, "print_J_modes: creating output directory %s...", buf);
  status = mkdir(buf, 0777);
  if (status == -1 && errno != EEXIST)
    {
      fprintf(stderr, "error: unable to create %s: %s\n", buf, strerror(errno));
      exit(1);
    }
  fprintf(stderr, "done\n");

  /* loop over modes t \in [0, T-1] */
  for (t = 0; t < T; ++t)
    {
      gsl_vector_complex_const_view v = gsl_matrix_complex_const_column(U, t);
      FILE *fp;
      size_t i = 1;

      sprintf(buf, "%s/band_%zu/J_mode_%02zu.txt", prefix, iband, t + 1);
      fprintf(stderr, "print_J_modes: writing %s...", buf);

      fp = fopen(buf, "w");

      fprintf(fp, "# Freq band: %zu\n", iband);
      fprintf(fp, "# Radius:    %g [km]\n", r);
      fprintf(fp, "# Altitude:  %g [km]\n", r - R_EARTH_KM);
      fprintf(fp, "# Field %zu: geocentric longitude (deg)\n", i++);
      fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
      fprintf(fp, "# Field %zu: Re J_r\n", i++);
      fprintf(fp, "# Field %zu: Re J_t\n", i++);
      fprintf(fp, "# Field %zu: Re J_p\n", i++);
      fprintf(fp, "# Field %zu: Im J_r\n", i++);
      fprintf(fp, "# Field %zu: Im J_t\n", i++);
      fprintf(fp, "# Field %zu: Im J_p\n", i++);

      /* prepare magfield_eval_complex workspace for this singular vector */
      fill_magfield(&v.vector, w);

      for (lon = -180.0; lon <= 180.0; lon += 1.0)
        {
          double phi = lon * M_PI / 180.0;

          for (lat = -89.5; lat <= 89.5; lat += 1.0)
            {
              double theta = M_PI / 2.0 - lat * M_PI / 180.0;
              complex double J[3];

              magfield_eval_complex_J(r * 1.0e3, theta, phi, J, w);

              fprintf(fp, "%8.4f %8.4f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                      lon,
                      lat,
                      creal(J[0]),
                      creal(J[1]),
                      creal(J[2]),
                      cimag(J[0]),
                      cimag(J[1]),
                      cimag(J[2]));
            }

          fprintf(fp, "\n");
        }

      fclose(fp);

      fprintf(stderr, "done\n");
    }

  magfield_eval_complex_free(w);

  return 0;
}

int
print_B_modes(const char * prefix, const size_t iband, const double r, const gsl_matrix_complex *U, magfield_params * params)
{
  int status;
  const size_t T = GSL_MIN(U->size2, 5);
  magfield_eval_complex_workspace * w = magfield_eval_complex_alloc(params);
  double lat, lon;
  size_t t;
  char buf[2048];

  sprintf(buf, "%s/band_%zu", prefix, iband);
  fprintf(stderr, "print_B_modes: creating output directory %s...", buf);
  status = mkdir(buf, 0777);
  if (status == -1 && errno != EEXIST)
    {
      fprintf(stderr, "error: unable to create %s: %s\n", buf, strerror(errno));
      exit(1);
    }
  fprintf(stderr, "done\n");

  /* loop over modes t \in [0, T-1] */
  for (t = 0; t < T; ++t)
    {
      gsl_vector_complex_const_view v = gsl_matrix_complex_const_column(U, t);
      FILE *fp;
      size_t i = 1;

      sprintf(buf, "%s/band_%zu/B_mode_%02zu.txt", prefix, iband, t + 1);
      fprintf(stderr, "print_B_modes: writing %s...", buf);

      fp = fopen(buf, "w");

      fprintf(fp, "# Freq band: %zu\n", iband);
      fprintf(fp, "# Radius:    %g [km]\n", r);
      fprintf(fp, "# Altitude:  %g [km]\n", r - R_EARTH_KM);
      fprintf(fp, "# Field %zu: geocentric longitude (deg)\n", i++);
      fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
      fprintf(fp, "# Field %zu: Re B_r\n", i++);
      fprintf(fp, "# Field %zu: Re B_t\n", i++);
      fprintf(fp, "# Field %zu: Re B_p\n", i++);
      fprintf(fp, "# Field %zu: Im B_r\n", i++);
      fprintf(fp, "# Field %zu: Im B_t\n", i++);
      fprintf(fp, "# Field %zu: Im B_p\n", i++);

      /* prepare magfield_eval_complex workspace for this singular vector */
      fill_magfield(&v.vector, w);

      for (lon = -180.0; lon <= 180.0; lon += 1.0)
        {
          double phi = lon * M_PI / 180.0;

          for (lat = -89.5; lat <= 89.5; lat += 1.0)
            {
              double theta = M_PI / 2.0 - lat * M_PI / 180.0;
              complex double B[3];

              magfield_eval_complex_B(r * 1.0e3, theta, phi, B, w);

              fprintf(fp, "%8.4f %8.4f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                      lon,
                      lat,
                      creal(B[0]),
                      creal(B[1]),
                      creal(B[2]),
                      cimag(B[0]),
                      cimag(B[1]),
                      cimag(B[2]));
            }

          fprintf(fp, "\n");
        }

      fclose(fp);

      fprintf(stderr, "done\n");
    }

  magfield_eval_complex_free(w);

  return 0;
}

int
main(int argc, char *argv[])
{
  int status;
  char * prefix = "modes_Gary";
  struct timeval tv0, tv1;
  gsl_matrix_complex *U;
  magfield_params params;
  double alt = 110.0;
  int print_J = 1;   /* print J or B */
  size_t iband = 11; /* desired frequency band */
  char buf[2048];
  
  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "prefix", required_argument, NULL, 'p' },
          { "modes_J", no_argument, NULL, 'J' },
          { "modes_B", no_argument, NULL, 'B' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:BJf:p:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            alt = atof(optarg);
            break;

          case 'J':
            print_J = 1;
            break;

          case 'B':
            print_J = 0;
            break;

          case 'f':
            iband = (size_t) atoi(optarg);
            break;

          case 'p':
            prefix = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s [-a altitude (km)] [-J (print J modes)] [-B (print B modes)] [-f freq_band_idx] [-p dir_prefix]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: freq band:  %zu\n", iband);
  fprintf(stderr, "main: altitude: %g [km]\n", alt);

  fprintf(stderr, "main: reading magfield parameters from %s...", PCA3D_STAGE1B_DATA);
  gettimeofday(&tv0, NULL);
  params = pca3d_read_params(PCA3D_STAGE1B_DATA);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_U, iband);

  fprintf(stderr, "main: reading U matrix from %s...", buf);
  gettimeofday(&tv0, NULL);
  U = pca3d_read_matrix_complex(buf);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: creating output directory %s...", prefix);
  status = mkdir(prefix, 0777);
  if (status == -1 && errno != EEXIST)
    {
      fprintf(stderr, "error: unable to create %s: %s\n", prefix, strerror(errno));
      exit(1);
    }
  fprintf(stderr, "done\n");

  if (print_J)
    print_J_modes(prefix, iband, alt + R_EARTH_KM, U, &params);
  else
    print_B_modes(prefix, iband, alt + R_EARTH_KM, U, &params);

  gsl_matrix_complex_free(U);

  return 0;
}
