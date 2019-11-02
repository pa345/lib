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
print_modes(const char * prefix, const double freq, const double r, const double r_B, const gsl_matrix_complex *U, magfield_params * params)
{
  const size_t T = GSL_MIN(U->size2, 20);
  magfield_eval_complex_workspace * w = magfield_eval_complex_alloc(params);
  double lat, lon;
  size_t t;

  /* loop over modes t \in [0, T-1] */
  for (t = 0; t < T; ++t)
    {
      gsl_vector_complex_const_view v = gsl_matrix_complex_const_column(U, t);
      FILE *fp;
      char buf[2048];
      size_t i = 1;

      sprintf(buf, "%s/%g_cpd/mode_%02zu.txt", prefix, freq, t + 1);
      fprintf(stderr, "print_modes: writing %s...", buf);

      fp = fopen(buf, "w");

      fprintf(fp, "# Frequency: %g [cpd]\n", freq);
#if !PRINT_B_ONLY
      fprintf(fp, "# Altitude for J: %g [km]\n", r - R_EARTH_KM);
#endif
      fprintf(fp, "# Radius for B: %g [km]\n", r_B);
      fprintf(fp, "# Field %zu: geocentric longitude (deg)\n", i++);
      fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
#if !PRINT_B_ONLY
      fprintf(fp, "# Field %zu: Re J_r\n", i++);
      fprintf(fp, "# Field %zu: Re J_t\n", i++);
      fprintf(fp, "# Field %zu: Re J_p\n", i++);
      fprintf(fp, "# Field %zu: Im J_r\n", i++);
      fprintf(fp, "# Field %zu: Im J_t\n", i++);
      fprintf(fp, "# Field %zu: Im J_p\n", i++);
#endif
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
              complex double J[3], B[3];

#if !PRINT_B_ONLY
              magfield_eval_complex_J(r, theta, phi, J, w);
              magfield_eval_complex_B(r_B, theta, phi, B, w);

              fprintf(fp, "%8.4f %8.4f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                      lon,
                      lat,
                      creal(J[0]),
                      creal(J[1]),
                      creal(J[2]),
                      cimag(J[0]),
                      cimag(J[1]),
                      cimag(J[2]),
                      creal(B[0]),
                      creal(B[1]),
                      creal(B[2]),
                      cimag(B[0]),
                      cimag(B[1]),
                      cimag(B[2]));
#else
              magfield_eval_complex_B(r_B, theta, phi, B, w);

              fprintf(fp, "%8.4f %8.4f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                      lon,
                      lat,
                      creal(B[0]),
                      creal(B[1]),
                      creal(B[2]),
                      cimag(B[0]),
                      cimag(B[1]),
                      cimag(B[2]));
#endif
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
  char * prefix = "modes_Gary";
  double freq = 1.0; /* desired frequency in cpd */
  pca3d_fft_data data;
  struct timeval tv0, tv1;
  gsl_matrix_complex *U;
  magfield_params params;
  double alt = 110.0;
  double alt_B = 0.0;
  size_t ifreq;      /* index of desired frequency */
  char buf[2048];
  
  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "prefix", required_argument, NULL, 'p' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:f:p:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            alt = atof(optarg);
            break;

          case 'b':
            alt_B = atof(optarg);
            break;

          case 'f':
            freq = atof(optarg);
            break;

          case 'p':
            prefix = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s [-a J altitude (km)] [-b B altitudde (km)] [-f freq (cpd)] [-p dir_prefix]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: frequency: %g [cpd]\n", freq);
  fprintf(stderr, "main: altitude:  %g [km]\n", alt);

  fprintf(stderr, "main: reading %s...", PCA3D_STAGE2B_FFT_DATA_LIGHT);
  gettimeofday(&tv0, NULL);
  data = pca3d_read_fft_data2(PCA3D_STAGE2B_FFT_DATA_LIGHT);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: reading magfield parameters from %s...", PCA3D_STAGE1B_DATA);
  gettimeofday(&tv0, NULL);
  params = pca3d_read_params(PCA3D_STAGE1B_DATA);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  ifreq = (size_t) (freq * data.window_size);
  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_U, ifreq);

  fprintf(stderr, "main: reading U matrix from %s...", buf);
  gettimeofday(&tv0, NULL);
  U = pca3d_read_matrix_complex(buf);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  print_modes(prefix, freq, alt + R_EARTH_KM, alt_B + R_EARTH_KM, U, &params);

  gsl_matrix_complex_free(U);

  return 0;
}
