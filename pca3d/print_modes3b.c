/*
 * print_modes3b.c
 *
 * Read eigenvector file (left singular vector file) containing modes U
 * (output by stage3b) * and print them in format suitable for plotting
 *
 * ./print_modes [-i U_data_file]
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <common/common.h>
#include <common/bsearch.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

/*
fill_magfield()
  Given an eigenmode, fill in the 'qtcoeff', 'pcoeff', and 'qcoeff' arrays
in magfield and compute radial splines to prepare for magfield_eval_J()
and magfield_eval_B() routines

Inputs: u - singular vector (3*N-by-T, where N = nr*nlm)
        w - magfield workspace
*/

int
fill_magfield(const gsl_vector_complex * u, magfield_workspace * w)
{
  const size_t N = w->nr * w->nlm;
  size_t ir;

  for (ir = 0; ir < w->nr; ++ir)
    {
      size_t l;

      for (l = w->lmin; l <= w->lmax; ++l)
        {
          int M = (int) GSL_MIN(l, w->mmax);
          int m;

          for (m = 0; m <= M; ++m)
            {
              size_t lmidx = magfield_lmidx(l, m, w->mmax);
              size_t uidx = CIDX2(ir, w->nr, lmidx, w->nlm);
              size_t midx = MAG_COEFIDX(ir, lmidx, w);

              gsl_complex qt = gsl_vector_complex_get(u, uidx);
              gsl_complex p = gsl_vector_complex_get(u, uidx + N);
              gsl_complex q = gsl_vector_complex_get(u, uidx + 2*N);

              gsl_vector_complex_set(w->qtcoeff, midx, qt);
              gsl_vector_complex_set(w->pcoeff, midx, p);
              gsl_vector_complex_set(w->qcoeff, midx, q);
            }
        }
    }

  /* compute spline coefficients */
  magfield_spline_qt_calc(w);
  magfield_spline_p_calc(w);
  magfield_spline_q_calc(w);

  w->field_computed = 1;

  return 0;
}

int
print_modes(const double freq, const double r, const gsl_matrix_complex *U, magfield_workspace * w)
{
  const size_t T = U->size2;
  double lat, lon;
  size_t t;

  /* loop over modes t \in [0, T-1] */
  for (t = 0; t < T; ++t)
    {
      gsl_vector_complex_const_view v = gsl_matrix_complex_const_column(U, t);
      FILE *fp;
      char buf[2048];
      size_t i = 1;

      sprintf(buf, "plots/J_%02zu.txt", t + 1);
      fprintf(stderr, "print_modes: writing %s...", buf);

      fp = fopen(buf, "w");

      fprintf(fp, "# Frequency: %g [cpd]\n", freq);
      fprintf(fp, "# Field %zu: geocentric longitude (deg)\n", i++);
      fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
      fprintf(fp, "# Field %zu: J_r\n", i++);
      fprintf(fp, "# Field %zu: J_t\n", i++);
      fprintf(fp, "# Field %zu: J_p\n", i++);

      /* fill magfield workspace for this singular vector */
      fill_magfield(&v.vector, w);

      for (lon = -180.0; lon <= 180.0; lon += 1.0)
      /*for (lon = 0.0; lon <= 360.0; lon += 1.0)*/
        {
          double phi = lon * M_PI / 180.0;

          for (lat = -89.9; lat <= 89.9; lat += 1.0)
            {
              double theta = M_PI / 2.0 - lat * M_PI / 180.0;
              double J[3];

              magfield_eval_J(r, theta, phi, J, w);

              fprintf(fp, "%8.4f %8.4f %12.4e %12.4e %12.4e\n",
                      lon,
                      lat,
                      J[0],
                      J[1],
                      J[2]);
            }

          fprintf(fp, "\n");
        }

      fclose(fp);

      fprintf(stderr, "done\n");
    }

  return 0;
}

int
main(int argc, char *argv[])
{
  const double freq = 1.0; /* desired frequency in cpd */
  pca3d_data pcadata;
  pca3d_fft_data data;
  struct timeval tv0, tv1;
  gsl_matrix_complex *U;
  double alt = 110.0;
  size_t ifreq;            /* index of desired frequency */
  char buf[2048];
  
  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            alt = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-i U_data_file> [-a altitude (km)]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: frequency: %g [cpd]\n", freq);
  fprintf(stderr, "main: altitude:  %g [km]\n", alt);

  fprintf(stderr, "main: reading %s...", PCA3D_STAGE2B_FFT_DATA);
  gettimeofday(&tv0, NULL);
  data = pca3d_read_fft_data2(PCA3D_STAGE2B_FFT_DATA);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: allocating magfield workspace from %s...", PCA3D_STAGE1B_DATA);
  gettimeofday(&tv0, NULL);
  pcadata = pca3d_read_data(PCA3D_STAGE1B_DATA);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  ifreq = (size_t) (freq * data.window_size);
  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_U, ifreq);

  fprintf(stderr, "main: reading U matrix from %s...", buf);
  gettimeofday(&tv0, NULL);
  U = pca3d_read_matrix_complex(buf);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  print_modes(freq, alt + R_EARTH_KM, U, pcadata.w);

  gsl_matrix_complex_free(U);

  return 0;
}
