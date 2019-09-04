/*
 * plot_modes.c
 *
 * Read eigenvector file (left singular vector file) containing modes U
 * and print them in format suitable for plotting
 *
 * ./plot_modes [-i U_data_file]
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

#include <mainlib/ml_common.h>
#include <mainlib/ml_bsearch.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

int
print_modes(const size_t ir, const gsl_matrix_complex *U, const tiegcm3d_fft_data *data)
{
  const size_t N = data->nr * data->nlat * data->nlon;
  const size_t T = U->size2;
  size_t t, ilat, ilon;

  /* loop over modes t \in [0, T-1] */
  for (t = 0; t < T; ++t)
    {
      gsl_vector_complex_const_view v = gsl_matrix_complex_const_column(U, t);
      FILE *fp;
      char buf[2048];
      size_t i = 1;

      sprintf(buf, "plots/U_%02zu.txt", t + 1);
      fprintf(stderr, "print_modes: writing %s...", buf);

      fp = fopen(buf, "w");

      fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
      fprintf(fp, "# Field %zu: geocentric longitude (deg)\n", i++);
      fprintf(fp, "# Field %zu: Re [ U_r ]\n", i++);
      fprintf(fp, "# Field %zu: Im [ U_r ]\n", i++);
      fprintf(fp, "# Field %zu: Re [ U_t ]\n", i++);
      fprintf(fp, "# Field %zu: Im [ U_t ]\n", i++);
      fprintf(fp, "# Field %zu: Re [ U_p ]\n", i++);
      fprintf(fp, "# Field %zu: Im [ U_p ]\n", i++);

      for (ilon = 0; ilon < data->nlon; ++ilon)
        {
          for (ilat = 0; ilat < data->nlat; ++ilat)
            {
              size_t idx = CIDX3(ir, data->nr, ilat, data->nlat, ilon, data->nlon);
              gsl_complex Ur = gsl_vector_complex_get(&v.vector, idx);
              gsl_complex Ut = gsl_vector_complex_get(&v.vector, idx + N);
              gsl_complex Up = gsl_vector_complex_get(&v.vector, idx + 2*N);

              fprintf(fp, "%8.4f %8.4f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                      data->glon[ilon],
                      data->glat[ilat],
                      GSL_REAL(Ur),
                      GSL_IMAG(Ur),
                      GSL_REAL(Ut),
                      GSL_IMAG(Ut),
                      GSL_REAL(Up),
                      GSL_IMAG(Up));
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
  tiegcm3d_fft_data data;
  struct timeval tv0, tv1;
  char *infile = NULL;
  gsl_matrix_complex *U;
  double alt = 110.0;
  size_t r_idx;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            alt = atof(optarg);
            break;

          case 'i':
            infile = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s <-i U_data_file> [-a altitude (km)]\n", argv[0]);
            exit(1);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i U_data_file> [-a altitude (km)]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: reading TIEGCM metadata from %s...", PCA3D_STAGE2A_FFT_DATA_LIGHT);
  gettimeofday(&tv0, NULL);
  data = pca3d_read_fft_data(PCA3D_STAGE2A_FFT_DATA_LIGHT);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: reading U matrix from %s...", infile);
  gettimeofday(&tv0, NULL);
  U = pca3d_read_matrix_complex(infile);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  r_idx = bsearch_double(data.r, alt + R_EARTH_KM, 0, data.nr - 1);
  fprintf(stderr, "main: r_idx = %zu (%.2f [km])\n", r_idx, data.r[r_idx] - R_EARTH_KM);

  print_modes(r_idx, U, &data);

  gsl_matrix_complex_free(U);

  return 0;
}
