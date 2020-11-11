/*
 * print_spectrum1b.c
 *
 * 1. Read SH coefficients computed from stage1b
 * 2. Print spectrum at given radius and time
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_blas.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_bsearch.h>

#include <magfield/magfield.h>
#include <magfield/magfield_eval.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

/*
main_proc()
  Compute SH decomposition of each radial shell of TIEGCM 3D current
values; store SH coefficients to disk

Inputs: filename - output file
        tidx     - time index
        alt      - altitude (km)
        w        - magfield workspace
        eval_p   - magfield eval workspace

Return: success/error
*/

int
print_spectrum(const char * filename, const size_t tidx, const double alt,
               magfield_workspace * w, magfield_eval_workspace * eval_p)
{
  int s = 0;
  const double r = R_EARTH_M + alt * 1.0e3;
  const size_t lmax = eval_p->lmax;
  FILE *fp = fopen(filename, "w");
  char buf[2048];
  size_t l;

  sprintf(buf, "%s_%03zu.dat", PCA3D_STAGE1B_SH_PREFIX, tidx + 1);
  fprintf(stderr, "print_spectrum: reading %s...", buf);
  magfield_read_SH(buf, w);
  fprintf(stderr, "done\n");

  fprintf(stderr, "print_spectrum: initializing magfield_eval...");
  magfield_eval_init(w->qtcoeff, w->qcoeff, w->pcoeff, eval_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "print_spectrum: writing %s...", filename);

  l = 1;
  fprintf(fp, "# SH coef file: %s\n", buf);
  fprintf(fp, "# Radius: %.2f [km] (altitude: %.2f [km])\n", r * 1.0e-3, alt);
  fprintf(fp, "# Field %zu: spherical harmonic degree\n", l++);
  fprintf(fp, "# Field %zu: power (nT^2)\n", l++);

  for (l = 1; l <= lmax; ++l)
    {
      double Rl = magfield_eval_spectrum(r, l, eval_p);
      fprintf(fp, "%4zu %.12e\n", l, Rl * 1.0e18);
    }

  fclose(fp);

  fprintf(stderr, "done\n");

  return s;
}

int
main(int argc, char *argv[])
{
  char *output_file = "spectrum.txt";
  struct timeval tv0, tv1;
  magfield_eval_workspace * magfield_eval_p = NULL;
  size_t tidx = 0;
  double alt = 110.0;
  pca3d_data pdata;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "output_file", required_argument, NULL, 'o' },
          { "altitude", required_argument, NULL, 'a' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:o:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            alt = atoi(optarg);
            break;

          case 'o':
            output_file = optarg;
            break;

          case 't':
            tidx = (size_t) atoi(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s [-a altitude (km)] [-o output_file] [-t tidx]\n", argv[0]);
            break;
        }
    }

  fprintf(stderr, "main: reading magfield parameters %s...", PCA3D_STAGE1B_DATA);
  pdata = pca3d_read_data(PCA3D_STAGE1B_DATA);
  fprintf(stderr, "done (%zu timestamps)\n", pdata.nt);

  fprintf(stderr, "main: lmax = %zu\n", pdata.w->lmax);
  fprintf(stderr, "main: mmax = %zu\n", pdata.w->mmax);

  fprintf(stderr, "main: allocating magfield eval workspace...");
  magfield_eval_p = magfield_eval_alloc(&(pdata.w->params));
  fprintf(stderr, "done\n");

  print_spectrum(output_file, tidx, alt, pdata.w, magfield_eval_p);

  return 0;
}
