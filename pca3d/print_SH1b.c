/*
 * print_SH.c
 *
 * 1. Read previously calculated SH coefficients for each timestep (stage1b)
 * 2. Print time series of q_l^m(r,t), etc
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_test.h>

#include <mainlib/ml_common.h>

#include <magfield/magfield.h>

#include "io.h"
#include "pca3d.h"
#include "window.h"

#include "common.c"

/*
main_proc()

Inputs: data         - pca3d data
        prefix       - input file prefix

Return: success/error
*/

int
main_proc(const pca3d_data * data, const char * prefix, const size_t ir, const size_t l, const int m)
{
  int s = 0;
  const size_t nt = data->nt;
  char buf[2048];
  size_t i;
  struct timeval tv0, tv1;
  magfield_workspace *w = data->w;
  const size_t lmidx = magfield_lmidx(l, m, w->mmax);
  const size_t idx = MAG_COEFIDX(ir, lmidx, w);

  i = 1;
  printf("# SH degree: %zu\n", l);
  printf("# SH order:  %d\n", m);
  printf("# Field %zu: timestamp\n", i++);
  printf("# Field %zu: Re [ qt_l^m(r) ] (nT)\n", i++);
  printf("# Field %zu: Im [ qt_l^m(r) ] (nT)\n", i++);
  printf("# Field %zu: Re [ p_l^m(r) ] (nT)\n", i++);
  printf("# Field %zu: Im [ p_l^m(r) ] (nT)\n", i++);
  printf("# Field %zu: Re [ q_l^m(r) ] (nT)\n", i++);
  printf("# Field %zu: Im [ q_l^m(r) ] (nT)\n", i++);

  for (i = 0; i < nt; ++i)
    {
      gsl_complex qt, p, q;

      sprintf(buf, "%s_%03zu.dat", prefix, i + 1);
      fprintf(stderr, "main_proc: reading %s...", buf);
      gettimeofday(&tv0, NULL);

      /* read previously calculated SH coefficients for this timestep */
      magfield_read_SH(buf, w);

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      qt = gsl_vector_complex_get(w->qtcoeff, idx);
      p = gsl_vector_complex_get(w->pcoeff, idx);
      q = gsl_vector_complex_get(w->qcoeff, idx);

      printf("%ld %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
             data->t[i],
             GSL_REAL(qt) * 1.0e9,
             GSL_IMAG(qt) * 1.0e9,
             GSL_REAL(p) * 1.0e9,
             GSL_IMAG(p) * 1.0e9,
             GSL_REAL(q) * 1.0e9,
             GSL_IMAG(q) * 1.0e9);
    }

  magfield_free(w);

  return s;
}

int
main(int argc, char *argv[])
{
  char *input_prefix = PCA3D_STAGE1B_SH_PREFIX;
  pca3d_data data;
  size_t l = 1;
  int m = 0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "l:m:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'l':
            l = (size_t) atoi(optarg);
            break;

          case 'm':
            m = atoi(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s [-l sh_degree] [-m sh_order]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: l = %zu\n", l);
  fprintf(stderr, "main: m = %d\n", m);

  fprintf(stderr, "main: allocating magfield workspace from %s...", PCA3D_STAGE1B_DATA);
  data = pca3d_read_data(PCA3D_STAGE1B_DATA);
  fprintf(stderr, "done\n");

  main_proc(&data, input_prefix, 0, l, m);

  return 0;
}
