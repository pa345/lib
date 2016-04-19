/*
 * stage2b.c
 *
 * 1. Read in spherical harmonic time series knm from stage1 matrix file (nnm-by-nt)
 * 2. Compute SVD of knm matrix
 * 3. Write singular values and singular vectors to disk
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>

#include <fftw3.h>
#include <lapacke/lapacke.h>

#include "common.h"
#include "lapack_wrapper.h"

#include "io.h"
#include "pca.h"

int
main(int argc, char *argv[])
{
  size_t nnm;                             /* number of spherical harmonic time series */
  size_t nt;                              /* number of time steps */
  char *infile = PCA_STAGE1_KNM;
  gsl_matrix *K;
  struct timeval tv0, tv1;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s <-i stage1_matrix_file>\n", argv[0]);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i stage1_matrix_file>\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  K = pca_read_matrix(infile);
  fprintf(stderr, "done (%zu-by-%zu matrix)\n", K->size1, K->size2);

  nnm = K->size1;
  nt = K->size2;

  /* compute 1/sqrt(nt) K for SVD computation */
  gsl_matrix_scale(K, 1.0 / sqrt((double) nt));

  {
    const char *sval_txt_file = "sval_time.txt";
    const char *sval_file = PCA_STAGE2B_SVAL;
    const char *U_file = PCA_STAGE2B_U;
    const char *V_file = PCA_STAGE2B_V;
    gsl_vector *S = gsl_vector_alloc(GSL_MIN(nnm, nt));
    gsl_matrix *U = gsl_matrix_alloc(nnm, nnm);
    gsl_matrix *V = gsl_matrix_alloc(nt, nt);
    int status;
    FILE *fp;

    fprintf(stderr, "main: performing SVD of K matrix (%zu-by-%zu)...", nnm, nt);
    gettimeofday(&tv0, NULL);
    status = lapack_svd(K, S, U, V);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds, status = %d)\n",
            time_diff(tv0, tv1), status);

    fprintf(stderr, "main: writing singular values in text format to %s...",
            sval_txt_file);
    fp = fopen(sval_txt_file, "w");
    gsl_vector_fprintf(fp, S, "%.12e");
    fclose(fp);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing singular values in binary format to %s...",
            sval_file);
    pca_write_vector(sval_file, S);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing left singular vectors in binary format to %s...",
            U_file);
    pca_write_matrix(U_file, U);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing right singular vectors in binary format to %s...",
            V_file);
    pca_write_matrix(V_file, V);
    fprintf(stderr, "done\n");

    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(S);
  }

  gsl_matrix_free(K);

  return 0;
}
