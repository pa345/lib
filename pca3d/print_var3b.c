/*
 * print_var3b.c
 *
 * Read singular value file and print the cumulative variance curve
 *
 * ./print_var3b [-i S_data_file]
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
#include <gsl/gsl_blas.h>

#include <common/common.h>

#include "io.h"
#include "pca3d.h"

int
print_variance(const char *filename, gsl_vector *eval)
{
  int s = 0;
  const size_t n = eval->size;
  const double lambda_max = gsl_vector_max(eval);
  FILE *fp = fopen(filename, "w");
  size_t i;
  double sum_all, cumsum = 0.0;

  /* normalize eigenvalues and compute sum */
  gsl_vector_scale(eval, 1.0 / lambda_max);
  sum_all = gsl_blas_dasum(eval);

  i = 1;
  fprintf(fp, "# Field %zu: eigenvalue number i\n", i++);
  fprintf(fp, "# Field %zu: normalized eigenvalue\n", i++);
  fprintf(fp, "# Field %zu: cumulative variance: ( sum_{j=1}^i lambda_j ) / ( sum_j lambda_j )\n", i++);

  for (i = 0; i < n; ++i)
    {
      double lambda = gsl_vector_get(eval, i);
      cumsum += lambda;
      fprintf(fp, "%zu %.12e %.12e\n", i + 1, lambda, cumsum / sum_all);
    }

  fclose(fp);

  return s;
}

int
main(int argc, char *argv[])
{
  char * outfile = "plots/variance.txt";
  gsl_vector * S = NULL;
  
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
            S = pca3d_read_vector(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-i S_data_file>\n", argv[0]);
            exit(1);
            break;
        }
    }

  if (S == NULL)
    {
      fprintf(stderr, "Usage: %s <-i S_data_file>\n", argv[0]);
      exit(1);
    }

  /* compute eigenvalues from singular values */
  gsl_vector_mul(S, S);

  fprintf(stderr, "main: writing %s...", outfile);
  print_variance(outfile, S);
  fprintf(stderr, "done\n");

  gsl_vector_free(S);

  return 0;
}
