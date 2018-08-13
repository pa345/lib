/*
 * print2.c
 *
 * Print lat/lon maps from stage1_knm and stage1_gnm files
 *
 * ./print2 <-k knm_input_file> <-g gnm_input_file> [-t time_idx]
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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <common/common.h>

#include "green.h"

#include "io.h"

/*
print_data()
  Print Bx/By/Bz grid for a fixed time

Inputs: filename - output file
        kg       - vector of knm or gnm coefficients
        type     - 'k' for knm, 'g' for gnm
*/

int
print_data(const char *filename, const gsl_vector *kg, const char type)
{
  int s = 0;
  size_t i;
  FILE *fp;
  double lon, lat;
  green_workspace *green_p;
  gsl_vector *X, *Y, *Z;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  green_p = green_alloc(60, 30, R_EARTH_KM);
  X = gsl_vector_alloc(green_p->nnm);
  Y = gsl_vector_alloc(green_p->nnm);
  Z = gsl_vector_alloc(green_p->nnm);

  assert(kg->size == green_p->nnm);

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: B_x (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z (nT)\n", i++);

  for (lon = 0.0; lon <= 360.0; lon += 5.0)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -89.9; lat <= 89.9; lat += 2.0)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double B[3];

          if (type == 'k')
            green_calc_ext(R_EARTH_KM, theta, phi, X->data, Y->data, Z->data, green_p);
          else
            green_calc_int(R_EARTH_KM, theta, phi, X->data, Y->data, Z->data, green_p);

          gsl_blas_ddot(kg, X, &B[0]);
          gsl_blas_ddot(kg, Y, &B[1]);
          gsl_blas_ddot(kg, Z, &B[2]);

          fprintf(fp, "%8.4f %8.4f %8.2f %8.2f %8.2f\n",
                  lon,
                  lat,
                  B[0],
                  B[1],
                  B[2]);
        }

      fprintf(fp, "\n");
    }

  green_free(green_p);

  fclose(fp);

  return s;
}

int
main(int argc, char *argv[])
{
  gsl_matrix *G = NULL;
  gsl_matrix *K = NULL;
  int time_idx = 0;
  char *outfile = "data.txt";
  gsl_vector_view v;
  char type = 'k';

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "k:g:o:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'k':
            fprintf(stderr, "main: reading knm coefficients from %s...", optarg);
            K = pca_read_matrix(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'g':
            fprintf(stderr, "main: reading gnm coefficients from %s...", optarg);
            G = pca_read_matrix(optarg);
            fprintf(stderr, "done\n");
            break;

          case 't':
            time_idx = atol(optarg);
            break;

          case 'o':
            outfile = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s <-k knm_file> <-g gnm_file> [-t time_idx] [-o output_file]\n", argv[0]);
            break;
        }
    }

  if (!K && !G)
    {
      fprintf(stderr, "Usage: %s <-k knm_file> <-g gnm_file> [-t time_idx] [-o output_file]\n", argv[0]);
      exit(1);
    }

  if (K)
    {
      v = gsl_matrix_column(K, time_idx);
      type = 'k';
    }
  else
    {
      v = gsl_matrix_column(G, time_idx);
      type = 'g';
    }

  fprintf(stderr, "main: writing data to %s...", outfile);
  print_data(outfile, &v.vector, type);
  fprintf(stderr, "done\n");

  if (K)
    gsl_matrix_free(K);

  if (G)
    gsl_matrix_free(G);

  return 0;
}
