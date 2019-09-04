/*
 * plot_maps.c
 *
 * Print lat/lon maps of magnetic field and current flow due to
 * first few PCs in time domain
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>

#include <mainlib/ml_common.h>

#include "pca.h"

int
main(int argc, char *argv[])
{
  pca_workspace *pca_workspace_p;
  char buf[2048];
  size_t i;
  double r = R_EARTH_KM + 0.0; /* radius for magnetic field maps */
  size_t ut = 12;              /* UT hour for maps */
  pca_source_t source_type = PCA_SRC_EXTERNAL;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:t:g", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            r = R_EARTH_KM + atof(optarg);
            break;

          case 't':
            ut = (size_t) atoi(optarg);
            break;

          case 'g':
            source_type = PCA_SRC_INDUCED;
            break;

          default:
            fprintf(stderr, "Usage: %s [-a altitude] [-t UT] [-g]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: radius = %g [km]\n", r);

  pca_workspace_p = pca_alloc(source_type);
  pca_set_UT(ut, pca_workspace_p);
  
  sprintf(buf, "maps/map_mean_%02zuUT.dat", ut);
  fprintf(stderr, "main: printing mean map to %s...", buf);
  pca_print_mean_map(buf, r, pca_workspace_p);
  fprintf(stderr, "done\n");

  for (i = 0; i < 40; ++i)
    {
      sprintf(buf, "maps/map_%02zuUT_%02zu.dat", ut, i + 1);

      fprintf(stderr, "main: printing principle component %zu map to %s...",
              i + 1, buf);

      pca_print_pc_map(buf, r, i, pca_workspace_p);

      fprintf(stderr, "done\n");
    }

  pca_free(pca_workspace_p);

  return 0;
}
