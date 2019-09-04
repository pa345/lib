/*
 * mfield_dBdt
 *
 * Compute dB/dt for various observatory locations
 * with DMSP model
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <gsl/gsl_math.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>

#include "mfield.h"

int
main(int argc, char *argv[])
{
  char *coeffile = NULL;
  mfield_workspace *w;
  double r, theta, phi, t;
  double dBdt[4];
  int c;

  while ((c = getopt(argc, argv, "c:")) != (-1))
    {
      switch (c)
        {
          case 'c':
            coeffile = optarg;
            break;

          default:
            fprintf(stderr, "usage: %s <-c coef_input_file>\n",
                    argv[0]);
            exit(1);
            break;
        }
    }

  if (!coeffile)
    {
      fprintf(stderr, "usage: %s <-c coef_input_file>\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: reading coefficients from %s...", coeffile);
  w = mfield_read(coeffile);
  if (!w)
    exit(1);
  fprintf(stderr, "done\n");

#if 1
  /* KOU */
  r = R_EARTH_KM + 10.0 * 1.0e-3;
  theta = M_PI / 2.0 - 5.21 * M_PI / 180.0;
  phi = 307.269 * M_PI / 180.0;
#endif

#if 0
  /* TAM */
  r = R_EARTH_KM + 1.373;
  theta = M_PI / 2.0 - 22.79 * M_PI / 180.0;
  phi = 5.53 * M_PI / 180.0;
#endif

#if 0
  /* MBO */
  r = R_EARTH_KM + 7.0 * 1.0e-3;
  theta = M_PI / 2.0 - 14.38 * M_PI / 180.0;
  phi = 343.03 * M_PI / 180.0;
#endif

  t = w->epoch;

  mfield_eval_dBdt(t, r, theta, phi, dBdt, w);
  printf("%f %f %f %f\n",
         t,
         dBdt[0],
         dBdt[1],
         dBdt[2]);

  mfield_free(w);

  return 0;
} /* main() */
