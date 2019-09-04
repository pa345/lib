/*
 * mfield_plot.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_msynth.h>

#include "mfield.h"

int
plot_field(msynth_workspace *w)
{
#if 1
  const double r = R_EARTH_KM;
#else
  const double r = 3485.0;
#endif
  const double t = 2010.0;
  double theta, phi;
  double dphi = 1.0 * M_PI / 180.0;
  double dtheta = 1.0 * M_PI / 180.0;

  for (phi = -M_PI; phi < M_PI; phi += dphi)
    {
      for (theta = 0.01; theta < M_PI; theta += dtheta)
        {
          double B[4], dBdt[4];

          msynth_eval(t, r, theta, phi, B, w);
          msynth_eval_dBdt(t, r, theta, phi, dBdt, w);

          printf("%f %f %e %e %e %e\n",
                 phi * 180.0 / M_PI,
                 90.0 - theta * 180.0 / M_PI,
                 B[3],
                 dBdt[3],
                 B[2],
                 dBdt[2]);
        }
      printf("\n");
    }

  return GSL_SUCCESS;
}

int
main(int argc, char *argv[])
{
  int c;
  char *coeffile = NULL;
  int wmm = 0;
  msynth_workspace *msynth_p = NULL;

  while ((c = getopt(argc, argv, "c:w")) != (-1))
    {
      switch (c)
        {
          case 'c':
            coeffile = optarg;
            break;

          case 'w':
            wmm = 1;
            break;

          default:
            break;
        }
    }

  if (coeffile)
    {
      fprintf(stderr, "main: loading coefficients from %s...", coeffile);
      msynth_p = msynth_read(coeffile);
      fprintf(stderr, "done\n");
    }
  else if (wmm)
    {
      fprintf(stderr, "main: loading WMM coefficients...");
      msynth_p = msynth_wmm_read(MSYNTH_WMM_FILE);
      /*msynth_p = msynth_wmm_read(MSYNTH_WMM_BGS_FILE);*/
      fprintf(stderr, "done\n");
    }
  else
    {
      fprintf(stderr, "usage: %s [-c ascii_coef_file] [-w]\n", argv[0]);
      exit(1);
    }

  plot_field(msynth_p);

  msynth_free(msynth_p);

  return 0;
}
