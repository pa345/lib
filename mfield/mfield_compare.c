/*
 * Compare mfield model with POMME or with another mfield model
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_pomme.h>

int
compare_pomme(msynth_workspace *msynth_p)
{
  pomme_workspace *pomme_p = pomme_alloc_default();
  double theta_deg, phi;
  const double R = 6371.2;
  const time_t t = 1277942400;
  const double tyr = get_year(t);
  size_t i = 0;

  pomme_set_deg(15, pomme_p);
  pomme_set_radius(R, pomme_p);

  printf("# Field %zu: longitude (degrees)\n", ++i);
  printf("# Field %zu: latitude (degrees)\n", ++i);
  printf("# Field %zu: X residual (nT)\n", ++i);
  printf("# Field %zu: Y residual (nT)\n", ++i);
  printf("# Field %zu: Z residual (nT)\n", ++i);
  printf("# Field %zu: F residual (nT)\n", ++i);

  for (theta_deg = 1.0e-2; theta_deg < 179.9; theta_deg += 1.0)
    {
      double theta = theta_deg * M_PI / 180.0;

      for (phi = 0.0; phi < 2.0 * M_PI; phi += 1.0 * M_PI / 180.0)
        {
          double B_pomme[4], B_model[4];
          size_t j;

          msynth_eval(tyr, R, theta, phi, B_model, msynth_p);
          pomme_calc_int(theta, phi, t, 0.0, B_pomme, pomme_p);

          for (j = 0; j < 4; ++j)
            B_pomme[j] *= 1.0e9;

          printf("%f %f %f %f %f %f\n",
                 wrap180(phi * 180.0 / M_PI),
                 90.0 - theta_deg,
                 B_model[0] - B_pomme[0],
                 B_model[1] - B_pomme[1],
                 B_model[2] - B_pomme[2],
                 B_model[3] - B_pomme[3]);
        }
    }

  pomme_free(pomme_p);

  return 0;
}

int
compare_models(const double t, msynth_workspace *msynth1, msynth_workspace *msynth2)
{
  double theta_deg, phi;
  const double R = 6371.2;
  size_t i = 0;

  printf("# Field %zu: longitude (degrees)\n", ++i);
  printf("# Field %zu: latitude (degrees)\n", ++i);
  printf("# Field %zu: X residual (nT)\n", ++i);
  printf("# Field %zu: Y residual (nT)\n", ++i);
  printf("# Field %zu: Z residual (nT)\n", ++i);
  printf("# Field %zu: F residual (nT)\n", ++i);
  printf("# Field %zu: dX/dt residual (nT)\n", ++i);
  printf("# Field %zu: dY/dt residual (nT)\n", ++i);
  printf("# Field %zu: dZ/dt residual (nT)\n", ++i);
  printf("# Field %zu: |dB/dt| residual (nT)\n", ++i);

  for (theta_deg = 1.0e-2; theta_deg < 179.9; theta_deg += 1.0)
    {
      double theta = theta_deg * M_PI / 180.0;

      for (phi = 0.0; phi < 2.0 * M_PI; phi += 1.0 * M_PI / 180.0)
        {
          double B1[4], B2[4], dBdt1[4], dBdt2[4];

          msynth_eval(t, R, theta, phi, B1, msynth1);
          msynth_eval_dBdt(t, R, theta, phi, dBdt1, msynth1);

          msynth_eval(t, R, theta, phi, B2, msynth2);
          msynth_eval_dBdt(t, R, theta, phi, dBdt2, msynth2);

          printf("%f %f %f %f %f %f %f %f %f %f\n",
                 wrap180(phi * 180.0 / M_PI),
                 90.0 - theta_deg,
                 B2[0] - B1[0],
                 B2[1] - B1[1],
                 B2[2] - B1[2],
                 B2[3] - B1[3],
                 dBdt2[0] - dBdt1[0],
                 dBdt2[1] - dBdt1[1],
                 dBdt2[2] - dBdt1[2],
                 dBdt2[3] - dBdt1[3]);
        }
    }

  return 0;
}

int
main(int argc, char *argv[])
{
  int c;
  msynth_workspace *msynth1 = NULL;
  msynth_workspace *msynth2 = NULL;
  double t = -1.0;

  while ((c = getopt(argc, argv, "c:a:b:it:")) != (-1))
    {
      switch (c)
        {
          case 'a':
          case 'c':
            fprintf(stderr, "main: reading %s...", optarg);
            msynth1 = msynth_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'b':
            fprintf(stderr, "main: reading %s...", optarg);
            msynth2 = msynth_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'i':
            fprintf(stderr, "main: reading IGRF coefficients...");
            msynth2 = msynth_igrf_read(MSYNTH_IGRF_FILE);
            fprintf(stderr, "done\n");
            break;

          case 't':
            t = atof(optarg);
            break;

          default:
            break;
        }
    }

  if (!msynth1)
    {
      fprintf(stderr, "Usage: %s [-c coef_ascii_file] [-a coef_ascii_file1 -b coef_ascii_file2] [-i] [-t decimal_year]\n", argv[0]);
      exit(1);
    }

  if (t < 0.0)
    t = msynth1->epochs[msynth1->n_epochs - 1];

  if (!msynth2)
    compare_pomme(msynth1);
  else
    compare_models(t, msynth1, msynth2);

  if (msynth1)
    msynth_free(msynth1);
  if (msynth2)
    msynth_free(msynth2);

  return 0;
}
