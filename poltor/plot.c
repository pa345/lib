/*
 * plot.c
 *
 * Plot the toroidal current output of invert_main
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

#include <flow/flow.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_common.h>

#include "poltor.h"

int
main_K(const double x[], double J[], void *params)
{
  int s;
  poltor_workspace *w = (poltor_workspace *) params;
  double theta = M_PI / 2.0 - x[1] * M_PI / 180.0;
  double phi = x[0] * M_PI / 180.0;
  double r = w->R + 110.0;
  double Jtmp[3] = { 0.0, 0.0, 0.0 };

  s = poltor_eval_K_tor(r, theta, phi, Jtmp, w);

  J[0] = Jtmp[2];  /* J_x = J_east = J_phi */
  J[1] = -Jtmp[1]; /* J_y = J_north = -J_theta */

  return s;
}

int
main_J_shell(const double x[], double J[], void *params)
{
  int s;
  poltor_workspace *w = (poltor_workspace *) params;
  const double r = 350.0;
  double theta = M_PI / 2.0 - x[1] * M_PI / 180.0;
  double phi = x[0] * M_PI / 180.0;
  double Jtmp[3] = { 0.0, 0.0, 0.0 };

  s = poltor_eval_J_shell(r, theta, phi, Jtmp, w);

  J[0] = Jtmp[2];  /* J_x = J_east = J_phi */
  J[1] = -Jtmp[1]; /* J_y = J_north = -J_theta */

  return s;
}

/* compute height-integrated current density in A/m */
int
main_J_shell_hi(const double x[], double J[], void *params)
{
  int s;
  poltor_workspace *w = (poltor_workspace *) params;
  const double rmin = 250.0;
  const double rmax = 500.0;
  const size_t nr = 500;
  const double dr = (rmax - rmin) / (nr - 1.0);
  double theta = M_PI / 2.0 - x[1] * M_PI / 180.0;
  double phi = x[0] * M_PI / 180.0;
  double Jtmp[3] = { 0.0, 0.0, 0.0 };
  size_t i;

  J[0] = 0.0;
  J[1] = 0.0;

  for (i = 0; i < nr; ++i)
    {
      double r = rmin + i * dr;

      s = poltor_eval_J_shell(r, theta, phi, Jtmp, w);

      J[0] += Jtmp[2];  /* J_x = J_east = J_phi */
      J[1] += -Jtmp[1]; /* J_y = J_north = -J_theta */
    }

  J[0] *= dr;
  J[1] *= dr;

  return s;
}

int
plot_flow(const double ut, poltor_workspace *w)
{
  int s = 0;
  flow_workspace *flow_p;
#if 1
  const double lon_min = -180.0;
  const double lon_max = 180.0;
  double scale = 6.0;
#else
  const double lon_min = 120.0;
  const double lon_max = 135.0;
  double scale = 1.0;
#endif
  const double lat_min = -60.0;
  const double lat_max = 60.0;

  flow_p = flow_alloc(lon_min, lon_max, lat_min, lat_max);

  flow_set_vfield(&main_K, w, flow_p);

  {
    char *quiver_file = "quiver.dat";

    fprintf(stderr, "plot_flow: writing quiver to %s...", quiver_file);
    flow_calc_quiver(quiver_file, 100, 100, &scale, flow_p);
    fprintf(stderr, "done (final scale = %f)\n", scale);
  }

  {
    fprintf(stderr, "plot_flow: plotting current flow for matlab...");
    flow_print_matlab(100, 100, "X", "Y", "U", "V", flow_p);
    fprintf(stderr, "done\n");
  }

#if 0
  {
    double lt = 12.0;
    double phi0, lat0;
    char *pathline_file = "flowpath.dat";
    FILE *fp;
    size_t i;

    fp = fopen(pathline_file, "w");
    if (!fp)
      {
        fprintf(stderr, "plot_flow: unable to open %s: %s\n",
                pathline_file, strerror(errno));
        return -1;
      }

    i = 1;
    fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
    fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);

    fprintf(stderr, "plot_flow: writing pathlines to %s...", pathline_file);

    /* find longitude corresponding to desired LT/UT */
    phi0 = (lt - ut) * 15.0;
    /*for (phi0 = lon_min; phi0 < lon_max; phi0 += 25.0)*/
      {
        for (lat0 = -60.0; lat0 < 60.0; lat0 += 5.0)
          {
            /* trace line forward */
            flow_calc_pathline(phi0, lat0, 1, flow_p);
            flow_print_pathline(fp, flow_p);

            /* trace line backward */
            flow_calc_pathline(phi0, lat0, -1, flow_p);
            flow_print_pathline(fp, flow_p);
          }
      }

    fprintf(stderr, "done\n");

    fclose(fp);
  }
#endif

  /* plot current stream function grid */
  {
    const double b = w->R + 110.0;
    const double d = w->R + 350.0;
    const double lat_min2 = -80.0;
    const double lat_max2 = 80.0;
#if 1
    const double lon_min2 = lon_min;
    const double lon_max2 = lon_max;
#else
    const double lon_min2 = 45.0;
    const double lon_max2 = 105.0;
#endif
    const double phi_min = lon_min2 * M_PI / 180.0;
    const double phi_max = lon_max2 * M_PI / 180.0;
    const double theta_min = M_PI / 2.0 - lat_max2 * M_PI / 180.0;
    const double theta_max = M_PI / 2.0 - lat_min2 * M_PI / 180.0;
    double theta, phi;
    char *chi_file = "chi.dat";
    FILE *fp;
    size_t i;

    fp = fopen(chi_file, "w");
    if (!fp)
      {
        fprintf(stderr, "plot_flow: unable to open %s: %s\n",
                chi_file, strerror(errno));
        return -1;
      }

    i = 1;
    fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
    fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
    fprintf(fp, "# Field %zu: chi (internal) (kA)\n", i++);
    fprintf(fp, "# Field %zu: chi (shell) (kA)\n", i++);

    fprintf(stderr, "plot_flow: writing current stream function to %s...", chi_file);

    for (phi = phi_min; phi < phi_max; phi += 2.0 * M_PI / 180.0)
      {
        for (theta = theta_min; theta < theta_max; theta += 2.0 * M_PI / 180.0)
          {
            double chi_int, chi_sh; /* current stream function */

            poltor_eval_chi_int(b, theta, phi, &chi_int, w);
            poltor_eval_chi_sh(d, theta, phi, &chi_sh, w);

            fprintf(fp, "%f %f %f %f\n",
                    phi * 180.0 / M_PI,
                    90.0 - theta * 180.0 / M_PI,
                    chi_int,
                    chi_sh);
          }
        fprintf(fp, "\n");
      }

    fprintf(stderr, "done\n");

    fclose(fp);
  }

  flow_free(flow_p);

  return s;
} /* plot_flow() */

int
plot_flow_shell(poltor_workspace *w)
{
  int s = 0;
  flow_workspace *flow_p;
  const double lon_min = -180.0;
  const double lon_max = 180.0;
  const double lat_min = -60.0;
  const double lat_max = 60.0;

  flow_p = flow_alloc(lon_min, lon_max, lat_min, lat_max);

  flow_set_vfield(&main_J_shell, w, flow_p);

  {
    char *quiver_file = "quiver_shell.dat";
    double scale = 1.0;

    fprintf(stderr, "plot_flow_shell: writing quiver to %s...", quiver_file);
    flow_calc_quiver(quiver_file, 50, 50, &scale, flow_p);
    fprintf(stderr, "done (final scale = %f)\n", scale);
  }

#if 0
  flow_set_vfield(&main_J_shell_hi, w, flow_p);

  {
    char *quiver_file = "quiver_shell_hi.dat";
    double scale = 1.0;

    fprintf(stderr, "plot_flow_shell: writing height-integrated quiver to %s...", quiver_file);
    flow_calc_quiver(quiver_file, 50, 50, &scale, flow_p);
    fprintf(stderr, "done (final scale = %f)\n", scale);
  }
#endif

  flow_free(flow_p);

  return s;
} /* plot_flow_shell() */

int
plot_B(const char *filename, const double r, poltor_workspace *w)
{
  int s = 0;
  const double phi_min = -M_PI;
  const double phi_max = M_PI;
  const double dphi = 1.0 * M_PI / 180.0;
  const double lat_min = -89.9;
  const double lat_max = 89.9;
  const double dlat = 0.5;
  double phi, lat;
  double B_int[4], B_ext[4], B_sh[4], B_tor[4];
  FILE *fp;
  size_t i;
  const time_t unix_time = 1427846400; /* Apr 1 2015 00:00:00 UTC */
  const double t = satdata_timet2epoch(unix_time);

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "plot_B: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal model X (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal model Y (nT)\n", i++);
  fprintf(fp, "# Field %zu: internal poloidal model Z (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal model X (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal model Y (nT)\n", i++);
  fprintf(fp, "# Field %zu: shell poloidal model Z (nT)\n", i++);

  for (phi = phi_min; phi <= phi_max; phi += dphi)
    {
      for (lat = lat_min; lat <= lat_max; lat += dlat)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;

          poltor_eval_B_all(t, r, theta, phi, B_int, B_ext, B_sh, B_tor, w);

          fprintf(fp, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",
                  phi * 180.0 / M_PI,
                  lat,
                  B_int[0],
                  B_int[1],
                  B_int[2],
                  B_sh[0],
                  B_sh[1],
                  B_sh[2]);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
} /* plot_B() */

int
print_coefficients(poltor_workspace *w)
{
  const size_t nmax_int = GSL_MIN(3, w->nmax_int);
  const size_t mmax_int = GSL_MIN(3, w->mmax_int);
  const size_t nmax_ext = GSL_MIN(2, w->nmax_ext);
  const size_t mmax_ext = GSL_MIN(2, w->mmax_ext);
  const size_t nmax_tor = GSL_MIN(3, w->nmax_tor);
  const size_t mmax_tor = GSL_MIN(3, w->mmax_tor);
  size_t n;

  /* print internal poloidal coefficients */
  fprintf(stderr, "Internal poloidal coefficients:\n");
  for (n = 1; n <= nmax_int; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax_int);
      int m;

      /*
       * only need to print positive m coefficients since
       * c_{n,-m} = c_{nm}
       */
      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_complex coef = poltor_get(cidx, w);
          double gnm = GSL_REAL(coef);

          fprintf(stderr, "g(%2zu,%2d) = %12g [nT]\n",
                  n, m, gnm);
        }
    }

  /* print external poloidal coefficients */
  fprintf(stderr, "External poloidal coefficients:\n");
  for (n = 1; n <= nmax_ext; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax_ext);
      int m;

      /*
       * only need to print positive m coefficients since
       * c_{n,-m} = c_{nm}
       */
      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          gsl_complex coef = poltor_get(cidx, w);
          double knm = GSL_REAL(coef);

          fprintf(stderr, "k(%2zu,%2d) = %12g [nT]\n",
                  n,
                  m,
                  knm);
        }
    }

  /* print toroidal coefficients */
  fprintf(stderr, "Shell toroidal coefficients:\n");
  for (n = 1; n <= nmax_tor; ++n)
    {
      int ni = (int) GSL_MIN(n, mmax_tor);
      int m;

      /*
       * only need to print positive m coefficients since
       * c_{n,-m} = c_{nm}
       */
      for (m = 0; m <= ni; ++m)
        {
          size_t cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
          gsl_complex coef = poltor_get(cidx, w);
          double phinm = GSL_REAL(coef);

          fprintf(stderr, "phi(%2zu,%2d) = %12g [nT]\n",
                  n,
                  m,
                  phinm);
        }
    }

  return 0;
} /* print_coefficients() */

int
print_slice(const char *filename, const double r, const double phi, poltor_workspace *w)
{
  FILE *fp;
  double theta;
  size_t i;
  const time_t unix_time = 1427846400; /* Apr 1 2015 00:00:00 UTC */
  const double t = satdata_timet2epoch(unix_time);

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_slice: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: toroidal K_theta on r = b shell (kA/km)\n", i++);
  fprintf(fp, "# Field %zu: toroidal K_phi on r = b shell (kA/km)\n", i++);
  fprintf(fp, "# Field %zu: B_x (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z (nT)\n", i++);

  for (theta = 0.01; theta < M_PI; theta += 1.0*M_PI/180.0)
    {
      double K[3] = { 0.0, 0.0, 0.0 };
      double B[4];

      poltor_eval_K_tor(0.0, theta, phi, K, w);
      poltor_eval_B(t, r, theta, phi, B, w);

      fprintf(fp, "%f %e %e %e %e %e\n",
              90.0 - theta * 180.0 / M_PI,
              K[1],
              K[2],
              B[0],
              B[1],
              B[2]);
    }

  fclose(fp);

  return 0;
} /* print_slice() */

int
main(int argc, char *argv[])
{
  int c;
  double universal_time = 11.0; /* UT in hours for data selection */
  char *spectrum_file = "poltor.s";
  poltor_workspace *poltor_p = NULL;

  while ((c = getopt(argc, argv, "t:i:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            poltor_p = poltor_read(optarg);
            break;

          case 't':
            universal_time = atof(optarg);
            break;

          default:
            break;
        }
    }

  if (!poltor_p)
    {
      fprintf(stderr, "Usage: %s [-t universal_time] [-i binary_coef_file]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: universal time = %.1f\n", universal_time);

  print_coefficients(poltor_p);

  fprintf(stderr, "main: printing spectrum to %s...", spectrum_file);
  poltor_print_spectrum(spectrum_file, poltor_p->c, poltor_p);
  fprintf(stderr, "done\n");

  plot_flow(universal_time, poltor_p);

  plot_flow_shell(poltor_p);

  /* write lat/lon maps of B field model */
  {
    const size_t nalt = 4;
    const double alt[] = { 250.0, 350.0, 450.0, 500.0 };
    char buf[2048];
    size_t i;

    for (i = 0; i < nalt; ++i)
      {
        sprintf(buf, "B_file_%.0f.dat", alt[i]);
        fprintf(stderr, "main: plotting fixed %.1f [km] altitude map of model B to %s...",
                alt[i], buf);
        plot_B(buf, poltor_p->R + alt[i], poltor_p);
        fprintf(stderr, "done\n");
      }
  }

  {
    const double alt = 350.0;
    const double lon = 130.0;
    const char *slice_file = "slice.dat";

    fprintf(stderr, "main: plotting latitude profile through %.1f [km] altitude, %.1f [deg] longitude to %s...",
            alt, lon, slice_file);
    print_slice(slice_file, 6371.2 + alt, lon * M_PI / 180.0, poltor_p);
    fprintf(stderr, "done\n");
  }

  poltor_free(poltor_p);

  return 0;
} /* main() */
