/*
 * stage1.c
 *
 * 1. Read tiegcm data file(s)
 * 2. For each time step t_k, invert B(t_k) grid for SH coefficients q_{nm}(t_k)
 * 3. Store SH coefficients in a matrix:
 *
 *      X_{ij} = q_i(t_j) where i = shidx(n,m)
 * 4. X matrix is output to a binary file
 *
 * ./stage1 <-i tiegcm_nc_file> [-o output_matrix_file]
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

#include <lapacke/lapacke.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include "apex.h"
#include "common.h"
#include "geo.h"
#include "green.h"
#include "magdata.h"
#include "oct.h"
#include "poltor.h"

#include "io.h"
#include "tiegcm.h"

#define USE_SYNTH_DATA              0

int
print_spectrum(const char *filename, const gsl_vector *c,
               const green_workspace *w)
{
  size_t nmax = w->nmax;
  size_t n;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_spectrum: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
  fprintf(fp, "# Field %zu: power (nT^2)\n", n++);

  for (n = 1; n <= nmax; ++n)
    {
      double sum = 0.0;
      int M = (int) n;
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m);
          double knm = gsl_vector_get(c, cidx);

          sum += knm * knm;
        }

      sum *= (n + 1.0);

      fprintf(fp, "%zu %.12e\n", n, sum);
    }

  fclose(fp);

  return 0;
}

/*
print_residuals()
  Print model residuals for a given timestamp

Inputs: filename - where to store residuals
        tidx     - time index
        A        - model matrix
        c        - model coefficients
        data     - TIEGCM data
*/

double
print_residuals(const char *filename, const size_t tidx,
                const gsl_matrix *A, const gsl_vector *c,
                const tiegcm_data *data)
{
  FILE *fp;
  gsl_vector *b = gsl_vector_alloc(A->size1); /* model prediction */
  size_t bidx = 0;
  size_t ilon, ilat, j;
  double rnorm = 0.0;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_residuals: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0.0;
    }

  /* compute b = A c */
  gsl_blas_dgemv(CblasNoTrans, 1.0, A, c, 0.0, b);

  j = 1;
  fprintf(fp, "# Field %zu: geographic longitude (degrees)\n", j++);
  fprintf(fp, "# Field %zu: geodetic latitude (degrees)\n", j++);
  fprintf(fp, "# Field %zu: TIEGCM B_x (nT)\n", j++);
  fprintf(fp, "# Field %zu: TIEGCM B_y (nT)\n", j++);
  fprintf(fp, "# Field %zu: TIEGCM B_z (nT)\n", j++);
  fprintf(fp, "# Field %zu: Modeled B_x (nT)\n", j++);
  fprintf(fp, "# Field %zu: Modeled B_y (nT)\n", j++);
  fprintf(fp, "# Field %zu: Modeled B_z (nT)\n", j++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM_BIDX(tidx, ilat, ilon, data);
          double B_model[3], B_data[3];

          B_data[0] = data->Bx[idx] * 1.0e9;
          B_data[1] = data->By[idx] * 1.0e9;
          B_data[2] = data->Bz[idx] * 1.0e9;

          for (j = 0; j < 3; ++j)
            {
              B_model[j] = gsl_vector_get(b, bidx++);

              /* update residual norm */
              rnorm = gsl_hypot(rnorm, B_data[j] - B_model[j]);
            }

          fprintf(fp, "%f %f %f %f %f %f %f %f\n",
                  data->glon[ilon],
                  data->glat[ilat],
                  data->Bx[idx] * 1.0e9,
                  data->By[idx] * 1.0e9,
                  data->Bz[idx] * 1.0e9,
                  B_model[0],
                  B_model[1],
                  B_model[2]);
        }
    }

  gsl_vector_free(b);
  fclose(fp);

  return rnorm;
}

/* plot current stream function grid */
int
print_chi(const char *filename, poltor_workspace *w)
{
  const double b = w->R + 110.0;
  const double lat_min = -80.0;
  const double lat_max = 80.0;
  const double lon_min = -180.0;
  const double lon_max = 180.0;
  const double phi_min = lon_min * M_PI / 180.0;
  const double phi_max = lon_max * M_PI / 180.0;
  const double theta_min = M_PI / 2.0 - lat_max * M_PI / 180.0;
  const double theta_max = M_PI / 2.0 - lat_min * M_PI / 180.0;
  double theta, phi;
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_chi: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: chi (external) (kA)\n", i++);

  fprintf(stderr, "print_chi: writing current stream function to %s...", filename);

  for (phi = phi_min; phi < phi_max; phi += 2.0 * M_PI / 180.0)
    {
      for (theta = theta_min; theta < theta_max; theta += 2.0 * M_PI / 180.0)
        {
          double chi; /* current stream function */

          /*poltor_eval_chi_ext(b, theta, phi, &chi, w);*/

          fprintf(fp, "%f %f %f\n",
                  phi * 180.0 / M_PI,
                  90.0 - theta * 180.0 / M_PI,
                  chi);
        }
      fprintf(fp, "\n");
    }

  fprintf(stderr, "done\n");

  fclose(fp);

  return 0;
}

magdata *
tiegcm_magdata(const size_t tidx, tiegcm_data *data)
{
  const size_t grid_size = data->nlon * data->nlat;
  magdata *mdata;
  magdata_datum datum;
  size_t ilat, ilon;
  apex_workspace *apex_p;

  mdata = magdata_alloc(grid_size, R_EARTH_KM);
  if (!mdata)
    return 0;

  apex_p = apex_alloc(2016);

  magdata_datum_init(&datum);

  datum.t = satdata_timet2epoch(data->t[tidx]);
  datum.flags = MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z;

  fprintf(stderr, "tiegcm_magdata: building magdata structure for time index %zu...", tidx);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = data->glon[ilon] * M_PI / 180.0;

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          double latd = data->glat[ilat] * M_PI / 180.0; /* geodetic latitude */
          double thetad = M_PI / 2.0 - latd;             /* geodetic colatitude */
          double r, latc, theta;                         /* geocentric radius, latitude and colatitude */
          double qdlat, alon, alat;
          size_t idx = TIEGCM_BIDX(tidx, ilat, ilon, data);

#if 0
          geodetic2geo(latd, 0.0, &latc, &r);
          theta = M_PI / 2.0 - latc;
#else
          theta = thetad;
          r = R_EARTH_KM;
#endif

          apex_transform_geodetic(thetad, phi, 0.0, &alon, &alat, &qdlat,
                                  NULL, NULL, NULL, apex_p);

          datum.r = r;
          datum.theta = theta;
          datum.phi = phi;
          datum.qdlat = qdlat;
          datum.B_nec[0] = data->Bx[idx] * 1.0e9;
          datum.B_nec[1] = data->By[idx] * 1.0e9;
          datum.B_nec[2] = data->Bz[idx] * 1.0e9;

          magdata_add(&datum, mdata);
        }
    }

  fprintf(stderr, "done\n");

  apex_free(apex_p);

  return mdata;
}

int
main_build_matrix(const magdata *mdata, green_workspace *green_p,
                  gsl_matrix *A)
{
  const size_t n = A->size1;
  const double eps = 1.0e-6;
  size_t rowidx = 0;
  size_t ilon, ilat;
  size_t i;

  for (i = 0; i < mdata->n; ++i)
    {
      double r = mdata->r[i];
      double theta = mdata->theta[i];
      double phi = mdata->phi[i];
      gsl_vector_view vx, vy, vz;

      if (theta < eps)
        theta = eps;
      if (theta > M_PI - eps)
        theta = M_PI - eps;

      vx = gsl_matrix_row(A, rowidx++);
      vy = gsl_matrix_row(A, rowidx++);
      vz = gsl_matrix_row(A, rowidx++);

      /* compute external Green's functions */
      green_calc_ext(r, theta, phi,
                     vx.vector.data,
                     vy.vector.data,
                     vz.vector.data,
                     green_p);
    }

  assert(rowidx == n);

  return 0;
}

/*
main_build_rhs()
  Construct RHS vector for a given time index
*/

int
main_build_rhs(const size_t tidx, const tiegcm_data *data,
               gsl_vector *b)
{
  const size_t n = b->size;
  size_t ilon, ilat;
  size_t rowidx = 0;

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM_BIDX(tidx, ilat, ilon, data);

          gsl_vector_set(b, rowidx++, data->Bx[idx] * 1.0e9);
          gsl_vector_set(b, rowidx++, data->By[idx] * 1.0e9);
          gsl_vector_set(b, rowidx++, data->Bz[idx] * 1.0e9);
        }
    }

  assert(rowidx == n);

  return 0;
}

int
lapack_lls(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * X)
{
  int s;
  lapack_int m = A->size1;
  lapack_int n = A->size2;
  lapack_int nrhs = B->size2;
  lapack_int lda = A->size1;
  lapack_int ldb = B->size1;
  lapack_int rank;
  lapack_int *jpvt = malloc(n * sizeof(lapack_int));
  gsl_matrix *work_A = gsl_matrix_alloc(A->size2, A->size1);
  gsl_matrix *work_B = gsl_matrix_alloc(B->size2, B->size1);
  double rcond = 1.0e-6;

  gsl_matrix_transpose_memcpy(work_A, A);
  gsl_matrix_transpose_memcpy(work_B, B);

  s = LAPACKE_dgelsy(LAPACK_COL_MAJOR,
                     m,
                     n,
                     nrhs,
                     work_A->data,
                     lda,
                     work_B->data,
                     ldb,
                     jpvt,
                     rcond,
                     &rank);

  /* store solution in X */
  {
    gsl_matrix_view m = gsl_matrix_submatrix(work_B, 0, 0, X->size2, X->size1);
    gsl_matrix_transpose_memcpy(X, &m.matrix);
  }

  gsl_matrix_free(work_A);
  gsl_matrix_free(work_B);
  free(jpvt);

  return s;
}

int
main_proc(const char *filename, const char *outfile_mat, tiegcm_data *data)
{
  int status;
  const char *res_file = "res.dat";
  const char *spectrum_file = "spectrum.s";
  const char *datamap_file = "datamap.dat";
  const size_t nmax = 60;
  green_workspace *green_p = green_alloc(nmax);
  const size_t n = 3 * data->nlon * data->nlat; /* number of residuals */
  const size_t p = green_p->nnm;                /* number of external coefficients */
  const size_t nrhs = data->nt;                 /* number of right hand sides */
  gsl_matrix *A = gsl_matrix_alloc(n, p);       /* least squares matrix */
  gsl_matrix *B = gsl_matrix_alloc(n, nrhs);    /* right hand sides */
  gsl_matrix *X = gsl_matrix_alloc(p, nrhs);    /* solution vectors */
  gsl_vector *r = gsl_vector_alloc(n);          /* residual vector */
  magdata *mdata;
  size_t k;
  FILE *fp;
  struct timeval tv0, tv1;

  fprintf(stderr, "main_proc: %zu observations per grid\n", n);
  fprintf(stderr, "main_proc: %zu SH model coefficients\n", p);
  fprintf(stderr, "main_proc: %zu timestamps\n", nrhs);

  /* store spatial locations in magdata structure - grid points are
   * the same for all timestamps t_k */
  mdata = tiegcm_magdata(0, data);

  /* print data map */
  fprintf(stderr, "main_proc: writing data map to %s...", datamap_file);
  magdata_map(datamap_file, mdata);
  fprintf(stderr, "done\n");

  /* construct least squares matrix (common for all timestamps) */
  fprintf(stderr, "main_proc: building least squares matrix A...");
  gettimeofday(&tv0, NULL);
  status = main_build_matrix(mdata, green_p, A);
  if (status)
    return status;
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fp = fopen(filename, "w");

  k = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", k++);
  fprintf(fp, "# Field %zu: q(1,0) (nT)\n", k++);
  fprintf(fp, "# Field %zu: q(1,1) (nT)\n", k++);
  fprintf(fp, "# Field %zu: q(2,0) (nT)\n", k++);
  fprintf(fp, "# Field %zu: q(2,1) (nT)\n", k++);
  fprintf(fp, "# Field %zu: q(2,2) (nT)\n", k++);

  fprintf(stderr, "main_proc: building rhs vectors...");
  gettimeofday(&tv0, NULL);

  /* construct right hand side vectors */
  for (k = 0; k < data->nt; ++k)
    {
      gsl_vector_view b = gsl_matrix_column(B, k);

      /* construct rhs vector for time t_k */
      main_build_rhs(k, data, &b.vector);

      if (k == 0)
        printv_octave(&b.vector, "b");
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* solve least squares system for all rhs vectors */
  fprintf(stderr, "main_proc: solving LS system with QR decomposition of A...");
  gettimeofday(&tv0, NULL);
  status = lapack_lls(A, B, X);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, s = %d)\n", time_diff(tv0, tv1), status);

  /* print spectrum of coefficients at time t_0 */
  {
    const size_t k = 0;
    gsl_vector_view x = gsl_matrix_column(X, k);
    fprintf(stderr, "main_proc: writing spectrum at t0 to %s...", spectrum_file);
    print_spectrum(spectrum_file, &x.vector, green_p);
    fprintf(stderr, "done\n");
  }

  /* print residuals at time t_0 */
  {
    const size_t k = 0;
    gsl_vector_view x = gsl_matrix_column(X, k);
    double rnorm;

    fprintf(stderr, "main_proc: writing residuals to %s...", res_file);
    rnorm = print_residuals(res_file, k, A, &x.vector, data);
    fprintf(stderr, "done (|| b - A x || = %.12e)\n", rnorm);
  }

  for (k = 0; k < data->nt; ++k)
    {
      size_t N;
#if 0
      gsl_vector_view b = gsl_matrix_column(B, k);
      gsl_vector_view x = gsl_matrix_column(X, k);

      /* compute r = b - A x */
      gsl_vector_memcpy(r, &b.vector);
      gsl_blas_dgemv(CblasNoTrans, -1.0, A, &x.vector, 1.0, r);

      fprintf(stderr, "main_proc: residual for timestamp (%zu/%zu): %.12e [nT]\n",
              k + 1, data->nt, gsl_blas_dnrm2(r));
#endif

      fprintf(fp, "%ld ", data->t[k]);

      for (N = 1; N <= 2; ++N)
        {
          int M = (int) N;
          int m;

          for (m = 0; m <= M; ++m)
            {
              size_t cidx = green_nmidx(N, m);
              double knm = gsl_matrix_get(X, cidx, k);

              fprintf(fp, "%f ", knm);
            }
        }

      putc('\n', fp);
      fflush(fp);
    }

  /* write matrix of solution vectors to output file */
  fprintf(stderr, "main_proc: writing solution matrix to %s...", outfile_mat);
  write_matrix(outfile_mat, X);
  fprintf(stderr, "done\n");

  green_free(green_p);
  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gsl_matrix_free(X);

  fclose(fp);

  fprintf(stderr, "main_proc: wrote qnm coefficients to %s\n", filename);

  return 0;
}

int
main(int argc, char *argv[])
{
  tiegcm_data *data = NULL;
  struct timeval tv0, tv1;
  char *infile = NULL;
  char *outfile = "qnm.txt";
  char *outfile_mat = "data/stage1_qnm.dat";

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:o:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'o':
            outfile_mat = optarg;
            break;

          default:
            break;
        }
    }

  while (optind < argc)
    {
      fprintf(stderr, "main: reading %s...", argv[optind]);
      gettimeofday(&tv0, NULL);

      data = tiegcm_read(argv[optind], data);
      if (!data)
        {
          fprintf(stderr, "main: error reading %s\n", argv[optind]);
          exit(1);
        }

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
              time_diff(tv0, tv1));

      ++optind;
    }

  main_proc(outfile, outfile_mat, data);

  tiegcm_free(data);

  return 0;
}
