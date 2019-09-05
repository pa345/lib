/*
 * stage1.c
 *
 * 1. Read tiegcm data file(s)
 * 2. For each time step t_k, invert B(t_k) grid for SH coefficients k_{nm}(t_k)
 * 3. Store k_{nm}(t) SH coefficients in a nnm-by-nt matrix:
 *
 *      X_{ij} = k_i(t_j) where i = shidx(n,m)
 * 5. X matrix is output to a binary file
 *
 * ./stage1 [-o binary_output_matrix_file] tiegcm1.nc tiegcm2.nc ...
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
#include <omp.h>

#include <lapacke/lapacke.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <mainlib/ml_bsearch.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_geo.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_apex.h>
#include <mainlib/ml_magdata.h>
#include <mainlib/ml_green.h>
#include <mainlib/ml_lapack_wrapper.h>

#include "poltor.h"

#include "io.h"
#include "pca.h"
#include "tiegcm.h"

#include "pca_common.c"

/* maximum external spherical harmonic degree and order for TIEGCM map inversion */
#define NMAX_EXT                    60
#define MMAX_EXT                    30

/* maximum internal spherical harmonic degree and order for TIEGCM map inversion */
#define NMAX_INT                    0
#define MMAX_INT                    0

#define TAPER_COEFFS                0

/* fine-tune which components of the surface B maps to fit */
#define FIT_X                       1
#define FIT_Y                       1
#define FIT_Z                       0

/*
print_residuals()
  Print model residuals for a given timestamp

Inputs: filename - where to store residuals
        tidx     - time index
        A        - model matrix
        c        - model coefficients
        data     - TIEGCM data
        green_ext - green workspace for external source
*/

double
print_residuals(const char *filename, const size_t tidx,
                const gsl_matrix *A, const gsl_vector *c,
                const tiegcm_data *data, green_workspace *green_ext)
{
  FILE *fp;
  gsl_vector *b = gsl_vector_alloc(A->size1); /* model prediction */
  size_t bidx = 0;
  size_t ilon, ilat, j;
  double rnorm = 0.0;
  const double eps = 1.0e-4;

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
  fprintf(fp, "# Grid for timestamp %ld\n", data->t[tidx]);
  fprintf(fp, "# Field %zu: geographic longitude (degrees)\n", j++);
  fprintf(fp, "# Field %zu: geodetic latitude (degrees)\n", j++);
  fprintf(fp, "# Field %zu: TIEGCM B_x (nT)\n", j++);
  fprintf(fp, "# Field %zu: TIEGCM B_y (nT)\n", j++);
  fprintf(fp, "# Field %zu: TIEGCM B_z (nT)\n", j++);
  fprintf(fp, "# Field %zu: Modeled B_x (nT)\n", j++);
  fprintf(fp, "# Field %zu: Modeled B_y (nT)\n", j++);
  fprintf(fp, "# Field %zu: Modeled B_z (nT)\n", j++);
  fprintf(fp, "# Field %zu: Modeled chi (kA)\n", j++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = data->glon[ilon] * M_PI / 180.0;

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM_BIDX(tidx, ilat, ilon, data);
          double B_model[3] = { 0.0, 0.0, 0.0 };
          double B_data[3], chi;
          double theta = M_PI / 2.0 - data->glat[ilat] * M_PI / 180.0;

          /* exclude poles */
          if (theta < eps)
            continue;
          if (theta > M_PI - eps)
            continue;

          B_data[0] = data->Bx[idx] * 1.0e9;
          B_data[1] = data->By[idx] * 1.0e9;
          B_data[2] = data->Bz[idx] * 1.0e9;

#if FIT_X
          B_model[0] = gsl_vector_get(b, bidx++);
          rnorm = gsl_hypot(rnorm, B_data[0] - B_model[0]);
#endif

#if FIT_Y
          B_model[1] = gsl_vector_get(b, bidx++);
          rnorm = gsl_hypot(rnorm, B_data[1] - B_model[1]);
#endif

#if FIT_Z
          B_model[2] = gsl_vector_get(b, bidx++);
          rnorm = gsl_hypot(rnorm, B_data[2] - B_model[2]);
#endif

          chi = green_eval_chi_ext(green_ext->R + 110.0, theta, phi, c, green_ext);

          fprintf(fp, "%8.4f %8.4f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
#if 0
                  wrap180(data->glon[ilon]),
#else
                  data->glon[ilon],
#endif
                  data->glat[ilat],
                  data->Bx[idx] * 1.0e9,
                  data->By[idx] * 1.0e9,
                  data->Bz[idx] * 1.0e9,
                  B_model[0],
                  B_model[1],
                  B_model[2],
                  chi);
        }

      fprintf(fp, "\n");
    }

  gsl_vector_free(b);
  fclose(fp);

  return rnorm;
}

/*
print_matlab()
  Print J_t, J_p components at given timestamp

Inputs: outdir   - where to store output
        tidx     - time index
        c        - model coefficients (knm)
        data     - TIEGCM data
        green_ext - green workspace for external source
*/

int
print_matlab(const char *outdir, const size_t tidx,
             const gsl_vector *c, green_workspace *green_ext)
{
  char buf[2048];
  const int max_threads = omp_get_max_threads();
  const size_t nphi = 360;
  const size_t ntheta = 180;
  const double b = green_ext->R + 110.0;
  const double nnm = green_ext->nnm;
  const double nmax = green_ext->nmax;
  const double mmax = green_ext->mmax;
  gsl_matrix *Jt = gsl_matrix_alloc(nphi, ntheta);
  gsl_matrix *Jp = gsl_matrix_alloc(nphi, ntheta);
  gsl_matrix *chi = gsl_matrix_alloc(nphi, ntheta);
  gsl_matrix *Bx = gsl_matrix_alloc(nphi, ntheta);
  gsl_matrix *By = gsl_matrix_alloc(nphi, ntheta);
  gsl_matrix *Bz = gsl_matrix_alloc(nphi, ntheta);
  gsl_vector *g = gsl_vector_alloc(c->size);
  gsl_vector *X = gsl_vector_alloc(nnm);
  gsl_vector *Y = gsl_vector_alloc(nnm);
  gsl_vector *Z = gsl_vector_alloc(nnm);
  green_workspace **green_p = malloc(max_threads * sizeof(green_workspace *));
  size_t i;

  for (i = 0; i < max_threads; ++i)
    green_p[i] = green_alloc(nmax, mmax, green_ext->R);

  /* convert knm to gnm */
  green_k2g(green_ext->R + 110.0, c, g, green_ext);

#pragma omp parallel for private(i)
  for (i = 0; i < nphi; ++i)
    {
      int thread_id = omp_get_thread_num();
      double phi = (i + 0.5) * M_PI / 180.0;
      size_t j;

      for (j = 0; j < ntheta; ++j)
        {
          double theta = (j + 0.5) * M_PI / 180.0;
          double K[3], B[3], chi1;

#if 1
          green_eval_sheet_int(b, theta, phi, g, K, green_p[thread_id]);
          gsl_matrix_set(Jt, i, j, -K[0]);
          gsl_matrix_set(Jp, i, j, K[1]);
#else
          chi1 = green_eval_chi_ext(b, theta, phi, c, green_p[thread_id]);
          gsl_matrix_set(chi, i, j, chi1);

          green_calc_ext(R_EARTH_KM, theta, phi, X->data, Y->data, Z->data, green_p[thread_id]);
          gsl_blas_ddot(X, c, &B[0]);
          gsl_blas_ddot(Y, c, &B[1]);
          gsl_blas_ddot(Z, c, &B[2]);

          gsl_matrix_set(Bx, i, j, B[0]);
          gsl_matrix_set(By, i, j, B[1]);
          gsl_matrix_set(Bz, i, j, B[2]);
#endif
        }
    }

#if 1
  sprintf(buf, "%s/Jt_%04zu.txt", outdir, tidx);
  print_octave(Jt, buf);

  sprintf(buf, "%s/Jp_%04zu.txt", outdir, tidx);
  print_octave(Jp, buf);

  sprintf(buf, "%s/chi_%04zu.txt", outdir, tidx);
  print_octave(chi, buf);
#else

  sprintf(buf, "%s/Bx_%04zu.txt", outdir, tidx);
  print_octave(Bx, buf);

  sprintf(buf, "%s/By_%04zu.txt", outdir, tidx);
  print_octave(By, buf);

  sprintf(buf, "%s/Bz_%04zu.txt", outdir, tidx);
  print_octave(Bz, buf);
#endif

  for (i = 0; i < max_threads; ++i)
    green_free(green_p[i]);

  free(green_p);
  gsl_matrix_free(Jt);
  gsl_matrix_free(Jp);
  gsl_matrix_free(chi);
  gsl_matrix_free(Bx);
  gsl_matrix_free(By);
  gsl_matrix_free(Bz);
  gsl_vector_free(g);
  gsl_vector_free(X);
  gsl_vector_free(Y);
  gsl_vector_free(Z);

  return 0;
}

magdata *
tiegcm_magdata(const size_t tidx, tiegcm_data *data)
{
  const double epoch = 2009.0;
  const size_t grid_size = data->nlon * data->nlat;
  magdata *mdata;
  magdata_datum datum;
  size_t ilat, ilon;
  apex_workspace *apex_p;

  mdata = magdata_alloc(grid_size, R_EARTH_KM);
  if (!mdata)
    return 0;

  apex_p = apex_alloc();

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

          apex_transform_geodetic(epoch, thetad, phi, 0.0, &alon, &alat, &qdlat,
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

size_t
main_build_matrix(const magdata *mdata, green_workspace *green_ext,
                  green_workspace *green_int, gsl_matrix *A)
{
  const size_t p_ext = green_ext->nnm;
  const size_t p_int = green_int->nnm;
  const double eps = 1.0e-4;
  gsl_vector *X_ext = gsl_vector_alloc(p_ext);
  gsl_vector *Y_ext = gsl_vector_alloc(p_ext);
  gsl_vector *Z_ext = gsl_vector_alloc(p_ext);
  gsl_vector *X_int, *Y_int, *Z_int;
  size_t rowidx = 0;
  size_t i;

  if (p_int > 0)
    {
      X_int = gsl_vector_alloc(p_int);
      Y_int = gsl_vector_alloc(p_int);
      Z_int = gsl_vector_alloc(p_int);
    }

  for (i = 0; i < mdata->n; ++i)
    {
      double r = mdata->r[i];
      double theta = mdata->theta[i];
      double phi = mdata->phi[i];
      gsl_vector_view v;

      /* exclude poles from fit */
      if (theta < eps)
        continue;
      if (theta > M_PI - eps)
        continue;

      /* compute external Green's functions */
      green_calc_ext(r, theta, phi, X_ext->data, Y_ext->data, Z_ext->data, green_ext);

      /* compute internal Green's functions */
      if (p_int > 0)
        green_calc_int(r, theta, phi, X_int->data, Y_int->data, Z_int->data, green_int);

#if FIT_X
      v = gsl_matrix_subrow(A, rowidx, 0, p_ext);
      gsl_vector_memcpy(&v.vector, X_ext);

      if (p_int > 0)
        {
          v = gsl_matrix_subrow(A, rowidx, p_ext, p_int);
          gsl_vector_memcpy(&v.vector, X_int);
        }

      ++rowidx;
#endif

#if FIT_Y
      v = gsl_matrix_subrow(A, rowidx, 0, p_ext);
      gsl_vector_memcpy(&v.vector, Y_ext);

      if (p_int > 0)
        {
          v = gsl_matrix_subrow(A, rowidx, p_ext, p_int);
          gsl_vector_memcpy(&v.vector, Y_int);
        }

      ++rowidx;
#endif

#if FIT_Z
      v = gsl_matrix_subrow(A, rowidx, 0, p_ext);
      gsl_vector_memcpy(&v.vector, Z_ext);

      if (p_int > 0)
        {
          v = gsl_matrix_subrow(A, rowidx, p_ext, p_int);
          gsl_vector_memcpy(&v.vector, Z_int);
        }

      ++rowidx;
#endif
    }

  gsl_vector_free(X_ext);
  gsl_vector_free(Y_ext);
  gsl_vector_free(Z_ext);

  if (p_int > 0)
    {
      gsl_vector_free(X_int);
      gsl_vector_free(Y_int);
      gsl_vector_free(Z_int);
    }

  return rowidx;
}

/*
main_build_rhs()
  Construct RHS vector for a given time index
*/

size_t
main_build_rhs(const size_t tidx, const tiegcm_data *data, gsl_vector *b)
{
  const double eps = 1.0e-4;
  size_t ilon, ilat;
  size_t rowidx = 0;
  green_workspace *green_p = green_alloc(1, 1, R_EARTH_KM);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM_BIDX(tidx, ilat, ilon, data);
          double theta = M_PI / 2.0 - data->glat[ilat] * M_PI / 180.0;

          /* exclude poles from fit */
          if (theta < eps)
            continue;
          if (theta > M_PI - eps)
            continue;

#if FIT_X
          gsl_vector_set(b, rowidx++, data->Bx[idx] * 1.0e9);
#endif

#if FIT_Y
          gsl_vector_set(b, rowidx++, data->By[idx] * 1.0e9);
#endif

#if FIT_Z
          gsl_vector_set(b, rowidx++, data->Bz[idx] * 1.0e9);
#endif
        }
    }

  green_free(green_p);

  return rowidx;
}

/*
convert_qnm()
  Convert k_{nm}(t) time series to q_{nm}(t)

Inputs: b       - radius of current shell (km)
        v       - vector of k_{nm}(t) for some time t, size nnm
        green_p - green workspace
*/

int
convert_qnm(const double b, gsl_vector * v, const green_workspace * green_p)
{
  const size_t nmax = green_p->nmax;
  const size_t mmax = green_p->mmax;
  const double ratio = b / green_p->R;
  double rterm = pow(ratio, -1.0); /* (b/R)^{n-2} */
  size_t n;

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;
      double nfac = (2.0 * n + 1.0) / (n + 1.0);

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m, green_p);
          double knm = gsl_vector_get(v, cidx);

          gsl_vector_set(v, cidx, nfac * rterm * knm);
        }

      /* (b/R)^{n-2} */
      rterm *= ratio;
    }

  return 0;
}

int
main_proc(const char *filename, const char *outfile_mat, const char *outdir_mat, tiegcm_data *data)
{
  int status;
  const char *res_file = "res.dat";
  const char *spectrum_file = "spectrum.s";
  const char *spectrum_taper_file = "spectrum_taper.s";
  const char *spectrum_azim_file = "spectrum_azim.s";
  const char *datamap_file = "datamap.dat";
  const size_t nmax_ext = NMAX_EXT;
  const size_t mmax_ext = GSL_MIN(nmax_ext, MMAX_EXT);
  const size_t nmax_int = NMAX_INT;
  const size_t mmax_int = GSL_MIN(nmax_int, MMAX_INT);
  const double R = R_EARTH_KM;
  green_workspace *green_ext = green_alloc(nmax_ext, mmax_ext, R);
  green_workspace *green_int = green_alloc(nmax_int, mmax_int, R);
  const size_t n = 3 * data->nlon * data->nlat; /* number of total residuals */
  const size_t p_ext = green_ext->nnm;          /* number of external coefficients */
  const size_t p_int = green_int->nnm;          /* number of internal coefficients */
  const size_t p = p_ext + p_int;               /* number of total coefficients */
  const size_t nrhs = data->nt;                 /* number of right hand sides */
  gsl_matrix *A = gsl_matrix_alloc(n, p);       /* least squares matrix */
  gsl_matrix *B = gsl_matrix_alloc(n, nrhs);    /* right hand sides */
  gsl_matrix *X = gsl_matrix_alloc(p, nrhs);    /* solution vectors */
  gsl_matrix *X_taper = gsl_matrix_alloc(p, nrhs);  /* tapered solution vectors to reduce ringing */
  gsl_vector *r = gsl_vector_alloc(n);          /* residual vector */
  magdata *mdata;
  size_t ndata;                                 /* number of residuals after excluding poles */
  gsl_matrix_view AA, BB;
  size_t k;
  FILE *fp;
  struct timeval tv0, tv1;
  int rank;

  fprintf(stderr, "main_proc: %zu longitude grid points\n", data->nlon);
  fprintf(stderr, "main_proc: %zu latitude grid points\n", data->nlat);
  fprintf(stderr, "main_proc: %zu observations per grid\n", n);
  fprintf(stderr, "main_proc: %zu external SH model coefficients\n", p_ext);
  fprintf(stderr, "main_proc: %zu internal SH model coefficients\n", p_int);
  fprintf(stderr, "main_proc: %zu total SH model coefficients\n", p);
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
  ndata = main_build_matrix(mdata, green_ext, green_int, A);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, %zu total data)\n", time_diff(tv0, tv1), ndata);

  fp = fopen(filename, "w");

  fprintf(stderr, "main_proc: building rhs vectors...");
  gettimeofday(&tv0, NULL);

  /* construct right hand side vectors */
  for (k = 0; k < data->nt; ++k)
    {
      gsl_vector_view b = gsl_matrix_column(B, k);
      size_t nrows;

      /* construct rhs vector for time t_k */
      nrows = main_build_rhs(k, data, &b.vector);

      assert(nrows == ndata);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* solve least squares system for all rhs vectors */
  AA = gsl_matrix_submatrix(A, 0, 0, ndata, p);
  BB = gsl_matrix_submatrix(B, 0, 0, ndata, nrhs);

  fprintf(stderr, "main_proc: solving LS systems with QR decomposition of A (%zu-by-%zu)...", ndata, p);
  gettimeofday(&tv0, NULL);
  status = lapack_lls(&AA.matrix, &BB.matrix, X, &rank);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, s = %d, rank = %d)\n",
          time_diff(tv0, tv1), status, rank);

  /* form residual matrix R = B - A*X */
  {
    gsl_matrix *R = gsl_matrix_alloc(ndata, nrhs);
    double rnorm_max = 0.0;

    gsl_matrix_memcpy(R, &BB.matrix);

    fprintf(stderr, "main_proc: computing residual matrix R = B - AX...");
    gettimeofday(&tv0, NULL);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &AA.matrix, X, -1.0, R);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

    for (k = 0; k < nrhs; ++k)
      {
        gsl_vector_view v = gsl_matrix_column(R, k);
        double norm = gsl_blas_dnrm2(&v.vector);
        if (norm > rnorm_max)
          rnorm_max = norm;
      }

    fprintf(stderr, "main_proc: maximum residual norm = %.12e\n", rnorm_max);

    gsl_matrix_free(R);
  }

  /* taper high degree coefficients to correct TIEGCM ringing */
  {
    gsl_matrix_memcpy(X_taper, X);

#if TAPER_COEFFS
    taper_knm(30, X_taper, green_ext);
#endif
  }

  /* print spectrum of external coefficients at time t_0 */
  {
    const size_t k = 0;
    gsl_vector_view x = gsl_matrix_column(X, k);

    printv_octave(&x.vector, "c");

    fprintf(stderr, "main_proc: writing spectrum at t0 to %s...", spectrum_file);
    green_print_spectrum(spectrum_file, &x.vector, green_ext);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main_proc: writing azimuth spectrum at t0 to %s...", spectrum_azim_file);
    green_print_spectrum_azim(spectrum_azim_file, &x.vector, green_ext);
    fprintf(stderr, "done\n");

    x = gsl_matrix_column(X_taper, k);
    fprintf(stderr, "main_proc: writing tapered spectrum at t0 to %s...", spectrum_taper_file);
    green_print_spectrum(spectrum_taper_file, &x.vector, green_ext);
    fprintf(stderr, "done\n");
  }

  /* print residuals for a given timestamp */
  {
#if 0
    const time_t unix_time = 1240660800; /* Apr 25 2009 12:00:00 UTC */
    const size_t k = bsearch_timet(data->t, unix_time, 0, data->nt - 1);
#elif 0
    const time_t unix_time = 1240617600; /* Apr 25 2009 12:00:00 UTC */
    const size_t k = bsearch_timet(data->t, unix_time, 0, data->nt - 1);
#elif 1
    const time_t unix_time = 1231808400; /* Jan 13 2009 01:00:00 UTC */
    const size_t k = bsearch_timet(data->t, unix_time, 0, data->nt - 1);
#elif 1
    const size_t k = 0;
#endif
    gsl_vector_view x = gsl_matrix_column(X_taper, k);
    double rnorm;

    fprintf(stderr, "main_proc: writing residuals for timestamp %ld to %s...", data->t[k], res_file);
    rnorm = print_residuals(res_file, k, &AA.matrix, &x.vector, data, green_ext);
    fprintf(stderr, "done (|| b - A x || = %.12e)\n", rnorm);
  }

#if 0
  /* convert k_{nm}(t) to q_{nm}(t) */
  fprintf(stderr, "main_proc: converting knm to qnm...");
  gettimeofday(&tv0, NULL);
  {
    for (k = 0; k < data->nt; ++k)
      {
        gsl_vector_view Xk = gsl_matrix_column(X, k);
        convert_qnm(b, &Xk.vector, green_p);
      }
  }
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
#endif

  k = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", k++);
  fprintf(fp, "# Field %zu: k(1,0) (nT)\n", k++);
  fprintf(fp, "# Field %zu: k(1,1) (nT)\n", k++);
  fprintf(fp, "# Field %zu: k(2,0) (nT)\n", k++);
  fprintf(fp, "# Field %zu: k(2,1) (nT)\n", k++);
  fprintf(fp, "# Field %zu: k(2,2) (nT)\n", k++);

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
              size_t cidx = green_nmidx(N, m, green_ext);
              double knm = gsl_matrix_get(X, cidx, k);

              fprintf(fp, "%f ", knm);
            }
        }

      putc('\n', fp);
      fflush(fp);
    }

  /* write matrix of solution vectors to output file */
  fprintf(stderr, "main_proc: writing solution matrix to %s...", outfile_mat);
  pca_write_matrix(outfile_mat, X);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main_proc: writing misc data to %s...", PCA_STAGE1_DATA);
  pca_write_data(PCA_STAGE1_DATA, nmax_ext, mmax_ext, data);
  fprintf(stderr, "done\n");

  if (outdir_mat)
    {
      char buf[2048];
      FILE *fp;

      sprintf(buf, "%s/ts", outdir_mat);
      fp = fopen(buf, "w");

      for (k = 0; k < data->nt; ++k)
        {
          gsl_vector_view x = gsl_matrix_column(X_taper, k);

          fprintf(stderr, "main_proc: writing matlab output %zu/%zu...", k, data->nt);
          print_matlab(outdir_mat, k, &x.vector, green_ext);
          fprintf(stderr, "done\n");

          fprintf(fp, "%ld\n", data->t[k]);
        }

      fclose(fp);
    }

  green_free(green_ext);
  gsl_matrix_free(A);
  gsl_matrix_free(B);
  gsl_matrix_free(X);
  gsl_matrix_free(X_taper);

  fclose(fp);

  fprintf(stderr, "main_proc: wrote knm coefficients to %s\n", filename);

  return 0;
}

int
main(int argc, char *argv[])
{
  tiegcm_data *data = NULL;
  struct timeval tv0, tv1;
  char *outfile = "knm.txt";
  char *outdir_mat = NULL;
  char *outfile_mat = PCA_STAGE1_KNM;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "o:m:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'o':
            outfile_mat = optarg;
            break;

          case 'm':
            outdir_mat = optarg;
            break;

          default:
            break;
        }
    }

  if (optind >= argc)
    {
      fprintf(stderr, "Usage: %s [-o binary_matrix_output_file] file1.nc file2.nc ...\n",
              argv[0]);
      exit(1);
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

  main_proc(outfile, outfile_mat, outdir_mat, data);

  tiegcm_free(data);

  return 0;
}
