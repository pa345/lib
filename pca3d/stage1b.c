/*
 * stage1b.c
 *
 * 1. Read TIEGCM 3D current grids for some time interval
 * 2. Perform spherical harmonic (Mie) decomposition on grids
 * 3. Store SH coefficients in output file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_blas.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_bsearch.h>

#include <magfield/magfield.h>
#include <magfield/magfield_eval.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

#define MAX_T        100000

#define TAPER_COEFFS 1

/*
taper_magfield()
  Taper high degree SH coefficients to try to reduce ringing.
We use a cosine taper set to 1 at lmin, and going to 0 at
lmax + 1. Don't make it go to 0 at lmax or the spectrum will
suddenly drop many orders of magnitude between lmax-1 and lmax
*/
int
taper_magfield(const size_t lmin, magfield_workspace * w)
{
  const size_t lmax = w->lmax;
  const double fac = M_PI / (2.0 * (lmax + 1 - lmin));
  size_t l;

  for (l = lmin + 1; l <= lmax; ++l)
    {
      double val = cos(fac * (l - lmin));
      double wl = val * val;
      int M = (int) GSL_MIN(l, w->mmax);
      int m;

      for (m = 0; m <= M; ++m)
        {
          size_t lmidx = magfield_lmidx(l, m, w->mmax);
          gsl_complex * pbelowlm = gsl_vector_complex_ptr(w->pbelow, lmidx);
          gsl_complex * pabovelm = gsl_vector_complex_ptr(w->pabove, lmidx);
          size_t i;

          GSL_REAL(*pbelowlm) *= wl; GSL_IMAG(*pbelowlm) *= wl;
          GSL_REAL(*pabovelm) *= wl; GSL_IMAG(*pabovelm) *= wl;

          for (i = 0; i < w->nr; ++i)
            {
              size_t idx = MAG_COEFIDX(i, lmidx, w);
              gsl_complex * qtlm = gsl_vector_complex_ptr(w->qtcoeff, idx);
              gsl_complex * qlm = gsl_vector_complex_ptr(w->qcoeff, idx);
              gsl_complex * plm = gsl_vector_complex_ptr(w->pcoeff, idx);

              GSL_REAL(*qtlm) *= wl; GSL_IMAG(*qtlm) *= wl;
              GSL_REAL(*qlm) *= wl; GSL_IMAG(*qlm) *= wl;
              GSL_REAL(*plm) *= wl; GSL_IMAG(*plm) *= wl;
            }
        }
    }

  return 0;
}

int
main_fill_grid(const size_t it, const tiegcm3d_data * data, magfield_workspace * w)
{
  int s = 0;
  size_t i, j, k;

  for (i = 0; i < w->nr; ++i)
    {
      for (j = 0; j < w->ntheta; ++j)
        {
          for (k = 0; k < w->nphi; ++k)
            {
              size_t idx = TIEGCM3D_IDX(it, i, j, k, data);

              magfield_grid_set(MAG_IDX_R, i, j, k, data->Jr[idx], w);
              magfield_grid_set(MAG_IDX_THETA, i, j, k, data->Jt[idx], w);
              magfield_grid_set(MAG_IDX_PHI, i, j, k, data->Jp[idx], w);
            }
        }
    }

  return s;
}

/*
print_residual()
  Print J grid for a fixed time and altitude from TIEGCM and magfield

Inputs: data - tiegcm data
*/

int
print_residual(const char *filename, const tiegcm3d_data *data, const int time_idx, const int ir,
               magfield_workspace * w, magfield_eval_workspace *eval_p, magfield_eval_workspace * eval_trunc_p)
{
  int s = 0;
  size_t i;
  size_t ilat, ilon;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_residual: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[time_idx], data->doy[time_idx] + data->ut[time_idx] / 24.0);
  fprintf(fp, "# Radius: %.2f (km) [%.2f km altitude]\n", data->r[ir] * 1.0e-3, data->r[ir] * 1.0e-3 - R_EARTH_KM);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: TIEGCM J_r (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: TIEGCM J_t (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: TIEGCM J_p (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: full magfield J_r (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: full magfield J_t (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: full magfield J_p (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: truncated magfield J_r (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: truncated magfield J_t (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: truncated magfield J_p (nA/m^2)\n", i++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = w->phi[ilon];

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM3D_IDX(time_idx, ir, ilat, ilon, data);
          double theta = w->theta[ilat];
          double J[4], J_trunc[4];

          magfield_eval_J(data->r[ir], theta, phi, J, eval_p);
          magfield_eval_J(data->r[ir], theta, phi, J_trunc, eval_trunc_p);

          fprintf(fp, "%8.4f %8.4f %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e\n",
                  data->glon[ilon],
                  data->glat[ilat],
                  data->Jr[idx] * 1.0e9,
                  data->Jt[idx] * 1.0e9,
                  data->Jp[idx] * 1.0e9,
                  J[0] * 1.0e9,
                  J[1] * 1.0e9,
                  J[2] * 1.0e9,
                  J_trunc[0] * 1.0e9,
                  J_trunc[1] * 1.0e9,
                  J_trunc[2] * 1.0e9);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/* copy poloidal/toroidal coefficients from src to dest, truncating if needed */
int
copy_coeffs(magfield_workspace * dest, const magfield_workspace * src)
{
  const size_t dest_lmax = dest->lmax;
  const size_t dest_mmax = dest->mmax;
  size_t i, l;

  assert(dest_lmax <= src->lmax);
  assert(dest_mmax <= src->mmax);
  assert(dest->nr == src->nr);

  for (l = dest->lmin; l <= dest_lmax; ++l)
    {
      int M = (int) GSL_MIN(l, dest_mmax);
      int m;

      for (m = 0; m <= M; ++m)
        {
          size_t dest_lmidx = magfield_lmidx(l, m, dest_mmax);
          size_t src_lmidx = magfield_lmidx(l, m, src->mmax);
          gsl_complex * dest_pabove = gsl_vector_complex_ptr(dest->pabove, dest_lmidx);
          gsl_complex * dest_pbelow = gsl_vector_complex_ptr(dest->pbelow, dest_lmidx);

          *dest_pabove = gsl_vector_complex_get(src->pabove, src_lmidx);
          *dest_pbelow = gsl_vector_complex_get(src->pbelow, src_lmidx);

          for (i = 0; i < src->nr; ++i)
            {
              size_t dest_idx = MAG_COEFIDX(i, dest_lmidx, dest);
              size_t src_idx = MAG_COEFIDX(i, src_lmidx, src);
              gsl_complex * dest_plmr = gsl_vector_complex_ptr(dest->pcoeff, dest_idx);
              gsl_complex * dest_qlmr = gsl_vector_complex_ptr(dest->qcoeff, dest_idx);
              gsl_complex * dest_qtlmr = gsl_vector_complex_ptr(dest->qtcoeff, dest_idx);

              *dest_plmr = gsl_vector_complex_get(src->pcoeff, src_idx);
              *dest_qlmr = gsl_vector_complex_get(src->qcoeff, src_idx);
              *dest_qtlmr = gsl_vector_complex_get(src->qtcoeff, src_idx);
            }
        }
    }

  return 0;
}

/*
main_proc()
  Compute SH decomposition of each radial shell of TIEGCM 3D current
values; store SH coefficients to disk

Inputs: prefix        - output file prefix
        tidx          - on input, starting time index
                        on output, ending time index
        tsteps        - (output) array of timestamps
        data          - TIEGCM data
        residual_file - (optional) write file of residuals between TIEGCM/magfield for first timestep
        w             - magfield workspace
        eval_p        - magfield eval workspace
        w_trunc       - magfield workspace with truncated SH degree/order
        eval_trunc_p  - magfield eval workspace for truncated SH degree/order

Return: success/error
*/

int
main_proc(const char * prefix, size_t * tidx, time_t * tsteps, const tiegcm3d_data * data, const char * residual_file,
          magfield_workspace * w, magfield_eval_workspace * eval_p, magfield_workspace * w_trunc,
          magfield_eval_workspace * eval_trunc_p)
{
  int s = 0;
  const size_t ttot = *tidx + data->nt;
  const double residual_alt = 110.0;
  int residual_ir = (int) bsearch_double(data->r, R_EARTH_M + residual_alt * 1.0e3, 0, data->nr - 1);
  char buf[2048];
  size_t i;

  /* ensure theta/phi grids from magfield match those from TIEGCM */

  for (i = 0; i < w->ntheta; ++i)
    gsl_test_rel((M_PI_2 - w->theta[i]) * 180.0 / M_PI, data->glat[i], 1.0e-10, "theta grid i=%zu", i);

  for (i = 0; i < w->nphi; ++i)
    gsl_test_rel(w->phi[i] * 180.0 / M_PI, data->glon[i], 1.0e-10, "phi grid i=%zu", i);

  for (i = 0; i < data->nt; ++i)
    {
      fprintf(stderr, "main_proc: filling current grid for timestep %zu/%zu\n", *tidx + 1, ttot);

      /* fill current grid for this timestep */
      main_fill_grid(i, data, w);

      /* compute SH decomposition */
      magfield_decomp(w);

      /* copy coefficients into w_trunc, truncating if needed */
      copy_coeffs(w_trunc, w);

#if TAPER_COEFFS
      taper_magfield(w_trunc->lmax - 30, w_trunc);
#endif

      /* write output file */
      sprintf(buf, "%s_%04zu.dat", prefix, *tidx + 1);
      fprintf(stderr, "main_proc: writing %s...", buf);
      magfield_write_SH(buf, w_trunc);
      fprintf(stderr, "done\n");

      tsteps[*tidx] = data->t[i];
      ++(*tidx);

      if (i == 0 && residual_file != NULL)
        {
          magfield_eval_init(w->qtcoeff, w->qcoeff, w->pcoeff, eval_p);
          magfield_eval_init(w_trunc->qtcoeff, w_trunc->qcoeff, w_trunc->pcoeff, eval_trunc_p);

          fprintf(stderr, "main_proc: writing residual for time index %zu to %s...", i, residual_file);
          print_residual(residual_file, data, i, residual_ir, w, eval_p, eval_trunc_p);
          fprintf(stderr, "done\n");
        }
    }

  return s;
}

int
main(int argc, char *argv[])
{
  const size_t lmin = 1;
  const double R = R_EARTH_M;
  size_t eval_lmax = 90;
  size_t eval_mmax = 30;
  char *file_prefix = PCA3D_STAGE1B_SH_PREFIX;
  char *residual_file = NULL;
  struct timeval tv0, tv1;
  magfield_workspace *w = NULL;
  magfield_workspace *w_trunc = NULL;
  magfield_eval_workspace * magfield_eval_p = NULL;
  magfield_eval_workspace * magfield_eval_trunc_p = NULL;
  size_t tidx = 0;
  time_t * tsteps = malloc(MAX_T * sizeof(time_t));

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "residual_file", required_argument, NULL, 'r' },
          { "lmax", required_argument, NULL, 'l' },
          { "mmax", required_argument, NULL, 'm' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "l:m:r:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'l':
            eval_lmax = (size_t) atoi(optarg);
            break;

          case 'm':
            eval_mmax = (size_t) atoi(optarg);
            break;

          case 'r':
            residual_file = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s [-l eval_lmax] [-m eval_mmax] [-r residual_file] file1.nc file2.nc ...\n", argv[0]);
            break;
        }
    }

  if (optind >= argc)
    {
      fprintf(stderr, "Usage: %s [-l eval_lmax] [-m eval_mmax] [-r residual_file] file1.nc file2.nc ...\n", argv[0]);
      exit(1);
    }

  while (optind < argc)
    {
      tiegcm3d_data *data;

      fprintf(stderr, "main: reading %s...", argv[optind]);
      gettimeofday(&tv0, NULL);

      data = tiegcm3d_read(argv[optind], NULL);
      if (!data)
        {
          fprintf(stderr, "main: error reading %s\n", argv[optind]);
          exit(1);
        }

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
              time_diff(tv0, tv1));

      if (w == NULL)
        {
          const size_t full_lmax = data->nlat - 1;
          const size_t full_mmax = data->nlon / 2 - 1;
          const size_t nr = data->nr;
          const size_t ntheta = data->nlat;
          const size_t nphi = data->nlon;
          magfield_params params;
          size_t i;

          if (eval_lmax == 0)
            eval_lmax = full_lmax;
          if (eval_mmax == 0)
            eval_mmax = full_mmax;

          params.lmin = lmin;
          params.lmax = full_lmax;
          params.mmax = full_mmax;
          params.nr = nr;
          params.ntheta = ntheta;
          params.nphi = nphi;
          params.R = R;
          params.grid_type = MAGFIELD_GAUSS;

          for (i = 0; i < nr; ++i)
            params.r[i] = data->r[i];

          fprintf(stderr, "main: allocating magfield workspace...");
          w = magfield_alloc(&params);
          magfield_eval_p = magfield_eval_alloc(&params);
          fprintf(stderr, "done\n");

          params.lmax = eval_lmax;
          params.mmax = eval_mmax;
          fprintf(stderr, "main: allocating truncated magfield workspace...");
          w_trunc = magfield_alloc(&params);
          magfield_eval_trunc_p = magfield_eval_alloc(&params);
          fprintf(stderr, "done\n");

          fprintf(stderr, "main: full lmax   = %zu\n", full_lmax);
          fprintf(stderr, "main: full mmax   = %zu\n", full_mmax);
          fprintf(stderr, "main: full nlm    = %zu\n", w->nlm);
          fprintf(stderr, "main: eval lmax   = %zu\n", eval_lmax);
          fprintf(stderr, "main: eval mmax   = %zu\n", eval_mmax);
          fprintf(stderr, "main: eval nlm    = %zu\n", w_trunc->nlm);
        }

      fprintf(stderr, "main: nr         = %zu\n", data->nr);
      fprintf(stderr, "main: nlat       = %zu\n", data->nlat);
      fprintf(stderr, "main: nlon       = %zu\n", data->nlon);

      main_proc(file_prefix, &tidx, tsteps, data, residual_file, w, magfield_eval_p, w_trunc, magfield_eval_trunc_p);

      tiegcm3d_free(data);

      ++optind;
    }

  /* write magfield parameter data and TIEGCM data to disk */
  {
    pca3d_data d;

    d.nt = tidx;
    d.t = tsteps;
    d.w = w_trunc;

    fprintf(stderr, "main: writing meta-data file %s...", PCA3D_STAGE1B_DATA);
    pca3d_write_data(PCA3D_STAGE1B_DATA, &d);
    fprintf(stderr, "done\n");
  }

  fprintf(stderr, "main: SH coefficients stored in %s\n", file_prefix);
  
  if (w)
    magfield_free(w);

  if (w_trunc)
    magfield_free(w_trunc);

  free(tsteps);

  return 0;
}
