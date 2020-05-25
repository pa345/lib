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
               magfield_workspace * w, magfield_eval_workspace *eval_p)
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
  fprintf(fp, "# Field %zu: magfield J_r (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: magfield J_t (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: magfield J_p (nA/m^2)\n", i++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = w->phi[ilon];

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM3D_IDX(time_idx, ir, ilat, ilon, data);
          double theta = w->theta[ilat];
          double J[4];

          magfield_eval_J(data->r[ir], theta, phi, J, eval_p);

          fprintf(fp, "%8.4f %8.4f %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e\n",
                  data->glon[ilon],
                  data->glat[ilat],
                  data->Jr[idx] * 1.0e9,
                  data->Jt[idx] * 1.0e9,
                  data->Jp[idx] * 1.0e9,
                  J[0] * 1.0e9,
                  J[1] * 1.0e9,
                  J[2] * 1.0e9);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
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

Return: success/error
*/

int
main_proc(const char * prefix, size_t * tidx, time_t * tsteps, const tiegcm3d_data * data, const char * residual_file,
          magfield_workspace * w, magfield_eval_workspace * eval_p)
{
  int s = 0;
  const size_t ttot = *tidx + data->nt;
  const double residual_alt = 310.0;
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

#if TAPER_COEFFS
      taper_magfield(30, w);
#endif

      /* write output file */
      sprintf(buf, "%s_%03zu.dat", prefix, *tidx + 1);
      fprintf(stderr, "main_proc: writing %s...", buf);
      magfield_write_SH(buf, w);
      fprintf(stderr, "done\n");

      tsteps[*tidx] = data->t[i];
      ++(*tidx);

      if (i == 0 && residual_file != NULL)
        {
          magfield_eval_init(w->qtcoeff, w->qcoeff, w->pcoeff, eval_p);

          fprintf(stderr, "main_proc: writing residual for time index %zu to %s...", i, residual_file);
          print_residual(residual_file, data, i, residual_ir, w, eval_p);
          fprintf(stderr, "done\n");
        }
    }

  return s;
}

int
main(int argc, char *argv[])
{
  const size_t lmin = 1;
  const size_t lmax = 60;
  const size_t mmax = 30;
  const double R = R_EARTH_M;
  char *file_prefix = PCA3D_STAGE1B_SH_PREFIX;
  char *residual_file = NULL;
  struct timeval tv0, tv1;
  magfield_workspace *w = NULL;
  magfield_eval_workspace * magfield_eval_p = NULL;
  size_t tidx = 0;
  time_t * tsteps = malloc(MAX_T * sizeof(time_t));

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "residual_file", required_argument, NULL, 'r' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "r:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'r':
            residual_file = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s [-r residual_file] file1.nc file2.nc ...\n", argv[0]);
            break;
        }
    }

  if (optind >= argc)
    {
      fprintf(stderr, "Usage: %s [-r residual_file] file1.nc file2.nc ...\n",
              argv[0]);
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
          const size_t nr = data->nr;
          const size_t ntheta = data->nlat;
          const size_t nphi = data->nlon;
          magfield_params params;
          size_t i;

          params.lmin = lmin;
          params.lmax = lmax;
          params.mmax = GSL_MIN(mmax, data->nlon / 2 - 1);
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

          fprintf(stderr, "main: lmax       = %zu\n", params.lmax);
          fprintf(stderr, "main: mmax       = %zu\n", params.mmax);
        }

      fprintf(stderr, "main: nr         = %zu\n", data->nr);
      fprintf(stderr, "main: nlat       = %zu\n", data->nlat);
      fprintf(stderr, "main: nlon       = %zu\n", data->nlon);

      main_proc(file_prefix, &tidx, tsteps, data, residual_file, w, magfield_eval_p);

      tiegcm3d_free(data);

      ++optind;
    }

  /* write magfield parameter data and TIEGCM data to disk */
  {
    pca3d_data d;

    d.nt = tidx;
    d.t = tsteps;
    d.w = w;

    fprintf(stderr, "main: writing meta-data file %s...", PCA3D_STAGE1B_DATA);
    pca3d_write_data(PCA3D_STAGE1B_DATA, &d);
    fprintf(stderr, "done\n");
  }

  fprintf(stderr, "main: SH coefficients stored in %s\n", file_prefix);
  
  if (w)
    magfield_free(w);

  free(tsteps);

  return 0;
}
