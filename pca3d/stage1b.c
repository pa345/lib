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

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_test.h>

#include <common/common.h>

#include <magfield/magfield.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

#define MAX_T        100000

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

              magfield_current_set(MAG_IDX_R, i, j, k, data->Jr[idx], w);
              magfield_current_set(MAG_IDX_THETA, i, j, k, data->Jt[idx], w);
              magfield_current_set(MAG_IDX_PHI, i, j, k, data->Jp[idx], w);
            }
        }
    }

  return s;
}

/*
main_proc()
  Compute SH decomposition of each radial shell of TIEGCM 3D current
values; store SH coefficients to disk

Inputs: prefix - output file prefix
        tidx   - on input, starting time index
                 on output, ending time index
        tsteps - (output) array of timestamps
        data   - TIEGCM data
        w      - magfield workspace

Return: success/error
*/

int
main_proc(const char * prefix, size_t * tidx, time_t * tsteps, const tiegcm3d_data * data, magfield_workspace * w)
{
  int s = 0;
  const size_t ttot = *tidx + data->nt;
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

      /* write output file */
      sprintf(buf, "%s_%03zu.dat", prefix, *tidx + 1);
      magfield_write_SH(buf, w);

      tsteps[*tidx] = data->t[i];
      ++(*tidx);
    }

  return s;
}

int
main(int argc, char *argv[])
{
  const size_t lmin = 1;
  const size_t lmax = 60;
  const double R = R_EARTH_KM;
  char *file_prefix = PCA3D_STAGE1B_SH_PREFIX;
  struct timeval tv0, tv1;
  magfield_workspace *w = NULL;
  size_t tidx = 0;
  time_t * tsteps = malloc(MAX_T * sizeof(time_t));

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          default:
            fprintf(stderr, "Usage: %s file1.nc file2.nc ...\n", argv[0]);
            break;
        }
    }

  if (optind >= argc)
    {
      fprintf(stderr, "Usage: %s file1.nc file2.nc ...\n",
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
          const size_t mmax = GSL_MIN(30, data->nlon / 2 - 1);
          const size_t nr = data->nr;
          const size_t ntheta = data->nlat;
          const size_t nphi = data->nlon;
          const double rmin = data->r[0];
          const double rmax = data->r[nr - 1];
          magfield_params params;

          params.lmin = lmin;
          params.lmax = lmax;
          params.mmax = mmax;
          params.nr = nr;
          params.ntheta = ntheta;
          params.nphi = nphi;
          params.rmin = rmin;
          params.rmax = rmax;
          params.R = R;
          params.grid_type = MAGFIELD_GAUSS;

          fprintf(stderr, "main: allocating magfield workspace...");
          w = magfield_alloc(&params);
          fprintf(stderr, "done\n");

          /* set radial points which use non-fixed spacing */
          magfield_set_r(data->r, w);
        }

      fprintf(stderr, "main: nr         = %zu\n", data->nr);
      fprintf(stderr, "main: nlat       = %zu\n", data->nlat);
      fprintf(stderr, "main: nlon       = %zu\n", data->nlon);

      main_proc(file_prefix, &tidx, tsteps, data, w);

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
