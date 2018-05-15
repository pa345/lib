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
        data   - TIEGCM data

Return: success/error
*/

int
main_proc(const char * prefix, const tiegcm3d_data * data)
{
  int s = 0;
  const size_t lmin = 1;
  const size_t lmax = 60;
  const size_t mmax = GSL_MIN(30, data->nlon / 2 - 1);
  const size_t nr = data->nr;
  const size_t ntheta = data->nlat;
  const size_t nphi = data->nlon;
  const double rmin = data->r[0];
  const double rmax = data->r[nr - 1];
  const double R = R_EARTH_KM;
  char buf[2048];
  magfield_workspace *w;
  magfield_params params;
  size_t i;

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

  fprintf(stderr, "main_proc: allocating magfield workspace...");
  w = magfield_alloc(&params);
  fprintf(stderr, "done\n");

  /* set radial points which use non-fixed spacing */
  magfield_set_r(data->r, w);

  /* ensure theta/phi grids from magfield match those from TIEGCM */

  for (i = 0; i < w->ntheta; ++i)
    gsl_test_rel((M_PI_2 - w->theta[i]) * 180.0 / M_PI, data->glat[i], 1.0e-10, "theta grid i=%zu", i);

  for (i = 0; i < w->nphi; ++i)
    gsl_test_rel(w->phi[i] * 180.0 / M_PI, data->glon[i], 1.0e-10, "phi grid i=%zu", i);

  /* write magfield parameter data and TIEGCM data to disk */
  {
    pca3d_data d;

    d.nt = data->nt;
    d.t = data->t;
    d.w = w;
    pca3d_write_data(PCA3D_STAGE1B_DATA, &d);
  }

  for (i = 0; i < data->nt; ++i)
    {
      fprintf(stderr, "main_proc: filling current grid for timestep %zu/%zu\n", i + 1, data->nt);

      /* fill current grid for this timestep */
      main_fill_grid(i, data, w);

      /* compute SH decomposition */
      magfield_decomp(w);

      /* write output file */
      sprintf(buf, "%s_%03zu.dat", prefix, i + 1);
      magfield_write_SH(buf, w);
    }

  magfield_free(w);

  return s;
}

int
main(int argc, char *argv[])
{
  char *infile = NULL;
  char *file_prefix = PCA3D_STAGE1B_SH_PREFIX;
  struct timeval tv0, tv1;
  tiegcm3d_data *data;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file>\n", argv[0]);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file>\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = tiegcm3d_read(infile, NULL);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
          time_diff(tv0, tv1));

  fprintf(stderr, "main: nr         = %zu\n", data->nr);
  fprintf(stderr, "main: nlat       = %zu\n", data->nlat);
  fprintf(stderr, "main: nlon       = %zu\n", data->nlon);

  main_proc(file_prefix, data);

  return 0;
}
