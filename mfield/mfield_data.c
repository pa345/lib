/*
 * mfield_data.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_att.h>
#include <mainlib/ml_common.h>

#include "mfield_data.h"

static int mfield_data_print_satellite(const char *dir_prefix, const gsl_vector *wts_spatial,
                                       const size_t nsource, const magdata * mptr, size_t * index);
static int mfield_data_print_observatory(const char *dir_prefix, const gsl_vector *wts_spatial,
                                         const magdata * mptr, size_t * index);
static int mfield_data_print_observatory_SV(const char *dir_prefix, const gsl_vector *wts_spatial,
                                            const magdata * mptr, size_t * index);

/*
mfield_data_alloc()
  Allocate mfield_data_workspace

Inputs: nsources - number of data sources (satellites,observatories)
        params   - data parameters

Return: pointer to workspace
*/

mfield_data_workspace *
mfield_data_alloc(const size_t nsources, const mfield_data_parameters *params)
{
  mfield_data_workspace *w;

  w = calloc(1, sizeof(mfield_data_workspace));
  if (!w)
    return 0;

  w->nsources = nsources;
  w->params = *params;

  w->t0 = malloc(nsources * sizeof(double));
  w->t1 = malloc(nsources * sizeof(double));
  if (!w->t0 || !w->t1)
    {
      mfield_data_free(w);
      return 0;
    }

  w->mdata = calloc(nsources, sizeof(magdata *));
  if (!w->mdata)
    {
      mfield_data_free(w);
      return 0;
    }

  w->t_scale = calloc(nsources, sizeof(double *));
  w->t_year = calloc(nsources, sizeof(double *));
  if (!w->t_scale || !w->t_year)
    {
      mfield_data_free(w);
      return 0;
    }

  w->rstat_workspace_p = gsl_rstat_alloc();
  if (!w->rstat_workspace_p)
    {
      mfield_data_free(w);
      return 0;
    }

  w->t_mu = -1.0;
  w->t_sigma = -1.0;
  w->t0_data = -1.0;
  w->t1_data = -1.0;

  return w;
}

void
mfield_data_free(mfield_data_workspace *w)
{
  if (w->t0)
    free(w->t0);

  if (w->t1)
    free(w->t1);

  if (w->rstat_workspace_p)
    gsl_rstat_free(w->rstat_workspace_p);

  if (w->mdata)
    {
      size_t i;

      for (i = 0; i < w->nsources; ++i)
        {
          magdata *mdata = w->mdata[i];

          if (mdata)
            magdata_free(mdata);
        }

      free(w->mdata);
    }

  if (w->t_scale)
    {
      size_t i;

      for (i = 0; i < w->nsources; ++i)
        {
          if (w->t_scale[i])
            free(w->t_scale[i]);
        }

      free(w->t_scale);
    }

  if (w->t_year)
    {
      size_t i;

      for (i = 0; i < w->nsources; ++i)
        {
          if (w->t_year[i])
            free(w->t_year[i]);
        }

      free(w->t_year);
    }

  free(w);
}

/*
mfield_data_filter_time()
  Flag any data points outside of [tmin,tmax] with
MAGDATA_FLG_DISCARD

Inputs: tmin - minimum time (decimal year)
        tmax - maximum time (decimal year)
        w    - workspace

Return: number of data flagged

Notes:
1) tmin/tmax can be set to -1 to exclude them from the test
*/

size_t
mfield_data_filter_time(const double tmin, const double tmax,
                        mfield_data_workspace *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          double t = satdata_epoch2year(mptr->t[j]);

#if 0 /*XXX*/
          if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY_SV && t > 2009.0)
#else
          if ((tmin > 0.0 && t < tmin) ||
              (tmax > 0.0 && t > tmax))
#endif
            {
              mptr->flags[j] |= MAGDATA_FLG_DISCARD;
              ++cnt;
            }
        }
    }

  return cnt;
}

/*
mfield_data_filter_align()
  We are not fitting alignment parameters (fit_align is 0),
so discard any data which is marked alignment only

Inputs: w    - workspace

Return: number of data flagged
*/

size_t
mfield_data_filter_align(mfield_data_workspace *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          if ((mptr->flags[j] & MAGDATA_FLG_FIT_ALIGN) &&
              !(mptr->flags[j] & MAGDATA_FLG_FIT_MF))
            {
              mptr->flags[j] |= MAGDATA_FLG_DISCARD;
              ++cnt;
            }
        }
    }

  return cnt;
}

/*
mfield_data_filter_fluxcal()
  When fitting fluxgate calibration parameters, scalar field data
requires also a VFM vector, so flag any data which has a scalar
measurement but not a vector

Inputs: w - workspace

Return: number of data flagged

Notes:
1) This function should only be called if fit_align and fit_fluxcal are true
*/

size_t
mfield_data_filter_fluxcal(mfield_data_workspace *w)
{
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_ExistScalar(mptr->flags[j]) && !MAGDATA_ExistVector(mptr->flags[j]))
            {
              mptr->flags[j] |= MAGDATA_FLG_DISCARD;
              ++cnt;
            }
        }
    }

  return cnt;
}

/*
mfield_data_filter_comp()
  Discard any data according to config file flags

Inputs: w    - workspace

Return: number of data flagged
*/

size_t
mfield_data_filter_comp(mfield_data_workspace *w)
{
  const mfield_data_parameters *params = &(w->params);
  const double qdlat_cutoff = params->fit_seplat ? params->qdlat_fit_cutoff : 100.0;
  size_t cnt = 0;
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      /* if this dataset is EEJ measurements of the magnetic equator, don't filter components */
      if (mptr->global_flags & MAGDATA_GLOBFLG_EEJ_MAGEQ)
        continue;

      for (j = 0; j < mptr->n; ++j)
        {
          double qdlat = mptr->qdlat[j];

#if 0
          /*XXX fit Swarm scalar data and observatory vector data */
          if (mptr->global_flags & MAGDATA_GLOBFLG_SATELLITE)
            mptr->flags[j] &= ~(MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z);
#endif

          if (fabs(qdlat) <= qdlat_cutoff)
            {
              /* select components for mid/low latitudes */

              if (!params->fit_X)
                mptr->flags[j] &= ~MAGDATA_FLG_X;

              if (!params->fit_Y)
                mptr->flags[j] &= ~MAGDATA_FLG_Y;

              if (!params->fit_Z)
                mptr->flags[j] &= ~MAGDATA_FLG_Z;

              if (!params->fit_F)
                mptr->flags[j] &= ~MAGDATA_FLG_F;

              if (!params->fit_DXDT)
                mptr->flags[j] &= ~MAGDATA_FLG_DXDT;

              if (!params->fit_DYDT)
                mptr->flags[j] &= ~MAGDATA_FLG_DYDT;

              if (!params->fit_DZDT)
                mptr->flags[j] &= ~MAGDATA_FLG_DZDT;

              if (!params->fit_DX_NS)
                mptr->flags[j] &= ~MAGDATA_FLG_DX_NS;

              if (!params->fit_DY_NS)
                mptr->flags[j] &= ~MAGDATA_FLG_DY_NS;

              if (!params->fit_DZ_NS)
                mptr->flags[j] &= ~MAGDATA_FLG_DZ_NS;

              if (!params->fit_DF_NS)
                mptr->flags[j] &= ~MAGDATA_FLG_DF_NS;

              if (!params->fit_DX_EW)
                mptr->flags[j] &= ~MAGDATA_FLG_DX_EW;

              if (!params->fit_DY_EW)
                mptr->flags[j] &= ~MAGDATA_FLG_DY_EW;

              if (!params->fit_DZ_EW)
                mptr->flags[j] &= ~MAGDATA_FLG_DZ_EW;

              if (!params->fit_DF_EW)
                mptr->flags[j] &= ~MAGDATA_FLG_DF_EW;
            }
          else
            {
              /* select components for high-latitudes */

              /* don't fit X/Y at high-latitudes, including gradients */
              mptr->flags[j] &= ~(MAGDATA_FLG_X | MAGDATA_FLG_Y);
              mptr->flags[j] &= ~(MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS);
              mptr->flags[j] &= ~(MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW);

              if (!params->fit_Z_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_Z;

              if (!params->fit_F_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_F;

              /*
               * dX/dt and dY/dt we still fit at high-latitudes (unless toggled off), since the observatory
               * time series data are still relatively clean
               */
              if (!params->fit_DXDT)
                mptr->flags[j] &= ~MAGDATA_FLG_DXDT;

              if (!params->fit_DYDT)
                mptr->flags[j] &= ~MAGDATA_FLG_DYDT;

              if (!params->fit_DZDT)
                mptr->flags[j] &= ~MAGDATA_FLG_DZDT;

              if (!params->fit_DZ_NS_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_DZ_NS;

              if (!params->fit_DF_NS_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_DF_NS;

              if (!params->fit_DZ_EW_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_DZ_EW;

              if (!params->fit_DF_EW_highlat)
                mptr->flags[j] &= ~MAGDATA_FLG_DF_EW;
            }
        }
    }

  return cnt;
}

/*
mfield_data_filter_observatory()
  When computing crustal biases, we need a certain minimum number of
measurements. Check each observatory to see that it has sufficient data,
otherwise flag the whole observatory dataset to remove it from the modeling.

Inputs: w - data workspace
*/

size_t
mfield_data_filter_observatory(mfield_data_workspace *w)
{
  const size_t min_data = 50;
  size_t i, j;
  size_t nflagged = 0;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);
      size_t ndata[3] = { 0, 0, 0 };

      if (!(mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY))
        continue;

      /* count number of data for each component */
      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_ExistX(mptr->flags[j]))
            ++ndata[0];

          if (MAGDATA_ExistY(mptr->flags[j]))
            ++ndata[1];

          if (MAGDATA_ExistZ(mptr->flags[j]))
            ++ndata[2];
        }

      /* if available data is below threshold, flag all data */
      if ((ndata[0] > 0 && ndata[0] < min_data) ||
          (ndata[1] > 0 && ndata[1] < min_data) ||
          (ndata[2] > 0 && ndata[2] < min_data))
        {
          ++nflagged;

          for (j = 0; j < mptr->n; ++j)
            mptr->flags[j] &= ~(MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z);
        }
    }

  return nflagged;
}

/*
mfield_data_compact()
  Remove magdata entries which have the discard flag set or will
not be used in the model
*/

int
mfield_data_compact(mfield_data_workspace * w)
{
  int s = 0;
  const size_t keep_flags = MAGDATA_FLG_FIT_MF | MAGDATA_FLG_FIT_ALIGN;
  magdata ** new_data = malloc(w->nsources * sizeof(magdata *));
  size_t i;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);
      new_data[i] = magdata_compact(keep_flags, mptr);
      magdata_free(mptr);
    }

  w->mdata = new_data;

  return s;
}

/*
mfield_data_init_align()
  Change alignment conventions at run-time for different satellites
*/

int
mfield_data_init_align(mfield_data_workspace *w)
{
  int s = 0;
  size_t i;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      if (mptr->global_flags & MAGDATA_GLOBFLG_SWARM)
        mptr->align_flags = ATT_TYPE_EULER_XYZ;
      else if (mptr->global_flags & MAGDATA_GLOBFLG_CHAMP)
        mptr->align_flags = ATT_TYPE_EULER_XYZ;
      else if (mptr->global_flags & MAGDATA_GLOBFLG_DMSP)
        mptr->align_flags = ATT_TYPE_MRP;
      else if (mptr->global_flags & MAGDATA_GLOBFLG_CRYOSAT)
        mptr->align_flags = ATT_TYPE_MRP;
    }

  return s;
}

/*
mfield_data_init()
  Compute mean and stddev of timestamps minus epoch for later time scaling

w_i = t_i - epoch
t_mu = mean(w_i)
t_sigma = stddev(w_i)

Inputs: w - workspace

Return: success/error

Notes:
1) w->t_mu and w->t_sigma are updated with timestamp mean/stddev in years

2) w->t0_data is initialized to the timestamp of the first data point (CDF_EPOCH)

3) w->t1_data is initialized to the timestamp of the last data point (CDF_EPOCH)

3) w->t0 and w->t1 are initialized to the first/last timestamps of each satellite

4) initialize w->t_scale and w->t_year arrays
*/

int
mfield_data_init(mfield_data_workspace *w)
{
  int s = 0;
  size_t i, j;

  gsl_rstat_reset(w->rstat_workspace_p);

  w->t0_data = 1.0e15;
  w->t1_data = -1.0e15;
  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      magdata_t(&(w->t0[i]), &(w->t1[i]), mptr);

      if (mptr->n != 0)
        {
          if (w->t0[i] > 0.0)
            w->t0_data = GSL_MIN(w->t0_data, w->t0[i]);

          if (w->t1[i] > 0.0)
            w->t1_data = GSL_MAX(w->t1_data, w->t1[i]);
        }

      for (j = 0; j < mptr->n; ++j)
        {
          double t;

          if (mptr->flags[j] & MAGDATA_FLG_DISCARD)
            continue;

          t = satdata_epoch2year(mptr->t[j]) - w->params.epoch;
          gsl_rstat_add(t, w->rstat_workspace_p);
        }
    }

  w->t_mu = gsl_rstat_mean(w->rstat_workspace_p);
  w->t_sigma = gsl_rstat_sd(w->rstat_workspace_p);

  if (w->t_sigma == 0.0)
    {
      /* this can happen for a fixed time grid like EMAG2 */
      w->t_mu = 0.0;
      w->t_sigma = 1.0;
    }

  /* initialize t_scale and t_year */
  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      w->t_scale[i] = malloc(GSL_MAX(mptr->n, 1) * sizeof(double));
      w->t_year[i] = malloc(GSL_MAX(mptr->n, 1) * sizeof(double));

      for (j = 0; j < mptr->n; ++j)
        {
          double * t_scale = w->t_scale[i];
          double * t_year = w->t_year[i];

          t_scale[j] = (mptr->t[j] - w->t_mu) / w->t_sigma;
          t_year[j] = epoch2year(mptr->t[j]);
        }
    }

  return s;
}

/*
mfield_data_epoch()
  Compute epoch of input data by averaging all timestamps
*/

double
mfield_data_epoch(mfield_data_workspace *w)
{
  /* initialize t_mu and t_sigma */
  mfield_data_init(w);

  return w->t_mu + w->params.epoch;
} /* mfield_data_epoch() */

int
mfield_data_map(const char *dir_prefix, const mfield_data_workspace *w)
{
  int s = 0;
  size_t i;
  char buf[2048];

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      sprintf(buf, "%s/datamap%zu", dir_prefix, i);
      fprintf(stderr, "mfield_data_map: printing spatial coverage of satellite %zu to %s (MF data)...", i, buf);
      magdata_map(buf, MAGDATA_FLG_FIT_MF, mptr);
      fprintf(stderr, "done\n");

      sprintf(buf, "%s/datamap%zu_align", dir_prefix, i);
      fprintf(stderr, "mfield_data_map: printing spatial coverage of satellite %zu to %s (alignment data)...", i, buf);
      magdata_map(buf, MAGDATA_FLG_FIT_ALIGN, mptr);
    }

  return s;
}

/*
mfield_data_print()
  Print out all data which will be used for main field modeling
*/

int
mfield_data_print(const char *dir_prefix, const gsl_vector *wts_spatial,
                  const mfield_data_workspace *w)
{
  int s = 0;
  size_t i;
  size_t idx = 0;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      if (mptr->global_flags & MAGDATA_GLOBFLG_SATELLITE)
        mfield_data_print_satellite(dir_prefix, wts_spatial, i, mptr, &idx);
      else if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
        mfield_data_print_observatory(dir_prefix, wts_spatial, mptr, &idx);
      else if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY_SV)
        mfield_data_print_observatory_SV(dir_prefix, wts_spatial, mptr, &idx);
    }

  assert(idx == wts_spatial->size);

  return s;
}

magdata *
mfield_data_ptr(const size_t idx, const mfield_data_workspace *w)
{
  if (idx >= w->nsources)
    {
      fprintf(stderr, "mfield_data_ptr: invalid index: %zu\n", idx);
      return 0;
    }

  return w->mdata[idx];
}

/*
mfield_data_print_satellite()
  Print satellite data used for main field modeling

Inputs: dir_prefix  - directory prefix
        wts_spatial - spatial weights
        nsource     - number of data source
        mptr        - magdata
        index       - (input/output) index into wts_spatial
*/

static int
mfield_data_print_satellite(const char *dir_prefix, const gsl_vector *wts_spatial,
                            const size_t nsource, const magdata * mptr, size_t * index)
{
  int s = 0;
  const char *fmtstr = "%ld %.8f %.4f %.4f %.4f %.4f %.3f %.4f %.4f\n";
  const char *fmtstr_grad = "%ld %.8f %.4f %.4f %.4f %.4f %.3f %.4f %.4f %.4f %.4f\n";
  const size_t n = 12; /* number of components to print */
  FILE *fp[12];
  char buf[2048];
  size_t j, k;
  size_t idx = *index;

  sprintf(buf, "%s/data%zu_X.dat", dir_prefix, nsource);
  fp[0] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_Y.dat", dir_prefix, nsource);
  fp[1] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_Z.dat", dir_prefix, nsource);
  fp[2] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_F.dat", dir_prefix, nsource);
  fp[3] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_DX_NS.dat", dir_prefix, nsource);
  fp[4] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_DY_NS.dat", dir_prefix, nsource);
  fp[5] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_DZ_NS.dat", dir_prefix, nsource);
  fp[6] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_DF_NS.dat", dir_prefix, nsource);
  fp[7] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_DX_EW.dat", dir_prefix, nsource);
  fp[8] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_DY_EW.dat", dir_prefix, nsource);
  fp[9] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_DZ_EW.dat", dir_prefix, nsource);
  fp[10] = fopen(buf, "w");

  sprintf(buf, "%s/data%zu_DF_EW.dat", dir_prefix, nsource);
  fp[11] = fopen(buf, "w");

  for (j = 0; j < n; ++j)
    {
      if (fp[j] == NULL)
        {
          fprintf(stderr, "mfield_data_print_satellite: fp[%zu] is NULL\n", j);
          return -1;
        }
    }

  /* header line */
  fprintf(fp[0], "# X vector data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[1], "# Y vector data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[2], "# Z vector data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[3], "# F scalar data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[4], "# DX gradient (N/S) vector data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[5], "# DY gradient (N/S) vector data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[6], "# DZ gradient (N/S) vector data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[7], "# DF gradient (N/S) scalar data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[8], "# DX gradient (E/W) vector data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[9], "# DY gradient (E/W) vector data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[10], "# DZ gradient (E/W) vector data for MF modeling (satellite %zu)\n", nsource);
  fprintf(fp[11], "# DF gradient (E/W) scalar data for MF modeling (satellite %zu)\n", nsource);

  for (j = 0; j < n; ++j)
    {
      k = 1;
      fprintf(fp[j], "# Field %zu: timestamp (UT seconds since 1970-01-01)\n", k++);
      fprintf(fp[j], "# Field %zu: time (decimal year)\n", k++);
      fprintf(fp[j], "# Field %zu: longitude (degrees)\n", k++);
      fprintf(fp[j], "# Field %zu: geocentric latitude (degrees)\n", k++);
      fprintf(fp[j], "# Field %zu: QD latitude (degrees)\n", k++);
      fprintf(fp[j], "# Field %zu: geocentric radius (km)\n", k++);
      fprintf(fp[j], "# Field %zu: spatial weight factor\n", k++);
    }

  fprintf(fp[0], "# Field %zu: X vector measurement (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y vector measurement (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z vector measurement (nT)\n", k);
  fprintf(fp[3], "# Field %zu: F scalar measurement (nT)\n", k);
  fprintf(fp[4], "# Field %zu: X vector measurement (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y vector measurement (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z vector measurement (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F scalar measurement (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X vector measurement (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y vector measurement (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z vector measurement (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F scalar measurement (nT)\n", k);
  ++k;

  fprintf(fp[0], "# Field %zu: X a priori model (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y a priori model (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z a priori model (nT)\n", k);
  fprintf(fp[3], "# Field %zu: F a priori model (nT)\n", k);
  fprintf(fp[4], "# Field %zu: X a priori model (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y a priori model (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z a priori model (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F a priori model (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X a priori model (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y a priori model (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z a priori model (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F a priori model (nT)\n", k);
  ++k;

  fprintf(fp[4], "# Field %zu: X vector measurement at N/S gradient point (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y vector measurement at N/S gradient point (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z vector measurement at N/S gradient point (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F scalar measurement at N/S gradient point (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X vector measurement at E/W gradient point (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y vector measurement at E/W gradient point (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z vector measurement at E/W gradient point (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F scalar measurement at E/W gradient point (nT)\n", k);
  ++k;

  fprintf(fp[4], "# Field %zu: X a priori model at N/S gradient point (nT)\n", k);
  fprintf(fp[5], "# Field %zu: Y a priori model at N/S gradient point (nT)\n", k);
  fprintf(fp[6], "# Field %zu: Z a priori model at N/S gradient point (nT)\n", k);
  fprintf(fp[7], "# Field %zu: F a priori model at N/S gradient point (nT)\n", k);
  fprintf(fp[8], "# Field %zu: X a priori model at E/W gradient point (nT)\n", k);
  fprintf(fp[9], "# Field %zu: Y a priori model at E/W gradient point (nT)\n", k);
  fprintf(fp[10], "# Field %zu: Z a priori model at E/W gradient point (nT)\n", k);
  fprintf(fp[11], "# Field %zu: F a priori model at E/W gradient point (nT)\n", k);
  ++k;

  for (j = 0; j < mptr->n; ++j)
    {
      double t = satdata_epoch2year(mptr->t[j]);
      time_t unix_time = satdata_epoch2timet(mptr->t[j]);
      double phi = wrap180(mptr->phi[j] * 180.0 / M_PI);
      double lat = 90.0 - mptr->theta[j] * 180.0 / M_PI;
      double qdlat = mptr->qdlat[j];
      double r = mptr->r[j];
      double B[4], B_grad[4];
      double B_model[4], B_grad_model[4];

      if (MAGDATA_Discarded(mptr->flags[j]))
        continue;

      B[0] = mptr->Bx_nec[j];
      B[1] = mptr->By_nec[j];
      B[2] = mptr->Bz_nec[j];
      B[3] = mptr->F[j];

      B_model[0] = mptr->Bx_model[j];
      B_model[1] = mptr->By_model[j];
      B_model[2] = mptr->Bz_model[j];
      B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

      B_grad[0] = mptr->Bx_nec_ns[j];
      B_grad[1] = mptr->By_nec_ns[j];
      B_grad[2] = mptr->Bz_nec_ns[j];
      B_grad[3] = mptr->F_ns[j];

      B_grad_model[0] = mptr->Bx_model_ns[j];
      B_grad_model[1] = mptr->By_model_ns[j];
      B_grad_model[2] = mptr->Bz_model_ns[j];
      B_grad_model[3] = gsl_hypot3(B_grad_model[0], B_grad_model[1], B_grad_model[2]);

      if ((j > 0) && (mptr->flags[j] & MAGDATA_FLG_TRACK_START))
        {
          size_t k;

          for (k = 0; k < n; ++k)
            fprintf(fp[k], "\n\n");
        }

      if (MAGDATA_ExistX(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);

          if (MAGDATA_FitMF(mptr->flags[j]))
            fprintf(fp[0], fmtstr, unix_time, t, phi, lat, qdlat, r, wj, B[0], B_model[0]);
        }

      if (MAGDATA_ExistY(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);

          if (MAGDATA_FitMF(mptr->flags[j]))
            fprintf(fp[1], fmtstr, unix_time, t, phi, lat, qdlat, r, wj, B[1], B_model[1]);
        }

      if (MAGDATA_ExistZ(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);

          if (MAGDATA_FitMF(mptr->flags[j]))
            fprintf(fp[2], fmtstr, unix_time, t, phi, lat, qdlat, r, wj, B[2], B_model[2]);
        }

      if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);
          fprintf(fp[3], fmtstr, unix_time, t, phi, lat, qdlat, r, wj, B[3], B_model[3]);
        }

      if (MAGDATA_ExistDX_NS(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);

          if (MAGDATA_FitMF(mptr->flags[j]))
            fprintf(fp[4], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[0], B_model[0], B_grad[0], B_grad_model[0]);
        }

      if (MAGDATA_ExistDY_NS(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);

          if (MAGDATA_FitMF(mptr->flags[j]))
            fprintf(fp[5], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[1], B_model[1], B_grad[1], B_grad_model[1]);
        }

      if (MAGDATA_ExistDZ_NS(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);

          if (MAGDATA_FitMF(mptr->flags[j]))
            fprintf(fp[6], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[2], B_model[2], B_grad[2], B_grad_model[2]);
        }

      if (MAGDATA_ExistDF_NS(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);
          fprintf(fp[7], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[3], B_model[3], B_grad[3], B_grad_model[3]);
        }

      if (MAGDATA_ExistDX_EW(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);

          if (MAGDATA_FitMF(mptr->flags[j]))
            fprintf(fp[8], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[0], B_model[0], B_grad[0], B_grad_model[0]);
        }

      if (MAGDATA_ExistDY_EW(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);

          if (MAGDATA_FitMF(mptr->flags[j]))
            fprintf(fp[9], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[1], B_model[1], B_grad[1], B_grad_model[1]);
        }

      if (MAGDATA_ExistDZ_EW(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);

          if (MAGDATA_FitMF(mptr->flags[j]))
            fprintf(fp[10], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[2], B_model[2], B_grad[2], B_grad_model[2]);
        }

      if (MAGDATA_ExistDF_EW(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);
          fprintf(fp[11], fmtstr_grad, unix_time, t, phi, lat, qdlat, r, wj, B[3], B_model[3], B_grad[3], B_grad_model[3]);
        }
    }

  for (j = 0; j < n; ++j)
    fclose(fp[j]);

  *index = idx;

  return s;
}

/*
mfield_data_print_observatory()
  Print observatory data used for main field modeling

Inputs: dir_prefix  - directory prefix
        wts_spatial - spatial weights
        mptr        - magdata
        index       - (input/output) index into wts_spatial
*/

static int
mfield_data_print_observatory(const char *dir_prefix, const gsl_vector *wts_spatial,
                              const magdata * mptr, size_t * index)
{
  int s = 0;
  const char *fmtstr = "%ld %8.4f %6.3f %10.4f\n";
  const size_t n = 3; /* number of components to print */
  FILE *fp[3];
  char buf[2048];
  size_t j, k;
  size_t idx = *index;

  sprintf(buf, "%s/obs/%s_X.dat", dir_prefix, mptr->name);
  fp[0] = fopen(buf, "w");

  sprintf(buf, "%s/obs/%s_Y.dat", dir_prefix, mptr->name);
  fp[1] = fopen(buf, "w");

  sprintf(buf, "%s/obs/%s_Z.dat", dir_prefix, mptr->name);
  fp[2] = fopen(buf, "w");

  for (j = 0; j < n; ++j)
    {
      if (fp[j] == NULL)
        {
          fprintf(stderr, "mfield_data_print_observatory: fp[%zu] is NULL\n", j);
          return -1;
        }
    }

  /* print header */
  fprintf(fp[0], "# %s observatory\n# X vector data for MF modeling\n", mptr->name);
  fprintf(fp[1], "# %s observatory\n# Y vector data for MF modeling\n", mptr->name);
  fprintf(fp[2], "# %s observatory\n# Z vector data for MF modeling\n", mptr->name);

  for (j = 0; j < n; ++j)
    {
      fprintf(fp[j], "# Radius:    %.4f [km]\n", mptr->r[0]);
      fprintf(fp[j], "# Longitude: %.4f [deg]\n", mptr->phi[0] * 180.0 / M_PI);
      fprintf(fp[j], "# Latitude:  %.4f [deg]\n", 90.0 - mptr->theta[0] * 180.0 / M_PI);

      k = 1;
      fprintf(fp[j], "# Field %zu: timestamp (UT seconds since 1970-01-01)\n", k++);
      fprintf(fp[j], "# Field %zu: QD latitude (degrees)\n", k++);
      fprintf(fp[j], "# Field %zu: spatial weight factor\n", k++);
    }

  fprintf(fp[0], "# Field %zu: X vector measurement (nT)\n", k);
  fprintf(fp[1], "# Field %zu: Y vector measurement (nT)\n", k);
  fprintf(fp[2], "# Field %zu: Z vector measurement (nT)\n", k);

  for (j = 0; j < mptr->n; ++j)
    {
      time_t unix_time = satdata_epoch2timet(mptr->t[j]);

      if (MAGDATA_Discarded(mptr->flags[j]))
        continue;

      if (!MAGDATA_FitMF(mptr->flags[j]))
        continue;

      if (MAGDATA_ExistX(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);
          fprintf(fp[0], fmtstr, unix_time, mptr->qdlat[j], wj, mptr->Bx_nec[j]);
        }

      if (MAGDATA_ExistY(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);
          fprintf(fp[1], fmtstr, unix_time, mptr->qdlat[j], wj, mptr->By_nec[j]);
        }

      if (MAGDATA_ExistZ(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);
          fprintf(fp[2], fmtstr, unix_time, mptr->qdlat[j], wj, mptr->Bz_nec[j]);
        }
    }

  for (j = 0; j < n; ++j)
    fclose(fp[j]);

  *index = idx;

  return s;
}

/*
mfield_data_print_observatory_SV()
  Print observatory data used for secular variation modeling

Inputs: dir_prefix  - directory prefix
        wts_spatial - spatial weights
        mptr        - magdata
        index       - (input/output) index into wts_spatial
*/

static int
mfield_data_print_observatory_SV(const char *dir_prefix, const gsl_vector *wts_spatial,
                                 const magdata * mptr, size_t * index)
{
  int s = 0;
  const char *fmtstr = "%ld %8.4f %6.3f %10.4f\n";
  const size_t n = 3; /* number of components to print */
  FILE *fp[3];
  char buf[2048];
  size_t j, k;
  size_t idx = *index;

  sprintf(buf, "%s/obs/%s_DXDT.dat", dir_prefix, mptr->name);
  fp[0] = fopen(buf, "w");

  sprintf(buf, "%s/obs/%s_DYDT.dat", dir_prefix, mptr->name);
  fp[1] = fopen(buf, "w");

  sprintf(buf, "%s/obs/%s_DZDT.dat", dir_prefix, mptr->name);
  fp[2] = fopen(buf, "w");

  for (j = 0; j < n; ++j)
    {
      if (fp[j] == NULL)
        {
          fprintf(stderr, "mfield_data_print_observatory: fp[%zu] is NULL\n", j);
          return -1;
        }
    }

  /* print header */
  fprintf(fp[0], "# %s observatory\n# dX/dt vector data for MF modeling\n", mptr->name);
  fprintf(fp[1], "# %s observatory\n# dY/dt vector data for MF modeling\n", mptr->name);
  fprintf(fp[2], "# %s observatory\n# dZ/dt vector data for MF modeling\n", mptr->name);

  for (j = 0; j < n; ++j)
    {
      fprintf(fp[j], "# Radius:    %.4f [km]\n", mptr->r[0]);
      fprintf(fp[j], "# Longitude: %.4f [deg]\n", mptr->phi[0] * 180.0 / M_PI);
      fprintf(fp[j], "# Latitude:  %.4f [deg]\n", 90.0 - mptr->theta[0] * 180.0 / M_PI);

      k = 1;
      fprintf(fp[j], "# Field %zu: timestamp (UT seconds since 1970-01-01)\n", k++);
      fprintf(fp[j], "# Field %zu: QD latitude (degrees)\n", k++);
      fprintf(fp[j], "# Field %zu: spatial weight factor\n", k++);
    }

  fprintf(fp[0], "# Field %zu: dX/dt vector measurement (nT/year)\n", k);
  fprintf(fp[1], "# Field %zu: dY/dt vector measurement (nT/year)\n", k);
  fprintf(fp[2], "# Field %zu: dZ/dt vector measurement (nT/year)\n", k);

  for (j = 0; j < mptr->n; ++j)
    {
      time_t unix_time = satdata_epoch2timet(mptr->t[j]);

      if (MAGDATA_Discarded(mptr->flags[j]))
        continue;

      if (!MAGDATA_FitMF(mptr->flags[j]))
        continue;

      if (MAGDATA_ExistDXDT(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);
          fprintf(fp[0], fmtstr, unix_time, mptr->qdlat[j], wj, mptr->dXdt_nec[j]);
        }

      if (MAGDATA_ExistDYDT(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);
          fprintf(fp[1], fmtstr, unix_time, mptr->qdlat[j], wj, mptr->dYdt_nec[j]);
        }

      if (MAGDATA_ExistDZDT(mptr->flags[j]))
        {
          double wj = gsl_vector_get(wts_spatial, idx++);
          fprintf(fp[2], fmtstr, unix_time, mptr->qdlat[j], wj, mptr->dZdt_nec[j]);
        }
    }

  for (j = 0; j < n; ++j)
    fclose(fp[j]);

  *index = idx;

  return s;
}
