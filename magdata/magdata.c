/*
 * magdata.c
 *
 * Calling sequence:
 * 1. magdata_alloc  - allocate data structure
 * 2. magdata_add    - add points to data structure
 * 3. magdata_init   - construct spatial weight histogram
 * 4. magdata_calc   - calculate spatial weights for all points
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <obsdata/obsdata.h>

#include <gsl/gsl_math.h>

#include <common/bsearch.h>
#include <common/common.h>
#include <common/interp.h>

#include "magdata.h"
#include "track.h"
#include "track_weight.h"

static int magdata_luhr(const double ne, const double Te, const double Ti,
                        const double B_nec[4], double b_luhr[4]);

/*
magdata_alloc()
  Allocate magdata workspace

Inputs: n     - maximum number of data to be added later
        R     - reference radius (km)
*/

magdata *
magdata_alloc(const size_t n, const double R)
{
  magdata *data;
  const size_t ntheta = 100;
  const size_t nphi = 100;

  data = calloc(1, sizeof(magdata));
  if (!data)
    return 0;

  data->weight_workspace_p = track_weight_alloc(ntheta, nphi);
  if (!data->weight_workspace_p)
    {
      fprintf(stderr, "magdata_alloc: failed to allocate workspaces\n");
      magdata_free(data);
      return 0;
    }

  magdata_realloc(n, R, data);

  data->n = 0;
  data->R = R;
  data->nx = 0;
  data->ny = 0;
  data->nz = 0;
  data->nf = 0;
  data->ndx = 0;
  data->ndy = 0;
  data->ndz = 0;
  data->ndf = 0;
  data->nvec = 0;
  data->ngrad = 0;
  data->nres = 0;

  data->rmin = 1.0e9;
  data->rmax = -1.0e9;

  data->euler_flags = 0;
  data->global_flags = 0;

  return data;
}

magdata *
magdata_realloc(const size_t n, const double R, magdata *data)
{
  if (!data)
    return magdata_alloc(n, R);

  if (n <= data->ntot)
    return data; /* nothing to do */

  data->t = realloc(data->t, n * sizeof(double));
  data->ts = realloc(data->ts, n * sizeof(double));
  data->r = realloc(data->r, n * sizeof(double));
  data->theta = realloc(data->theta, n * sizeof(double));
  data->phi = realloc(data->phi, n * sizeof(double));
  data->qdlat = realloc(data->qdlat, n * sizeof(double));
  data->MLT = realloc(data->MLT, n * sizeof(double));

  data->Bx_nec = realloc(data->Bx_nec, n * sizeof(double));
  data->By_nec = realloc(data->By_nec, n * sizeof(double));
  data->Bz_nec = realloc(data->Bz_nec, n * sizeof(double));
  data->Bx_vfm = realloc(data->Bx_vfm, n * sizeof(double));
  data->By_vfm = realloc(data->By_vfm, n * sizeof(double));
  data->Bz_vfm = realloc(data->Bz_vfm, n * sizeof(double));
  data->Bx_model = realloc(data->Bx_model, n * sizeof(double));
  data->By_model = realloc(data->By_model, n * sizeof(double));
  data->Bz_model = realloc(data->Bz_model, n * sizeof(double));
  data->F = realloc(data->F, n * sizeof(double));
  data->dXdt_nec = realloc(data->dXdt_nec, n * sizeof(double));
  data->dYdt_nec = realloc(data->dYdt_nec, n * sizeof(double));
  data->dZdt_nec = realloc(data->dZdt_nec, n * sizeof(double));
  data->q = realloc(data->q, 4 * n * sizeof(double));
  data->lt = realloc(data->lt, n * sizeof(double));
  data->lt_eq = realloc(data->lt_eq, n * sizeof(double));

  data->satdir = realloc(data->satdir, n * sizeof(int));
  data->ne = realloc(data->ne, n * sizeof(double));

  data->t_ns = realloc(data->t_ns, n * sizeof(double));
  data->ts_ns = realloc(data->ts_ns, n * sizeof(double));
  data->r_ns = realloc(data->r_ns, n * sizeof(double));
  data->theta_ns = realloc(data->theta_ns, n * sizeof(double));
  data->phi_ns = realloc(data->phi_ns, n * sizeof(double));
  data->qdlat_ns = realloc(data->qdlat_ns, n * sizeof(double));
  data->Bx_nec_ns = realloc(data->Bx_nec_ns, n * sizeof(double));
  data->By_nec_ns = realloc(data->By_nec_ns, n * sizeof(double));
  data->Bz_nec_ns = realloc(data->Bz_nec_ns, n * sizeof(double));
  data->Bx_vfm_ns = realloc(data->Bx_vfm_ns, n * sizeof(double));
  data->By_vfm_ns = realloc(data->By_vfm_ns, n * sizeof(double));
  data->Bz_vfm_ns = realloc(data->Bz_vfm_ns, n * sizeof(double));
  data->Bx_model_ns = realloc(data->Bx_model_ns, n * sizeof(double));
  data->By_model_ns = realloc(data->By_model_ns, n * sizeof(double));
  data->Bz_model_ns = realloc(data->Bz_model_ns, n * sizeof(double));
  data->F_ns = realloc(data->F_ns, n * sizeof(double));
  data->q_ns = realloc(data->q_ns, 4 * n * sizeof(double));
  data->lt_ns = realloc(data->lt_ns, n * sizeof(double));
  data->lt_eq_ns = realloc(data->lt_eq_ns, n * sizeof(double));

  data->flags = realloc(data->flags, n * sizeof(size_t));
  data->index = realloc(data->index, n * sizeof(size_t));
  data->weights = realloc(data->weights, n * sizeof(double));

  if (!data->t || !data->ts || !data->r || !data->theta || !data->phi || !data->qdlat || !data->flags ||
      !data->Bx_nec || !data->By_nec || !data->Bz_nec || !data->Bx_vfm || !data->By_vfm ||
      !data->Bz_vfm || !data->Bx_model || !data->By_model || !data->Bz_model || !data->lt_eq || !data->lt ||
      !data->F || !data->q || !data->weights || !data->t_ns || !data->ts_ns || !data->r_ns || !data->index ||
      !data->theta_ns || !data->phi_ns || !data->qdlat_ns || !data->satdir || !data->ne ||
      !data->Bx_nec_ns || !data->By_nec_ns || !data->Bz_nec_ns ||
      !data->Bx_vfm_ns || !data->By_vfm_ns || !data->Bz_vfm_ns ||
      !data->Bx_model_ns || !data->By_model_ns || !data->Bz_model_ns || !data->F_ns || !data->q_ns || !data->lt_eq_ns ||
      !data->lt_ns || !data->dXdt_nec || !data->dYdt_nec || !data->dZdt_nec || !data->MLT)
    {
      magdata_free(data);
      fprintf(stderr, "magdata_realloc: error reallocating variables\n");
      return 0;
    }

  data->ntot = n;

  return data;
} /* magdata_realloc() */

void
magdata_free(magdata *data)
{
  if (data->t)
    free(data->t);

  if (data->ts)
    free(data->ts);

  if (data->r)
    free(data->r);

  if (data->theta)
    free(data->theta);

  if (data->phi)
    free(data->phi);

  if (data->qdlat)
    free(data->qdlat);

  if (data->MLT)
    free(data->MLT);

  if (data->Bx_nec)
    free(data->Bx_nec);

  if (data->By_nec)
    free(data->By_nec);

  if (data->Bz_nec)
    free(data->Bz_nec);

  if (data->Bx_vfm)
    free(data->Bx_vfm);

  if (data->By_vfm)
    free(data->By_vfm);

  if (data->Bz_vfm)
    free(data->Bz_vfm);

  if (data->Bx_model)
    free(data->Bx_model);

  if (data->By_model)
    free(data->By_model);

  if (data->Bz_model)
    free(data->Bz_model);

  if (data->dXdt_nec)
    free(data->dXdt_nec);

  if (data->dYdt_nec)
    free(data->dYdt_nec);

  if (data->dZdt_nec)
    free(data->dZdt_nec);

  if (data->ne)
    free(data->ne);

  if (data->t_ns)
    free(data->t_ns);

  if (data->ts_ns)
    free(data->ts_ns);

  if (data->r_ns)
    free(data->r_ns);

  if (data->theta_ns)
    free(data->theta_ns);

  if (data->phi_ns)
    free(data->phi_ns);

  if (data->qdlat_ns)
    free(data->qdlat_ns);

  if (data->Bx_nec_ns)
    free(data->Bx_nec_ns);

  if (data->By_nec_ns)
    free(data->By_nec_ns);

  if (data->Bz_nec_ns)
    free(data->Bz_nec_ns);

  if (data->Bx_vfm_ns)
    free(data->Bx_vfm_ns);

  if (data->By_vfm_ns)
    free(data->By_vfm_ns);

  if (data->Bz_vfm_ns)
    free(data->Bz_vfm_ns);

  if (data->Bx_model_ns)
    free(data->Bx_model_ns);

  if (data->By_model_ns)
    free(data->By_model_ns);

  if (data->Bz_model_ns)
    free(data->Bz_model_ns);

  if (data->F_ns)
    free(data->F_ns);

  if (data->q_ns)
    free(data->q_ns);

  if (data->lt_ns)
    free(data->lt_ns);

  if (data->lt_eq_ns)
    free(data->lt_eq_ns);

  if (data->F)
    free(data->F);

  if (data->q)
    free(data->q);

  if (data->lt)
    free(data->lt);

  if (data->lt_eq)
    free(data->lt_eq);

  if (data->weights)
    free(data->weights);

  if (data->flags)
    free(data->flags);

  if (data->index)
    free(data->index);

  if (data->satdir)
    free(data->satdir);

  if (data->weight_workspace_p)
    track_weight_free(data->weight_workspace_p);

  free(data);
}

int
magdata_set_euler(const size_t flags, magdata *data)
{
  data->euler_flags = flags;
  return GSL_SUCCESS;
}

int
magdata_datum_init(magdata_datum *datum)
{
  int s = 0;
  size_t i;

  datum->name[0] = '\0';
  datum->t = 0.0;
  datum->r = 0.0;
  datum->theta = 0.0;
  datum->phi = 0.0;
  datum->qdlat = 0.0;
  datum->MLT = 0.0;
  datum->F = 0.0;
  datum->ne = 0.0;
  datum->satdir = 0;

  datum->t_ns = 0.0;
  datum->r_ns = 0.0;
  datum->theta_ns = 0.0;
  datum->phi_ns = 0.0;
  datum->qdlat_ns = 0.0;
  datum->F_ns = 0.0;

  datum->flags = 0;

  for (i = 0; i < 3; ++i)
    {
      datum->B_nec[i] = 0.0;
      datum->B_vfm[i] = 0.0;
      datum->B_model[i] = 0.0;
      datum->dBdt_nec[i] = 0.0;
      datum->B_nec_ns[i] = 0.0;
      datum->B_vfm_ns[i] = 0.0;
      datum->B_model_ns[i] = 0.0;
    }

  for (i = 0; i < 4; ++i)
    {
      datum->q[i] = 0.0;
      datum->q_ns[i] = 0.0;
    }

  return s;
} /* magdata_datum_init() */

/*
magdata_add()
  Add data point to magdata data structure

Inputs: datum - point data
        data  - data structure

Return: success/error
*/

int
magdata_add(const magdata_datum *datum, magdata *data)
{
  int s = 0;
  const size_t n = data->n;
  double *q = &(data->q[4 * n]);
  double *q_ns = &(data->q_ns[4 * n]);
  size_t i;

  if (n >= data->ntot)
    {
      fprintf(stderr, "magdata_add: error: maximum records reached: %zu\n", n);
      return -1;
    }

  strncpy(data->name, datum->name, MAGDATA_NAME_LENGTH);
  data->t[n] = datum->t;
  data->ts[n] = 0.0; /* filled in later */
  data->r[n] = datum->r;
  data->theta[n] = datum->theta;
  data->phi[n] = wrappi(datum->phi);
  data->qdlat[n] = datum->qdlat;
  data->MLT[n] = datum->MLT;
  data->flags[n] = datum->flags;
  data->index[n] = 0; /* filled in later */
  data->Bx_nec[n] = datum->B_nec[0];
  data->By_nec[n] = datum->B_nec[1];
  data->Bz_nec[n] = datum->B_nec[2];
  data->Bx_vfm[n] = datum->B_vfm[0];
  data->By_vfm[n] = datum->B_vfm[1];
  data->Bz_vfm[n] = datum->B_vfm[2];
  data->Bx_model[n] = datum->B_model[0];
  data->By_model[n] = datum->B_model[1];
  data->Bz_model[n] = datum->B_model[2];
  data->F[n] = datum->F;
  data->dXdt_nec[n] = datum->dBdt_nec[0];
  data->dYdt_nec[n] = datum->dBdt_nec[1];
  data->dZdt_nec[n] = datum->dBdt_nec[2];
  data->ne[n] = datum->ne;
  data->satdir[n] = datum->satdir;
  data->weights[n] = 1.0; /* computed in magdata_calc() */
  data->lt[n] = datum->lt;
  data->lt_eq[n] = datum->lt_eq;

  for (i = 0; i < 4; ++i)
    {
      q[i] = datum->q[i];
      q_ns[i] = datum->q_ns[i];
    }

  /*
   * add along-track gradient information; if not available the fields
   * of datum should be 0
   */
  data->t_ns[n] = datum->t_ns;
  data->ts_ns[n] = 0.0; /* filled in later */
  data->r_ns[n] = datum->r_ns;
  data->theta_ns[n] = datum->theta_ns;
  data->phi_ns[n] = datum->phi_ns;
  data->qdlat_ns[n] = datum->qdlat_ns;
  data->F_ns[n] = datum->F_ns;
  data->Bx_nec_ns[n] = datum->B_nec_ns[0];
  data->By_nec_ns[n] = datum->B_nec_ns[1];
  data->Bz_nec_ns[n] = datum->B_nec_ns[2];
  data->Bx_vfm_ns[n] = datum->B_vfm_ns[0];
  data->By_vfm_ns[n] = datum->B_vfm_ns[1];
  data->Bz_vfm_ns[n] = datum->B_vfm_ns[2];
  data->Bx_model_ns[n] = datum->B_model_ns[0];
  data->By_model_ns[n] = datum->B_model_ns[1];
  data->Bz_model_ns[n] = datum->B_model_ns[2];
  data->lt_ns[n] = datum->lt_ns;
  data->lt_eq_ns[n] = datum->lt_eq_ns;

  /* update rmin/rmax */
  data->rmin = GSL_MIN(data->rmin, data->r[n]);
  data->rmax = GSL_MAX(data->rmax, data->r[n]);

  ++(data->n);

  return s;
} /* magdata_add() */

/*
magdata_init()
  After data has been added to magdata struct,
pass through data and initialize spatial weighting histogram.
*/

int
magdata_init(magdata *data)
{
  int s = 0;
  size_t i;

  track_weight_reset(data->weight_workspace_p);

  for (i = 0; i < data->n; ++i)
    {
      if (data->flags[i] & MAGDATA_FLG_DISCARD)
        continue;

      /*
       * treat each (X,Y,Z) and (DX,DY,DZ) measurement separately
       * for spatial weighting
       */

      if (data->flags[i] & MAGDATA_FLG_X)
        track_weight_add_data(data->theta[i], data->phi[i], data->weight_workspace_p);
      if (data->flags[i] & MAGDATA_FLG_Y)
        track_weight_add_data(data->theta[i], data->phi[i], data->weight_workspace_p);
      if (data->flags[i] & MAGDATA_FLG_Z)
        track_weight_add_data(data->theta[i], data->phi[i], data->weight_workspace_p);

      if (data->flags[i] & MAGDATA_FLG_DX_NS)
        track_weight_add_data(data->theta[i], data->phi[i], data->weight_workspace_p);
      if (data->flags[i] & MAGDATA_FLG_DY_NS)
        track_weight_add_data(data->theta[i], data->phi[i], data->weight_workspace_p);
      if (data->flags[i] & MAGDATA_FLG_DZ_NS)
        track_weight_add_data(data->theta[i], data->phi[i], data->weight_workspace_p);
    }

  return s;
} /* magdata_init() */

/*
magdata_unit_weights()
  Set all spatial weights to 1.0

Notes:
1) on output, data->weights is updated with unit weights
*/

int
magdata_unit_weights(magdata *data)
{
  int s = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    data->weights[i] = 1.0;

  return s;
} /* magdata_unit_weights() */

/*
magdata_calc()
  Calculate spatial weights of data

Notes:
1) on output, data->weights is updated with spatial weights
*/

int
magdata_calc(magdata *data)
{
  int s = 0;
  size_t i;

  /*
   * compute data weights based on density/area using previously
   * constructed histogram
   */
  s = track_weight_calc(data->weight_workspace_p);

  data->nx = 0;
  data->ny = 0;
  data->nz = 0;
  data->nf = 0;
  data->ndx = 0;
  data->ndy = 0;
  data->ndz = 0;
  data->ndf = 0;
  data->nvec = 0;
  data->ngrad = 0;
  data->nres = 0;

  /* compute individual weights for each data point */
  for (i = 0; i < data->n; ++i)
    {
      s += track_weight_get(data->phi[i], data->theta[i], &(data->weights[i]),
                            data->weight_workspace_p);

      if (data->flags[i] & MAGDATA_FLG_X)
        ++(data->nx);

      if (data->flags[i] & MAGDATA_FLG_Y)
        ++(data->ny);

      if (data->flags[i] & MAGDATA_FLG_Z)
        ++(data->nz);

      if (data->flags[i] & MAGDATA_FLG_F)
        ++(data->nf);

      if (data->flags[i] & MAGDATA_FLG_DX_NS)
        ++(data->ndx);

      if (data->flags[i] & MAGDATA_FLG_DY_NS)
        ++(data->ndy);

      if (data->flags[i] & MAGDATA_FLG_DZ_NS)
        ++(data->ndz);

      if (data->flags[i] & MAGDATA_FLG_DF_NS)
        ++(data->ndf);

      ++(data->nvec);

      if (data->flags[i] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS))
        ++(data->ngrad);
    }

  data->nres = data->nx + data->ny + data->nz + data->nf +
               data->ndx + data->ndy + data->ndz + data->ndf;

  return s;
} /* magdata_calc() */

/*
magdata_print()
  Output data points in ASCII format

Inputs: prefix - file prefix where to store data
        data   - data

Return: success/error
*/

int
magdata_print(const char *prefix, const magdata *data)
{
  int s = 0;
  const char *fmtstr = "%ld %6.3f %8.3f %10.4f %10.4f %10.4f %10.4f %2d %.2f\n";
  size_t i, j;
  FILE *fp[11];
  size_t n = 11; /* number of files to write */
  char buf[2048];

  sprintf(buf, "%s_X.dat", prefix);
  fp[0] = fopen(buf, "w");

  sprintf(buf, "%s_Y.dat", prefix);
  fp[1] = fopen(buf, "w");

  sprintf(buf, "%s_Z.dat", prefix);
  fp[2] = fopen(buf, "w");

  sprintf(buf, "%s_F.dat", prefix);
  fp[3] = fopen(buf, "w");

  sprintf(buf, "%s_DX_NS.dat", prefix);
  fp[4] = fopen(buf, "w");

  sprintf(buf, "%s_DY_NS.dat", prefix);
  fp[5] = fopen(buf, "w");

  sprintf(buf, "%s_DZ_NS.dat", prefix);
  fp[6] = fopen(buf, "w");

  sprintf(buf, "%s_DF_NS.dat", prefix);
  fp[7] = fopen(buf, "w");

  sprintf(buf, "%s_DX_EW.dat", prefix);
  fp[8] = fopen(buf, "w");

  sprintf(buf, "%s_DY_EW.dat", prefix);
  fp[9] = fopen(buf, "w");

  sprintf(buf, "%s_DZ_EW.dat", prefix);
  fp[10] = fopen(buf, "w");

  for (i = 0; i < n; ++i)
    {
      if (!fp[i])
        {
          fprintf(stderr, "magdata_print: unable to open file %zu: %s\n",
                  i, strerror(errno));
          return -1;
        }
    }

  for (i = 0; i < n; ++i)
    {
      j = 1;
      fprintf(fp[i], "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", j++);
      fprintf(fp[i], "# Field %zu: local time (hours)\n", j++);
      fprintf(fp[i], "# Field %zu: season (day of year)\n", j++);
      fprintf(fp[i], "# Field %zu: longitude (degrees)\n", j++);
      fprintf(fp[i], "# Field %zu: geocentric latitude (degrees)\n", j++);
      fprintf(fp[i], "# Field %zu: QD latitude (degrees)\n", j++);
      fprintf(fp[i], "# Field %zu: geocentric altitude (km)\n", j++);
      fprintf(fp[i], "# Field %zu: satellite direction (+/- 1)\n", j++);
    }

  fprintf(fp[0], "# Field %zu: X residual (nT)\n", j);
  fprintf(fp[1], "# Field %zu: Y residual (nT)\n", j);
  fprintf(fp[2], "# Field %zu: Z residual (nT)\n", j);
  fprintf(fp[3], "# Field %zu: F residual (nT)\n", j);
  fprintf(fp[4], "# Field %zu: DX N/S residual (nT)\n", j);
  fprintf(fp[5], "# Field %zu: DY N/S residual (nT)\n", j);
  fprintf(fp[6], "# Field %zu: DZ N/S residual (nT)\n", j);
  fprintf(fp[7], "# Field %zu: DF N/S residual (nT)\n", j);
  fprintf(fp[8], "# Field %zu: DX E/W residual (nT)\n", j);
  fprintf(fp[9], "# Field %zu: DY E/W residual (nT)\n", j);
  fprintf(fp[10], "# Field %zu: DZ E/W residual (nT)\n", j);

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double lt = data->lt[i];
      double doy = get_season(unix_time);
      double phi = wrap180(data->phi[i] * 180.0 / M_PI);
      double lat = 90.0 - data->theta[i] * 180.0 / M_PI;
      double alt = data->r[i] - data->R;
      double B[4], dB[4];

      if (MAGDATA_Discarded(data->flags[i]))
        continue;

      magdata_residual(i, B, data);
      magdata_residual_dB_ns(i, dB, data);

      if (MAGDATA_ExistX(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        {
          fprintf(fp[0], fmtstr, unix_time, lt, doy, phi, lat, data->qdlat[i], alt, data->satdir[i], B[0]);
        }

      if (MAGDATA_ExistY(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        {
          fprintf(fp[1], fmtstr, unix_time, lt, doy, phi, lat, data->qdlat[i], alt, data->satdir[i], B[1]);
        }

      if (MAGDATA_ExistZ(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        {
          fprintf(fp[2], fmtstr, unix_time, lt, doy, phi, lat, data->qdlat[i], alt, data->satdir[i], B[2]);
        }

      if (MAGDATA_ExistScalar(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        {
          fprintf(fp[3], fmtstr, unix_time, lt, doy, phi, lat, data->qdlat[i], alt, data->satdir[i], B[3]);
        }

      if (MAGDATA_ExistDX_NS(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        {
          fprintf(fp[4], fmtstr, unix_time, lt, doy, phi, lat, data->qdlat[i], alt, data->satdir[i], data->satdir[i] * dB[0]);
        }

      if (MAGDATA_ExistDY_NS(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        {
          fprintf(fp[5], fmtstr, unix_time, lt, doy, phi, lat, data->qdlat[i], alt, data->satdir[i], data->satdir[i] * dB[1]);
        }

      if (MAGDATA_ExistDZ_NS(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        {
          fprintf(fp[6], fmtstr, unix_time, lt, doy, phi, lat, data->qdlat[i], alt, data->satdir[i], data->satdir[i] * dB[2]);
        }

      if (MAGDATA_ExistDF_NS(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        {
          fprintf(fp[7], fmtstr, unix_time, lt, doy, phi, lat, data->qdlat[i], alt, data->satdir[i], data->satdir[i] * dB[3]);
        }

      /* separate individual tracks with newlines */
      if (i < data->n - 1 && data->flags[i + 1] & MAGDATA_FLG_TRACK_START)
        {
          for (j = 0; j < n; ++j)
            fprintf(fp[j], "\n\n");
        }
    }

  for (i = 0; i < n; ++i)
    fclose(fp[i]);

  return s;
} /* magdata_print() */

/*
magdata_map()
  Output lat/lon map of data coverage

Inputs: prefix - file prefix for where to store data map
        data   - lat/lon data

Return: success/error
*/

int
magdata_map(const char *prefix, const magdata *data)
{
  int s = 0;
  size_t i;
  FILE *fp[11];
  size_t n = 11; /* number of files to write */
  char buf[2048];

  sprintf(buf, "%s_X.dat", prefix);
  fp[0] = fopen(buf, "w");

  sprintf(buf, "%s_Y.dat", prefix);
  fp[1] = fopen(buf, "w");

  sprintf(buf, "%s_Z.dat", prefix);
  fp[2] = fopen(buf, "w");

  sprintf(buf, "%s_F.dat", prefix);
  fp[3] = fopen(buf, "w");

  sprintf(buf, "%s_DX_NS.dat", prefix);
  fp[4] = fopen(buf, "w");

  sprintf(buf, "%s_DY_NS.dat", prefix);
  fp[5] = fopen(buf, "w");

  sprintf(buf, "%s_DZ_NS.dat", prefix);
  fp[6] = fopen(buf, "w");

  sprintf(buf, "%s_DF_NS.dat", prefix);
  fp[7] = fopen(buf, "w");

  sprintf(buf, "%s_DX_EW.dat", prefix);
  fp[8] = fopen(buf, "w");

  sprintf(buf, "%s_DY_EW.dat", prefix);
  fp[9] = fopen(buf, "w");

  sprintf(buf, "%s_DZ_EW.dat", prefix);
  fp[10] = fopen(buf, "w");

  for (i = 0; i < n; ++i)
    {
      if (!fp[i])
        {
          fprintf(stderr, "magdata_map: unable to open file %zu: %s\n",
                  i, strerror(errno));
          return -1;
        }
    }

  fprintf(fp[0], "# Spatial distribution of X vector data\n");
  fprintf(fp[1], "# Spatial distribution of Y vector data\n");
  fprintf(fp[2], "# Spatial distribution of Z vector data\n");
  fprintf(fp[3], "# Spatial distribution of F scalar data\n");
  fprintf(fp[4], "# Spatial distribution of N/S gradient DX vector data\n");
  fprintf(fp[5], "# Spatial distribution of N/S gradient DY vector data\n");
  fprintf(fp[6], "# Spatial distribution of N/S gradient DZ vector data\n");
  fprintf(fp[7], "# Spatial distribution of N/S scalar DF vector data\n");
  fprintf(fp[8], "# Spatial distribution of E/W gradient DX vector data\n");
  fprintf(fp[9], "# Spatial distribution of E/W gradient DY vector data\n");
  fprintf(fp[10], "# Spatial distribution of E/W gradient DZ vector data\n");

  for (i = 0; i < n; ++i)
    {
      size_t j = 1;
      fprintf(fp[i], "# Field %zu: time (decimal year)\n", j++);
      fprintf(fp[i], "# Field %zu: longitude (degrees)\n", j++);
      fprintf(fp[i], "# Field %zu: geocentric latitude (degrees)\n", j++);
      fprintf(fp[i], "# Field %zu: QD latitude (degrees)\n", j++);
      fprintf(fp[i], "# Field %zu: geocentric radius (km)\n", j++);
    }

#if 0
  fprintf(fp, "# Field %zu: vector data used in Euler angle fitting (1 or 0)\n", i++);
#endif

  for (i = 0; i < data->n; ++i)
    {
      double t = satdata_epoch2year(data->t[i]);
      double phi = wrap180(data->phi[i] * 180.0 / M_PI);
      double lat = 90.0 - data->theta[i] * 180.0 / M_PI;

      if (MAGDATA_Discarded(data->flags[i]))
        continue;

      if (MAGDATA_ExistX(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[0], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistY(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[1], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistZ(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[2], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistScalar(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[3], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistDX_NS(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[4], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistDY_NS(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[5], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistDZ_NS(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[6], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistDF_NS(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[7], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistDX_EW(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[8], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistDY_EW(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[9], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);

      if (MAGDATA_ExistDZ_EW(data->flags[i]) && MAGDATA_FitMF(data->flags[i]))
        fprintf(fp[10], "%f %.4f %.4f %.4f %.4f\n", t, phi, lat, data->qdlat[i], data->r[i]);
    }

  for (i = 0; i < n; ++i)
    fclose(fp[i]);

  return s;
} /* magdata_map() */

/*
magdata_residual()
  Compute magnetic residual in NEC frame

Inputs: idx  - data index
        B    - (output) magnetic residual in NEC frame (nT)
               B[0] = X residual (nT)
               B[1] = Y residual (nT)
               B[2] = Z residual (nT)
               B[3] = F residual (nT)
        data - magdata

Return: success/error

Notes:
1) B = B_i - B_{model}(r_i,t_i)
2) F = F_i - |B_{model}(r_i,t_i)|
*/

int
magdata_residual(const size_t idx, double B[4], const magdata *data)
{
  int s = 0;
  double B_model[4];

  B_model[0] = data->Bx_model[idx];
  B_model[1] = data->By_model[idx];
  B_model[2] = data->Bz_model[idx];
  B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

  B[0] = data->Bx_nec[idx] - B_model[0];
  B[1] = data->By_nec[idx] - B_model[1];
  B[2] = data->Bz_nec[idx] - B_model[2];
  B[3] = data->F[idx] - B_model[3];

  return s;
} /* magdata_residual() */

/*
magdata_residual_ns()
  Compute magnetic residual of along-track point

Inputs: idx  - data index
        B    - (output) magnetic residual in NEC frame (nT)
               B[0] = X residual (nT)
               B[1] = Y residual (nT)
               B[2] = Z residual (nT)
               B[3] = F residual (nT)
        data - magdata

Return: success/error

Notes:
1) B = B_i - B_{model}(r_i,t_i)
2) F = F_i - |B_{model}(r_i,t_i)|
*/

int
magdata_residual_ns(const size_t idx, double B[4], const magdata *data)
{
  int s = 0;
  double B_model[4];

  B_model[0] = data->Bx_model_ns[idx];
  B_model[1] = data->By_model_ns[idx];
  B_model[2] = data->Bz_model_ns[idx];
  B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

  B[0] = data->Bx_nec_ns[idx] - B_model[0];
  B[1] = data->By_nec_ns[idx] - B_model[1];
  B[2] = data->Bz_nec_ns[idx] - B_model[2];
  B[3] = data->F_ns[idx] - B_model[3];

  return s;
} /* magdata_residual_ns() */

/*
magdata_residual_dB_ns()
  Compute along-track magnetic field difference

Inputs: idx  - data index
        dB   - (output) magnetic dB residual in NEC frame (nT)
               dB[0] = dX residual (nT)
               dB[1] = dY residual (nT)
               dB[2] = dZ residual (nT)
               dB[3] = dF residual (nT)
        data - magdata

Return: success/error

Notes:
1) dB = [B_2 - B_{model}^{2}(r_2,t_2)] -
        [B_1 - B_{model}^{1}(r_1,t_1)]
2) dF = [F_2 - |B_{model}^{2}(r_2,t_2)|] -
        [F_1 - |B_{model}^{1}(r_1,t_1)|]
*/

int
magdata_residual_dB_ns(const size_t idx, double dB[4], const magdata *data)
{
  int s = 0;
  double B[4], B_ns[4];
  size_t i;

  /* compute magnetic residuals at each point */
  magdata_residual(idx, B, data);
  magdata_residual_ns(idx, B_ns, data);

  for (i = 0; i < 4; ++i)
    dB[i] = B_ns[i] - B[i];

  return s;
} /* magdata_residual_dB_ns() */

/*
magdata_prior()
  Store a priori model (B_main + B_crust + ...) in an array

Inputs: idx     - data index
        B_prior - (output) magnetic prior model in NEC frame (nT)
                  B_prior[0] = X model (nT)
                  B_prior[1] = Y model (nT)
                  B_prior[2] = Z model (nT)
                  B_prior[3] = F model (nT)
        data    - magdata

Return: success/error
*/

int
magdata_prior(const size_t idx, double B_prior[4], const magdata *data)
{
  int s = 0;

  B_prior[0] = data->Bx_model[idx];
  B_prior[1] = data->By_model[idx];
  B_prior[2] = data->Bz_model[idx];
  B_prior[3] = gsl_hypot3(B_prior[0], B_prior[1], B_prior[2]);

  return s;
}

/*
magdata_prior_grad()
  Store a priori model (B_main + B_crust + ...) for gradient point in an array

Inputs: idx     - data index
        B_prior - (output) magnetic prior model for gradient point in NEC frame (nT)
                  B_prior[0] = X model (nT)
                  B_prior[1] = Y model (nT)
                  B_prior[2] = Z model (nT)
                  B_prior[3] = F model (nT)
        data    - magdata

Return: success/error
*/

int
magdata_prior_grad(const size_t idx, double B_prior[4], const magdata *data)
{
  int s = 0;

  B_prior[0] = data->Bx_model_ns[idx];
  B_prior[1] = data->By_model_ns[idx];
  B_prior[2] = data->Bz_model_ns[idx];
  B_prior[3] = gsl_hypot3(B_prior[0], B_prior[1], B_prior[2]);

  return s;
}

/*
magdata_t()
  Find first/last timestamp in magdata structure, ignoring
discarded data

Inputs: t0   - (output) timestamp of first non-flagged data point (CDF_EPOCH)
        t1   - (output) timestamp of last non-flagged data point (CDF_EPOCH)
        data - data

Return: success/error
*/

int
magdata_t(double *t0, double *t1, const magdata *data)
{
  int s = 0;
  size_t i;

  *t0 = -1.0;
  *t1 = -1.0;

  for (i = 0; i < data->n; ++i)
    {
      if (!(data->flags[i] & MAGDATA_FLG_DISCARD))
        {
          *t0 = data->t[i];
          break;
        }
    }

  for (i = data->n; i > 0 && i--; )
    {
      if (!(data->flags[i] & MAGDATA_FLG_DISCARD))
        {
          *t1 = data->t[i];
          break;
        }
    }

  if (*t0 < 0.0 || *t1 < 0.0)
    s = -1; /* didn't find it */

  return s;
} /* magdata_t() */

int
magdata_write(const char *filename, magdata *data)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "magdata_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fwrite(&(data->n), sizeof(size_t), 1, fp);
  fwrite(&(data->R), sizeof(double), 1, fp);
  fwrite(&(data->rmin), sizeof(double), 1, fp);
  fwrite(&(data->rmax), sizeof(double), 1, fp);
  fwrite(&(data->nvec), sizeof(size_t), 1, fp);
  fwrite(&(data->nres), sizeof(size_t), 1, fp);
  fwrite(&(data->global_flags), sizeof(size_t), 1, fp);
  fwrite(&(data->euler_flags), sizeof(size_t), 1, fp);

  fwrite(data->name, sizeof(char), MAGDATA_NAME_LENGTH, fp);
  fwrite(data->t, sizeof(double), data->n, fp);
  fwrite(data->r, sizeof(double), data->n, fp);
  fwrite(data->theta, sizeof(double), data->n, fp);
  fwrite(data->phi, sizeof(double), data->n, fp);
  fwrite(data->qdlat, sizeof(double), data->n, fp);
  fwrite(data->MLT, sizeof(double), data->n, fp);
  fwrite(data->F, sizeof(double), data->n, fp);
  fwrite(data->F_ns, sizeof(double), data->n, fp);

  /* don't write vector arrays for EMAG dataset */
  if (!(data->global_flags & MAGDATA_GLOBFLG_SCALAR_GRID))
    {
      fwrite(data->Bx_nec, sizeof(double), data->n, fp);
      fwrite(data->By_nec, sizeof(double), data->n, fp);
      fwrite(data->Bz_nec, sizeof(double), data->n, fp);
      fwrite(data->Bx_vfm, sizeof(double), data->n, fp);
      fwrite(data->By_vfm, sizeof(double), data->n, fp);
      fwrite(data->Bz_vfm, sizeof(double), data->n, fp);
      fwrite(data->Bx_model, sizeof(double), data->n, fp);
      fwrite(data->By_model, sizeof(double), data->n, fp);
      fwrite(data->Bz_model, sizeof(double), data->n, fp);
      fwrite(data->dXdt_nec, sizeof(double), data->n, fp);
      fwrite(data->dYdt_nec, sizeof(double), data->n, fp);
      fwrite(data->dZdt_nec, sizeof(double), data->n, fp);
      fwrite(data->q, sizeof(double), 4 * data->n, fp);
      fwrite(data->satdir, sizeof(int), data->n, fp);
      fwrite(data->lt, sizeof(double), data->n, fp);
      fwrite(data->lt_eq, sizeof(double), data->n, fp);

      fwrite(data->t_ns, sizeof(double), data->n, fp);
      fwrite(data->r_ns, sizeof(double), data->n, fp);
      fwrite(data->theta_ns, sizeof(double), data->n, fp);
      fwrite(data->phi_ns, sizeof(double), data->n, fp);
      fwrite(data->qdlat_ns, sizeof(double), data->n, fp);
      fwrite(data->Bx_nec_ns, sizeof(double), data->n, fp);
      fwrite(data->By_nec_ns, sizeof(double), data->n, fp);
      fwrite(data->Bz_nec_ns, sizeof(double), data->n, fp);
      fwrite(data->Bx_vfm_ns, sizeof(double), data->n, fp);
      fwrite(data->By_vfm_ns, sizeof(double), data->n, fp);
      fwrite(data->Bz_vfm_ns, sizeof(double), data->n, fp);
      fwrite(data->Bx_model_ns, sizeof(double), data->n, fp);
      fwrite(data->By_model_ns, sizeof(double), data->n, fp);
      fwrite(data->Bz_model_ns, sizeof(double), data->n, fp);
      fwrite(data->q_ns, sizeof(double), 4 * data->n, fp);
      fwrite(data->lt_ns, sizeof(double), data->n, fp);
      fwrite(data->lt_eq_ns, sizeof(double), data->n, fp);
    }

  fwrite(data->flags, sizeof(size_t), data->n, fp);
  fwrite(data->weights, sizeof(double), data->n, fp);

  fclose(fp);

  return s;
} /* magdata_write() */

/*
magdata_read()
  Read binary magdata file written by magdata_write()

Inputs: filename - data file
        data     - data will be appended to this data structure (or NULL)

Notes:
1) After reading, the routines magdata_init() and magdata_calc()
should be called to update the spatial weighting histogram, flag outliers
and update nvec and nres counts
*/

magdata *
magdata_read(const char *filename, magdata *data)
{
  FILE *fp;
  size_t n;
  size_t ntot = 0;
  double R;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "magdata_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  if (data)
    ntot = data->ntot;

  fread(&n, sizeof(size_t), 1, fp);
  fread(&R, sizeof(double), 1, fp);

  ntot += n;
  data = magdata_realloc(ntot, R, data);
  if (!data)
    return NULL;

  fread(&(data->rmin), sizeof(double), 1, fp);
  fread(&(data->rmax), sizeof(double), 1, fp);
  fread(&(data->nvec), sizeof(size_t), 1, fp);
  fread(&(data->nres), sizeof(size_t), 1, fp);
  fread(&(data->global_flags), sizeof(size_t), 1, fp);
  fread(&(data->euler_flags), sizeof(size_t), 1, fp);

  fread(&(data->name[0]), sizeof(char), MAGDATA_NAME_LENGTH, fp);
  fread(&(data->t[data->n]), sizeof(double), n, fp);
  fread(&(data->r[data->n]), sizeof(double), n, fp);
  fread(&(data->theta[data->n]), sizeof(double), n, fp);
  fread(&(data->phi[data->n]), sizeof(double), n, fp);
  fread(&(data->qdlat[data->n]), sizeof(double), n, fp);
  fread(&(data->MLT[data->n]), sizeof(double), n, fp);
  fread(&(data->F[data->n]), sizeof(double), n, fp);
  fread(&(data->F_ns[data->n]), sizeof(double), n, fp);

  if (!(data->global_flags & MAGDATA_GLOBFLG_SCALAR_GRID))
    {
      fread(&(data->Bx_nec[data->n]), sizeof(double), n, fp);
      fread(&(data->By_nec[data->n]), sizeof(double), n, fp);
      fread(&(data->Bz_nec[data->n]), sizeof(double), n, fp);
      fread(&(data->Bx_vfm[data->n]), sizeof(double), n, fp);
      fread(&(data->By_vfm[data->n]), sizeof(double), n, fp);
      fread(&(data->Bz_vfm[data->n]), sizeof(double), n, fp);
      fread(&(data->Bx_model[data->n]), sizeof(double), n, fp);
      fread(&(data->By_model[data->n]), sizeof(double), n, fp);
      fread(&(data->Bz_model[data->n]), sizeof(double), n, fp);
      fread(&(data->dXdt_nec[data->n]), sizeof(double), n, fp);
      fread(&(data->dYdt_nec[data->n]), sizeof(double), n, fp);
      fread(&(data->dZdt_nec[data->n]), sizeof(double), n, fp);
      fread(&(data->q[4 * data->n]), sizeof(double), 4 * n, fp);
      fread(&(data->satdir[data->n]), sizeof(int), n, fp);
      fread(&(data->lt[data->n]), sizeof(double), n, fp);
      fread(&(data->lt_eq[data->n]), sizeof(double), n, fp);

      fread(&(data->t_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->r_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->theta_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->phi_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->qdlat_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->Bx_nec_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->By_nec_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->Bz_nec_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->Bx_vfm_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->By_vfm_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->Bz_vfm_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->Bx_model_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->By_model_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->Bz_model_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->q_ns[4 * data->n]), sizeof(double), 4 * n, fp);
      fread(&(data->lt_ns[data->n]), sizeof(double), n, fp);
      fread(&(data->lt_eq_ns[data->n]), sizeof(double), n, fp);
    }

  fread(&(data->flags[data->n]), sizeof(size_t), n, fp);
  fread(&(data->weights[data->n]), sizeof(double), n, fp);

  fclose(fp);

  data->n += n;

  return data;
} /* magdata_read() */

size_t
magdata_ndiscard(const magdata *data)
{
  size_t n = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      if (MAGDATA_Discarded(data->flags[i]))
        ++n;
    }

  return n;
} /* magdata_ndiscard() */

size_t
magdata_neuler(const magdata *data)
{
  size_t n = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      if (MAGDATA_FitEuler(data->flags[i]))
        ++n;
    }

  return n;
} /* magdata_neuler() */

int
magdata_clear(magdata *data)
{
  int s = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      data->flags[i] &= ~MAGDATA_FLG_DISCARD;
    }

  return s;
} /* magdata_clear() */

/*
magdata_flag_t()
  Flag all data points *not* inside [t0,t1] with
MAGDATA_FLG_DISCARD

Inputs: t0   - starting time (CDF_EPOCH)
        t1   - ending time (CDF_EPOCH)
        data - data
*/

int
magdata_flag_t(const double t0, const double t1,
               magdata *data)
{
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      if (data->t[i] < t0 || data->t[i] > t1)
        data->flags[i] |= MAGDATA_FLG_DISCARD;
    }

  return 0;
}

/*
magdata_flag_scalar()
  Flag all scalar data points with MAGDATA_FLG_DISCARD

Inputs: data - data
*/

int
magdata_flag_scalar(magdata *data)
{
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      if (MAGDATA_ExistScalar(data->flags[i]))
        data->flags[i] |= MAGDATA_FLG_DISCARD;
    }

  return 0;
}

/*
magdata_copy_track()
  Copy a single satellite track into magdata structure, discarding
bad data, and flagging when scalar/vector measurements are
available. For vector measurements, both VFM and NEC vectors
must be available or point is discarded.

  The start of the track is flagged with MAGDATA_FLG_TRACK_START
for later printing purposes.

Inputs: params    - parameters
        track_idx - track index
        data      - satellite data
        track_p   - track workspace
        mdata     - (output) where to store data
        ntype     - (output) counts of different data stored in mdata
                    ntype[0]   - number of scalar measurements
                    ntype[1]   - number of vector measurements (both VFM and NEC)
                    ntype[2]   - number of along-track scalar measurements
                    ntype[3]   - number of along-track vector measurements (both VFM and NEC)
                    ntype[4-5] - unused

Return: success/error

Notes:
1) ntype should be initialized by the calling function

2) mdata->global_flags is modified to include MAGDATA_GLOBFLG_SATELLITE
*/

int
magdata_copy_track(const magdata_params *params, const size_t track_idx,
                   const satdata_mag *data, const track_workspace *track_p,
                   magdata *mdata, size_t ntype[6])
{
  int s = 0;
  size_t i, k;
  track_data *tptr = &(track_p->tracks[track_idx]);
  const size_t start_idx = tptr->start_idx;
  const size_t end_idx = tptr->end_idx;
  const size_t grad_idx = (size_t) params->grad_dt_ns;
  const double grad_dt_min = params->grad_dt_ns - 2.0;
  const double grad_dt_max = params->grad_dt_ns + 2.0;
  magdata_datum datum;
  int flagged_start = 0;

  for (i = start_idx; i <= end_idx; ++i)
    {
      size_t flags = MAGDATA_FLG_F;

      /* ignore bad data */
      if (!SATDATA_AvailableData(data->flags[i]))
        continue;

      /* initialize to 0 */
      magdata_datum_init(&datum);

      if (!flagged_start)
        {
          flags |= MAGDATA_FLG_TRACK_START;
          flagged_start = 1;
        }

      /*
       * here there is at minimum a scalar measurement available; check
       * for vector measurement (both VFM and NEC)
       */
      if (SATDATA_ExistVector(data->flags[i]))
        {
          flags |= MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z;

          datum.B_nec[0] = SATDATA_VEC_X(data->B, i);
          datum.B_nec[1] = SATDATA_VEC_Y(data->B, i);
          datum.B_nec[2] = SATDATA_VEC_Z(data->B, i);
          datum.B_vfm[0] = SATDATA_VEC_X(data->B_VFM, i);
          datum.B_vfm[1] = SATDATA_VEC_Y(data->B_VFM, i);
          datum.B_vfm[2] = SATDATA_VEC_Z(data->B_VFM, i);
        }

      datum.F = data->F[i];

      for (k = 0; k < 4; ++k)
        datum.q[k] = data->q[4 * i + k];

      if (params->model_main)
        {
          datum.B_model[0] += SATDATA_VEC_X(data->B_main, i);
          datum.B_model[1] += SATDATA_VEC_Y(data->B_main, i);
          datum.B_model[2] += SATDATA_VEC_Z(data->B_main, i);
        }

      if (params->model_crust)
        {
          datum.B_model[0] += SATDATA_VEC_X(data->B_crust, i);
          datum.B_model[1] += SATDATA_VEC_Y(data->B_crust, i);
          datum.B_model[2] += SATDATA_VEC_Z(data->B_crust, i);
        }

      if (params->model_ext)
        {
          datum.B_model[0] += SATDATA_VEC_X(data->B_ext, i);
          datum.B_model[1] += SATDATA_VEC_Y(data->B_ext, i);
          datum.B_model[2] += SATDATA_VEC_Z(data->B_ext, i);
        }

      {
        size_t j = GSL_MIN(i + grad_idx, data->n - 1);
        double dt = (data->t[j] - data->t[i]) / 1000.0; /* in s */
        int status = 0;

        /* check if along-track point should be rejected */
        if (SATDATA_BadData(data->flags[j]) || (data->flags[j] & SATDATA_FLG_FILTER) ||
            (dt < grad_dt_min || dt > grad_dt_max))
          status = -1;

        if (status == 0)
          {
            /* set flag to indicate scalar gradient information available */
            flags |= MAGDATA_FLG_DF_NS;

            for (k = 0; k < 4; ++k)
              datum.q_ns[k] = data->q[4 * j + k];

            /* store along-track scalar measurement */
            datum.F_ns = data->F[j];

            /* check for vector measurement at along-track point */
            if (MAGDATA_ExistVector(flags) && SATDATA_ExistVector(data->flags[j]))
              {
                assert(data->flags[j] == 0 || data->flags[j] == SATDATA_FLG_DOWNSAMPLE ||
                       data->flags[j] == SATDATA_FLG_FILTER);

                flags |= MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS;

                datum.B_nec_ns[0] = SATDATA_VEC_X(data->B, j);
                datum.B_nec_ns[1] = SATDATA_VEC_Y(data->B, j);
                datum.B_nec_ns[2] = SATDATA_VEC_Z(data->B, j);
                datum.B_vfm_ns[0] = SATDATA_VEC_X(data->B_VFM, j);
                datum.B_vfm_ns[1] = SATDATA_VEC_Y(data->B_VFM, j);
                datum.B_vfm_ns[2] = SATDATA_VEC_Z(data->B_VFM, j);
              }

            if (params->model_main)
              {
                datum.B_model_ns[0] += SATDATA_VEC_X(data->B_main, j);
                datum.B_model_ns[1] += SATDATA_VEC_Y(data->B_main, j);
                datum.B_model_ns[2] += SATDATA_VEC_Z(data->B_main, j);
              }

            if (params->model_crust)
              {
                datum.B_model_ns[0] += SATDATA_VEC_X(data->B_crust, j);
                datum.B_model_ns[1] += SATDATA_VEC_Y(data->B_crust, j);
                datum.B_model_ns[2] += SATDATA_VEC_Z(data->B_crust, j);
              }

            if (params->model_ext)
              {
                datum.B_model_ns[0] += SATDATA_VEC_X(data->B_ext, j);
                datum.B_model_ns[1] += SATDATA_VEC_Y(data->B_ext, j);
                datum.B_model_ns[2] += SATDATA_VEC_Z(data->B_ext, j);
              }

            datum.t_ns = data->t[j];
            datum.r_ns = data->r[j];
            datum.theta_ns = M_PI / 2.0 - data->latitude[j] * M_PI / 180.0;
            datum.phi_ns = data->longitude[j] * M_PI / 180.0;
            datum.qdlat_ns = data->qdlat[j];
            datum.lt_ns = satdata_localtime(datum.t_ns, datum.phi_ns);
            datum.lt_eq_ns = tptr->lt_eq;
          }
      }

      datum.t = data->t[i];
      datum.r = data->r[i];
      datum.theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      datum.phi = data->longitude[i] * M_PI / 180.0;
      datum.qdlat = data->qdlat[i];
      datum.MLT = data->MLT[i];
      datum.ne = 0.0; /* filled in later */
      datum.satdir = satdata_mag_satdir(i, data);
      datum.flags = flags;
      datum.lt = satdata_localtime(datum.t, datum.phi);
      datum.lt_eq = tptr->lt_eq;

      s = magdata_add(&datum, mdata);
      if (s)
        return s;

      /* update counts */
      if (MAGDATA_ExistScalar(flags))
        ++(ntype[0]);
      if (MAGDATA_ExistVector(flags))
        ++(ntype[1]);
      if (flags & MAGDATA_FLG_DF_NS)
        ++(ntype[2]);
      if (flags & MAGDATA_FLG_DZ_NS)
        ++(ntype[3]);
    }

  mdata->global_flags |= MAGDATA_GLOBFLG_SATELLITE;

  return 0;
}

/*
magdata_copy_track_EW()
  Copy a single satellite track into magdata structure, discarding
bad data, and flagging when scalar/vector measurements are
available. For vector measurements, both VFM and NEC vectors
must be available or point is discarded.

  The start of the track is flagged with MAGDATA_FLG_TRACK_START
for later printing purposes.

Inputs: params    - parameters
        track_idx - track index
        data      - satellite data
        track_p   - track workspace
        data2     - satellite data for second satellite
        track_p2  - track workspace for second satellite
        mdata     - (output) where to store data
        ntype     - (output) counts of different data stored in mdata
                    ntype[0-3] - unused
                    ntype[4]   - number of E/W scalar measurements
                    ntype[5]   - number of E/W vector measurements (both VFM and NEC)

Return: success/error

Notes:
1) ntype should be initialized by the calling function

2) mdata->global_flags is modified to include MAGDATA_GLOBFLG_SATELLITE
*/

int
magdata_copy_track_EW(const magdata_params *params, const size_t track_idx,
                      const satdata_mag *data, const track_workspace *track_p,
                      const satdata_mag *data2, const track_workspace *track_p2,
                      magdata *mdata, size_t ntype[6])
{
  int s = 0;
  size_t i, k;
  track_data *tptr = &(track_p->tracks[track_idx]);
  track_data *tptr2;
  const size_t start_idx = tptr->start_idx;
  const size_t end_idx = tptr->end_idx;
  size_t track_idx2;
  magdata_datum datum;
  int flagged_start = 0;

  /* first locate the track of the second satellite close in time to the first */
  s = track_find_t(tptr->t_eq, &track_idx2, track_p2);

  if (s != GSL_SUCCESS)
    {
      /* track for second satellite not found, nothing to do */
      fprintf(stderr, "magdata_copy_track_EW: track for satellite 2 not found, timestamp: %f\n", tptr->t_eq);
      return 0;
    }

  tptr2 = &(track_p2->tracks[track_idx2]);

  if (tptr2->flags)
    {
      /*fprintf(stderr, "magdata_copy_track_EW: track for satellite 2 is flagged\n");*/
      return 0;
    }

  /* check time difference between tracks */
  {
    double dt = fabs(tptr->t_eq - tptr2->t_eq) * 1.0e-3; /* time difference in seconds */
    if (dt > params->grad_dt_ew)
      {
        /* time difference too large */
        fprintf(stderr, "magdata_copy_track_EW: track for satellite 2 not found due to dt [%.1f sec] [limit = %.1f sec]\n",\
                dt, params->grad_dt_ew);
        return 0;
      }
  }

  for (i = start_idx; i <= end_idx; ++i)
    {
      size_t flags = 0;
      size_t j;

      /* ignore bad data */
      if (!SATDATA_AvailableData(data->flags[i]))
        continue;

#if 1
      /* attempt to find a measurement for satellite 2 with approximately the
       * same latitude as satellite 1 */
      j = bsearch_double(data2->latitude, data->latitude[i],
                         tptr2->start_idx, tptr2->end_idx);
#else
      /* attempt to find a measurement for satellite 2 with approximately the
       * same time as satellite 1 */
      j = bsearch_double(data2->t, data->t[i],
                         tptr2->start_idx, tptr2->end_idx);
#endif

      /* now check that the latitude separation and time separation are within
       * allowed tolerances */
      if (fabs(data->latitude[i] - data2->latitude[j]) > params->grad_dlat_max)
        continue;
      else if (fabs(data->t[i] - data2->t[j]) > params->grad_dt_ew * 1000.0)
        continue;

      /*
       * check if satellite 2 measurement should be rejected - don't use
       * SATDATA_AvailableData() here because this point could have a downsample
       * flag applied which doesn't affect us using it as a gradient
       */
      if (SATDATA_BadData(data2->flags[j]) || (data2->flags[j] & SATDATA_FLG_FILTER))
        continue;

      /* initialize to 0 */
      magdata_datum_init(&datum);

      if (!flagged_start)
        {
          flags |= MAGDATA_FLG_TRACK_START;
          flagged_start = 1;
        }

      /*
       * here there is at minimum a scalar measurement available; check
       * for vector measurement (both VFM and NEC)
       */
      if (SATDATA_ExistVector(data->flags[i]))
        {
          datum.B_nec[0] = SATDATA_VEC_X(data->B, i);
          datum.B_nec[1] = SATDATA_VEC_Y(data->B, i);
          datum.B_nec[2] = SATDATA_VEC_Z(data->B, i);
          datum.B_vfm[0] = SATDATA_VEC_X(data->B_VFM, i);
          datum.B_vfm[1] = SATDATA_VEC_Y(data->B_VFM, i);
          datum.B_vfm[2] = SATDATA_VEC_Z(data->B_VFM, i);
        }

      datum.F = data->F[i];

      for (k = 0; k < 4; ++k)
        datum.q[k] = data->q[4 * i + k];

      if (params->model_main)
        {
          datum.B_model[0] += SATDATA_VEC_X(data->B_main, i);
          datum.B_model[1] += SATDATA_VEC_Y(data->B_main, i);
          datum.B_model[2] += SATDATA_VEC_Z(data->B_main, i);
        }

      if (params->model_crust)
        {
          datum.B_model[0] += SATDATA_VEC_X(data->B_crust, i);
          datum.B_model[1] += SATDATA_VEC_Y(data->B_crust, i);
          datum.B_model[2] += SATDATA_VEC_Z(data->B_crust, i);
        }

      if (params->model_ext)
        {
          datum.B_model[0] += SATDATA_VEC_X(data->B_ext, i);
          datum.B_model[1] += SATDATA_VEC_Y(data->B_ext, i);
          datum.B_model[2] += SATDATA_VEC_Z(data->B_ext, i);
        }

      /* set flag to indicate scalar gradient information available */
      flags |= MAGDATA_FLG_DF_EW;

      for (k = 0; k < 4; ++k)
        datum.q_ns[k] = data2->q[4 * j + k];

      /* store east-west scalar measurement */
      datum.F_ns = data2->F[j];

      /* if satellite 1 has vector measurement, check for satellite 2 vector measurement */
      if (SATDATA_ExistVector(data->flags[i]) && SATDATA_ExistVector(data2->flags[j]))
        {
          assert(data2->flags[j] == 0 || data2->flags[j] == SATDATA_FLG_DOWNSAMPLE);

          flags |= MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW;

          datum.B_nec_ns[0] = SATDATA_VEC_X(data2->B, j);
          datum.B_nec_ns[1] = SATDATA_VEC_Y(data2->B, j);
          datum.B_nec_ns[2] = SATDATA_VEC_Z(data2->B, j);
          datum.B_vfm_ns[0] = SATDATA_VEC_X(data2->B_VFM, j);
          datum.B_vfm_ns[1] = SATDATA_VEC_Y(data2->B_VFM, j);
          datum.B_vfm_ns[2] = SATDATA_VEC_Z(data2->B_VFM, j);
        }

      if (params->model_main)
        {
          datum.B_model_ns[0] += SATDATA_VEC_X(data2->B_main, j);
          datum.B_model_ns[1] += SATDATA_VEC_Y(data2->B_main, j);
          datum.B_model_ns[2] += SATDATA_VEC_Z(data2->B_main, j);
        }

      if (params->model_crust)
        {
          datum.B_model_ns[0] += SATDATA_VEC_X(data2->B_crust, j);
          datum.B_model_ns[1] += SATDATA_VEC_Y(data2->B_crust, j);
          datum.B_model_ns[2] += SATDATA_VEC_Z(data2->B_crust, j);
        }

      if (params->model_ext)
        {
          datum.B_model_ns[0] += SATDATA_VEC_X(data2->B_ext, j);
          datum.B_model_ns[1] += SATDATA_VEC_Y(data2->B_ext, j);
          datum.B_model_ns[2] += SATDATA_VEC_Z(data2->B_ext, j);
        }

      datum.t_ns = data2->t[j];
      datum.r_ns = data2->r[j];
      datum.theta_ns = M_PI / 2.0 - data2->latitude[j] * M_PI / 180.0;
      datum.phi_ns = data2->longitude[j] * M_PI / 180.0;
      datum.qdlat_ns = data2->qdlat[j];
      datum.lt_ns = satdata_localtime(datum.t_ns, datum.phi_ns);
      datum.lt_eq_ns = tptr2->lt_eq;

      datum.t = data->t[i];
      datum.r = data->r[i];
      datum.theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      datum.phi = data->longitude[i] * M_PI / 180.0;
      datum.qdlat = data->qdlat[i];
      datum.MLT = data->MLT[i];
      datum.ne = 0.0; /* filled in later */
      datum.satdir = satdata_mag_satdir(i, data);
      datum.flags = flags;
      datum.lt = satdata_localtime(datum.t, datum.phi);
      datum.lt_eq = tptr->lt_eq;

#if 0
      if (MAGDATA_ExistDZ_EW(flags))
        {
          printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                 datum.t,
                 datum.qdlat,
                 datum.qdlat_ns,
                 datum.B_nec[0],
                 datum.B_nec[1],
                 datum.B_nec[2],
                 datum.B_model[0],
                 datum.B_model[1],
                 datum.B_model[2],
                 datum.B_nec_ns[0],
                 datum.B_nec_ns[1],
                 datum.B_nec_ns[2],
                 datum.B_model_ns[0],
                 datum.B_model_ns[1],
                 datum.B_model_ns[2]);
        }
#endif

      s = magdata_add(&datum, mdata);
      if (s)
        return s;

      /* update counts */
      if (flags & MAGDATA_FLG_DF_EW)
        ++(ntype[4]);
      if (flags & MAGDATA_FLG_DZ_EW)
        ++(ntype[5]);
    }

  mdata->global_flags |= MAGDATA_GLOBFLG_SATELLITE;

  return 0;
}

/*
magdata_copy_station()
  Copy a single observatory station into magdata structure, discarding
bad data, and flagging when scalar/vector measurements are
available.

Inputs: params    - parameters
        station   - observatory station data
        mdata     - (output) where to store data
        ntype     - (output) counts of different data stored in mdata
                    ntype[0]   - number of scalar measurements
                    ntype[1]   - number of vector measurements (both VFM and NEC)
                    ntype[2]   - number of along-track scalar measurements
                    ntype[3]   - number of along-track vector measurements (both VFM and NEC)
                    ntype[4-5] - unused

Return: success/error

Notes:
1) ntype should be initialized by the calling function

2) mdata->global_flags is modified to include MAGDATA_GLOBFLG_OBSERVATORY
*/

int
magdata_copy_station(const magdata_params *params, const obsdata_station * station,
                     magdata *mdata, size_t ntype[6])
{
  magdata_datum datum;
  size_t i;

  /* initialize */
  magdata_datum_init(&datum);

  strncpy(datum.name, station->name, MAGDATA_NAME_LENGTH);

  for (i = 0; i < station->n_sv; ++i)
    {
      int s;
      size_t flags = 0;
      size_t idx;

      if (OBSDATA_ExistX(station->flags_sv[i]))
        flags |= MAGDATA_FLG_DXDT;

      if (OBSDATA_ExistY(station->flags_sv[i]))
        flags |= MAGDATA_FLG_DYDT;

      if (OBSDATA_ExistZ(station->flags_sv[i]))
        flags |= MAGDATA_FLG_DZDT;

      datum.t = station->t_sv[i];
      datum.r = station->radius;
      datum.theta = M_PI / 2.0 - station->latitude * M_PI / 180.0;
      datum.phi = station->longitude * M_PI / 180.0;
      datum.flags = flags;
      datum.lt = get_localtime(epoch2timet(datum.t), datum.phi);
      datum.dBdt_nec[0] = station->X_sv[i];
      datum.dBdt_nec[1] = station->Y_sv[i];
      datum.dBdt_nec[2] = station->Z_sv[i];

      idx = bsearch_double(station->t, datum.t, 0, station->n - 1);

      datum.qdlat = interp1d(station->t[idx], station->t[idx + 1],
                             station->qdlat[idx], station->qdlat[idx + 1],
                             datum.t);

      s = magdata_add(&datum, mdata);
      if (s)
        return s;
    }

  mdata->global_flags |= MAGDATA_GLOBFLG_OBSERVATORY;

  return 0;
}

/*
magdata_mag2sat()
  Convert magdata to satdata struct

Inputs: mdata - magdata

Return: satdata_mag struct
*/

satdata_mag *
magdata_mag2sat(const magdata *mdata)
{
  satdata_mag *data;
  size_t i, j;
  size_t idx = 0;

  data = satdata_mag_alloc(mdata->n);

  for (i = 0; i < mdata->n; ++i)
    {
      if (!MAGDATA_FitMF(mdata->flags[i]))
        continue;

      if (!MAGDATA_ExistZ(mdata->flags[i]))
        continue;

      data->t[idx] = mdata->t[i];
      data->r[idx] = mdata->r[i];
      data->latitude[idx] = 90.0 - mdata->theta[i] * 180.0 / M_PI;
      data->longitude[idx] = mdata->phi[i] * 180.0 / M_PI;
      data->qdlat[idx] = mdata->qdlat[i];
      data->MLT[idx] = mdata->MLT[i];
      SATDATA_VEC_X(data->B, idx) = mdata->Bx_nec[i];
      SATDATA_VEC_Y(data->B, idx) = mdata->By_nec[i];
      SATDATA_VEC_Z(data->B, idx) = mdata->Bz_nec[i];
      SATDATA_VEC_X(data->B_VFM, idx) = mdata->Bx_vfm[i];
      SATDATA_VEC_Y(data->B_VFM, idx) = mdata->By_vfm[i];
      SATDATA_VEC_Z(data->B_VFM, idx) = mdata->Bz_vfm[i];
      data->F[idx] = mdata->F[i];
      data->Flags_F[idx] = 0;
      data->Flags_B[idx] = 0;
      data->Flags_q[idx] = 0;

      for (j = 0; j < 4; ++j)
        data->q[4 * idx + j] = mdata->q[4 * i + j];

      ++idx;
    }

  data->n = idx;

  return data;
}

/*
magdata_replace_phi_LT()
  Replace longitude with corresponding LT angle:

phi_new = pi / 12 * (LT - LT0)

Inputs: lt0  - local time corresponding to 0 longitude (hours)
        data - magdata structure
*/

int
magdata_replace_phi_LT(const double lt0, magdata *data)
{
  int s = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      double lt;
      
      lt = data->lt[i];
      data->phi[i] = wrappi(M_PI / 12.0 * (lt - lt0));

      lt = data->lt_ns[i];
      data->phi_ns[i] = wrappi(M_PI / 12.0 * (lt - lt0));
    }

  return s;
}

/*
magdata_luhr()
  Calculate diamagnetic correction from Luhr et al, 2003

Inputs: ne     - electron density (cm^{-3})
        Te     - electron temperature (K)
        Ti     - ion temperature (K)
        B_nec  - magnetic field measurement in NEC (nT)
                 B_nec[0] = X (nT)
                 B_nec[1] = Y (nT)
                 B_nec[2] = Z (nT)
                 B_nec[3] = F (nT)
        b_luhr - (output) correction vector in NEC (nT)
                 b_luhr[0] = X (nT)
                 b_luhr[1] = Y (nT)
                 b_luhr[2] = Z (nT)
                 b_luhr[3] = F (nT)
*/

static int
magdata_luhr(const double ne, const double Te, const double Ti,
             const double B_nec[4], double b_luhr[4])
{
  int s = 0;
  const double F = B_nec[3];
  double b[3]; /* B / |B| */
  size_t i;

  /* compute unit vector in total field direction */
  for (i = 0; i < 3; ++i)
    b[i] = B_nec[i] / F;

  /* Luhr's formula */
  b_luhr[3] = (ne / F) * (Te + Ti) * MAGDATA_KB_MU0;

  /* the diamagnetic effect is in the minus field direction */
  for (i = 0; i < 3; ++i)
    b_luhr[i] = -b_luhr[3] * b[i];

  return s;
} /* magdata_luhr() */
