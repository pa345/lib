/*
 * track_synth.c
 *
 * Synthesize main, crustal, external field values along
 * satellite track
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#include <indices/indices.h>
#include <satdata/satdata.h>
#include <gsl/gsl_math.h>

#include <apex/apex.h>
#include <common/common.h>
#include <magcoord/magcoord.h>
#include <msynth/msynth.h>
#include <pomme/pomme.h>

#include "track.h"

/*
track_synth_int()
  Synthesize main and crustal fields along satellite track, with
parallelized implementation

Inputs: data           - satellite data output
        msynth_core_p  - msynth core workspace (degrees 1 to 15)
        msynth_crust_p - msynth crustal field workspace
*/

int
track_synth_int(satdata_mag *data, msynth_workspace *msynth_core_p, msynth_workspace *msynth_crust_p)
{
  int s = 0;
  size_t i;
  const size_t max_threads = (size_t) omp_get_max_threads();
  msynth_workspace **crust_p = malloc(max_threads * sizeof(msynth_workspace *));
  msynth_workspace **core_p = malloc(max_threads * sizeof(msynth_workspace *));

  for (i = 0; i < max_threads; ++i)
    {
      core_p[i] = msynth_copy(msynth_core_p);
      crust_p[i] = msynth_copy(msynth_crust_p);

      msynth_set(1, 15, core_p[i]);
      msynth_set(16, msynth_crust_p->eval_nmax, crust_p[i]);
    }

#pragma omp parallel for private(i)
  for (i = 0; i < data->n; ++i)
    {
      int thread_id = omp_get_thread_num();
      double tyr = satdata_epoch2year(data->t[i]);
      double r = data->r[i];
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double B_core[4], B_crust[4];

      /* compute internal and external fields */
      msynth_eval(tyr, r, theta, phi, B_core, core_p[thread_id]);
      msynth_eval(tyr, r, theta, phi, B_crust, crust_p[thread_id]);

      /* store vector core field */
      SATDATA_VEC_X(data->B_main, i) = B_core[0];
      SATDATA_VEC_Y(data->B_main, i) = B_core[1];
      SATDATA_VEC_Z(data->B_main, i) = B_core[2];

      /* store vector crustal field */
      SATDATA_VEC_X(data->B_crust, i) = B_crust[0];
      SATDATA_VEC_Y(data->B_crust, i) = B_crust[1];
      SATDATA_VEC_Z(data->B_crust, i) = B_crust[2];
    }

  for (i = 0; i < max_threads; ++i)
    {
      msynth_free(core_p[i]);
      msynth_free(crust_p[i]);
    }

  free(crust_p);
  free(core_p);

  return s;
}

/*
track_synth_core()
  Synthesize main field along satellite track, with
parallelized implementation

Inputs: data      - satellite data output
        msynth_p  - msynth core workspace (degrees 1 to 15)
*/

int
track_synth_core(satdata_mag *data, msynth_workspace *msynth_p)
{
  int s = 0;
  const size_t nmax = msynth_p->eval_nmax;
  size_t i;
  const size_t max_threads = (size_t) omp_get_max_threads();
  msynth_workspace **core_p = malloc(max_threads * sizeof(msynth_workspace *));

  for (i = 0; i < max_threads; ++i)
    {
      core_p[i] = msynth_copy(msynth_p);
      msynth_set(1, nmax, core_p[i]);
    }

#pragma omp parallel for private(i)
  for (i = 0; i < data->n; ++i)
    {
      int thread_id = omp_get_thread_num();
      double tyr = satdata_epoch2year(data->t[i]);
      double r = data->r[i];
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double B_core[4];

      /* compute core field */
      msynth_eval(tyr, r, theta, phi, B_core, core_p[thread_id]);

      /* store vector core field */
      SATDATA_VEC_X(data->B_main, i) = B_core[0];
      SATDATA_VEC_Y(data->B_main, i) = B_core[1];
      SATDATA_VEC_Z(data->B_main, i) = B_core[2];
    }

  for (i = 0; i < max_threads; ++i)
    msynth_free(core_p[i]);

  free(core_p);

  return s;
}

/*
track_synth_QD()
  Compute QD latitudes along track

Inputs: data - (input/output) satellite data output
*/

int
track_synth_QD(satdata_mag *data)
{
  int s = 0;
  size_t i;
  apex_workspace *apex_p = apex_alloc();

  /* compute QD latitudes (apex_transform is not thread safe so we can't parallelize this loop) */
  for (i = 0; i < data->n; ++i)
    {
      double tyr = satdata_epoch2year(data->t[i]);
      double r = data->r[i];
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double alon, alat, qdlat;

      /* compute QD latitude (apex_transform is not thread-safe) */
      apex_transform(tyr, theta, phi, r, &alon, &alat, &qdlat, NULL, NULL, NULL, apex_p);

      /* store QD latitude */
      data->qdlat[i] = qdlat;
    }

  apex_free(apex_p);

  return s;
}

/*
track_synth_MLT()
  Compute magnetic local time (MLT) along track

Inputs: data - (input/output) satellite data output
*/

int
track_synth_MLT(satdata_mag *data)
{
  int s = 0;
  size_t i;
  magcoord_workspace *magcoord_p = magcoord_alloc();

  /* compute MLT (apex_transform is not thread safe so we can't parallelize this loop) */
  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double r = data->r[i];
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;

      /* compute MLT (this routine is based on apex code and is not thread-safe) */
      data->MLT[i] = magcoord_MLT(unix_time, r, theta, phi, magcoord_p);
    }

  magcoord_free(magcoord_p);

  return s;
}

/*
track_synth_tilt()
  Compute dipole tilt angle along track

Inputs: data - (input/output) satellite data output
*/

int
track_synth_tilt(satdata_mag *data)
{
  int s = 0;
  const size_t max_threads = (size_t) omp_get_max_threads();
  magcoord_workspace **magcoord_p = malloc(max_threads * sizeof(magcoord_workspace *));
  size_t i;

  for (i = 0; i < max_threads; ++i)
    magcoord_p[i] = magcoord_alloc();

#pragma omp parallel for private(i)
  for (i = 0; i < data->n; ++i)
    {
      int thread_id = omp_get_thread_num();
      time_t unix_time = satdata_epoch2timet(data->t[i]);

      /* compute dipole tilt */
      data->tilt[i] = magcoord_tilt(unix_time, magcoord_p[thread_id]);
    }

  for (i = 0; i < max_threads; ++i)
    magcoord_free(magcoord_p[i]);

  free(magcoord_p);

  return s;
}

/*
track_synth_pomme()
  Synthesize POMME external field along satellite track.

Inputs: data - satellite data input/output

Notes:
1) On output, data->B_ext is filled with POMME values
*/

int
track_synth_pomme(satdata_mag *data)
{
  int s = 0;
  size_t i, j;
  const size_t max_threads = (size_t) omp_get_max_threads();
  pomme_workspace **ext_p = malloc(max_threads * sizeof(pomme_workspace *));
  estist_workspace *estist_workspace_p = estist_alloc(ESTIST_IDX_FILE);

  for (i = 0; i < max_threads; ++i)
    {
      ext_p[i] = pomme_alloc_default();
      pomme_set_radius(R_EARTH_KM, ext_p[i]);
    }

#pragma omp parallel for private(i)
  for (i = 0; i < data->n; ++i)
    {
      int thread_id = omp_get_thread_num();
      time_t t = satdata_epoch2timet(data->t[i]);
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double B_ext[4];
      double E_st = 0.0, I_st = 0.0;
      double IMF_By = 0.0, Em = 0.5;
      double f107 = 120.0;

      /*
       * satellite changed flight direction, new track; re-compute
       * indices for this track; this helps ensure POMME external
       * model is smooth/continuous over each track
       */
      s = pomme_get_indices(0, t, &E_st, &I_st, &IMF_By, &Em, &f107, ext_p[thread_id]);
      if (s)
        {
          fprintf(stderr, "track_synth_pomme: pomme_get_indices failed: %d\n", s);
          continue;
        }

      /* Est/Ist now interpolates so we can calculate it at each point */
      estist_get(t, &E_st, &I_st, estist_workspace_p);
      s = pomme_calc_ext_indices(theta, phi, t, data->altitude[i],
                                 E_st, I_st, IMF_By, Em, f107, B_ext, ext_p[thread_id]);
      if (s)
        {
          fprintf(stderr, "track_synth_pomme: pomme_calc_ext_indices failed: %d\n", s);
          continue;
        }

      /* convert to nT */
      for (j = 0; j < 3; ++j)
        B_ext[j] *= 1.0e9;

      /* store vector external field */
      SATDATA_VEC_X(data->B_ext, i) = B_ext[0];
      SATDATA_VEC_Y(data->B_ext, i) = B_ext[1];
      SATDATA_VEC_Z(data->B_ext, i) = B_ext[2];
    }

  for (i = 0; i < max_threads; ++i)
    pomme_free(ext_p[i]);

  free(ext_p);
  estist_free(estist_workspace_p);

  return s;
}
