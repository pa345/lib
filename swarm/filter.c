/*
 * filter.c
 *
 * Filter Swarm data and flag
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <mainlib/ml_solarpos.h>
#include <mainlib/ml_track.h>

#define MAX_LATITUDE     (60.0)
#define MAX_ZENITH       (100.0)

/* flag any data within [lt_min,lt_max] */
size_t
swarm_filter_lt(const double lt_min, const double lt_max,
                satdata_mag *data)
{
  size_t i;
  size_t nlt = 0;
  solarpos_workspace *work_sp = solarpos_alloc();

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time;

      /* ignore already flagged data to improve speed */
      if (data->flags[i])
        continue;

      unix_time = satdata_epoch2timet(data->t[i]);

      /*
       * at high latitudes, calculate solar zenith angle for a better
       * measure of darkness/sunlight than local time
       */
      if (fabs(data->latitude[i]) > MAX_LATITUDE)
        {
          double lat = data->latitude[i] * M_PI / 180.0;
          double lon = wrappi(data->longitude[i] * M_PI / 180.0);
          double zenith;

          solarpos_calc_zenith(unix_time, lat, lon, &zenith, work_sp);
          zenith *= 180.0 / M_PI;
          assert(zenith >= 0.0);

          /* small zenith angle means sunlight so flag it */
          if (zenith < MAX_ZENITH)
            {
              ++nlt;
              data->flags[i] |= SATDATA_FLG_LT;
            }
        }
      else
        {
          double lt = get_localtime(unix_time, data->longitude[i] * M_PI / 180.0);

          if (lt >= lt_min && lt <= lt_max)
            {
              ++nlt;
              data->flags[i] |= SATDATA_FLG_LT;
            }
        }
    }

  solarpos_free(work_sp);

  return nlt;
}

static size_t
swarm_filter_instrument(satdata_mag *data)
{
  size_t i;
  size_t nflag = 0;

  for (i = 0; i < data->n; ++i)
    {
      if (data->Flags_F[i] >= 64)
        {
          data->flags[i] |= SATDATA_FLG_INSTERR;
          ++nflag;
        }
      else if (data->Flags_B[i] == 255)
        {
          data->flags[i] |= SATDATA_FLG_INSTERR;
          ++nflag;
        }
      else if (data->Flags_q[i] >= 31)
        {
          data->flags[i] |= SATDATA_FLG_INSTERR;
          ++nflag;
        }
    }

  return nflag;
}

static size_t
swarm_filter_asmvfm(const double thresh, satdata_mag *data)
{
  size_t i;
  size_t nflag = 0;

  for (i = 0; i < data->n; ++i)
    {
      double F_asm = data->F[i];
      double B[3], F_vfm;

      B[0] = SATDATA_VEC_X(data->B, i);
      B[1] = SATDATA_VEC_Y(data->B, i);
      B[2] = SATDATA_VEC_Z(data->B, i);
      F_vfm = gsl_hypot3(B[0], B[1], B[2]);

      if (fabs(F_asm - F_vfm) > thresh)
        {
          data->flags[i] |= SATDATA_FLG_INSTERR;
          ++nflag;
        }
    }

  return nflag;
}

static size_t
swarm_filter_qdlat(const double qdlat_min, const double qdlat_max,
                   satdata_mag *data)
{
  size_t i;
  size_t nflag = 0;

  for (i = 0; i < data->n; ++i)
    {
      if (fabs(data->qdlat[i]) < qdlat_min ||
          fabs(data->qdlat[i]) > qdlat_max)
        {
          data->flags[i] |= SATDATA_FLG_FILTER;
          ++nflag;
        }
    }

  return nflag;
} /* swarm_filter_qdlat() */

/*
swarm_filter()
  Flag Swarm data according to various criteria

Inputs: lt_min - minimum local time
        lt_max - maximum local time
        data   - satellite data

Notes:
1) Data inside the interval [lt_min,lt_max] is flagged. To
prevent local time flagging, set lt_min = lt_max
*/

int
swarm_filter(const double lt_min, const double lt_max,
             satdata_mag *data)
{
  int s = 0;
  size_t nflag;

  satdata_filter_wmm(1, data);

  fprintf(stderr, "swarm_filter: checking instrument flags...");
  nflag = swarm_filter_instrument(data);
  fprintf(stderr, "done (%zu data flagged)\n", nflag);

  if (lt_min < lt_max)
    {
      fprintf(stderr, "swarm_filter: flagging points inside LT window [%g,%g]...",
              lt_min, lt_max);
      nflag = swarm_filter_lt(lt_min, lt_max, data);
      fprintf(stderr, "done (%zu data flagged)\n", nflag);
    }

  fprintf(stderr, "swarm_filter: checking discrepencies between ASM/VFM...");
  nflag = swarm_filter_asmvfm(2.0, data);
  fprintf(stderr, "done (%zu data flagged)\n", nflag);
  
  nflag = satdata_nflagged(data);
  fprintf(stderr, "swarm_filter: total flagged points: %zu/%zu (%.1f%%) (%zu remaining)\n",
          nflag,
          data->n,
          (double)nflag / (double)data->n * 100.0,
          data->n - nflag);

  return s;
}
