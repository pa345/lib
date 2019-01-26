/*
 * jump.c
 *
 * Routines for detecting and correcting step jumps in time series
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_statistics.h>

#include <satdata/satdata.h>
#include <common/common.h>
#include <common/quat.h>
#include <track/track.h>

#include "jump.h"
#include "peak.h"

#include "yuen.c"

#define JUMP_DEBUG 0

typedef struct
{
  int idx;      /* index of possible jump in satdata_mag */
  double t;     /* timestamp of jump (CDF_EPOCH) */
  double delta; /* size of jump: x_{idx+1} - x_{idx} */
} jump_data;

typedef jump_data deque_type;

#include "deque.c"

static size_t jump_correct(const size_t comp, const size_t start_idx, const size_t end_idx, const gsl_vector *signal,
                           const medtest_stats_type * stats, const gsl_vector_int *ioutlier, gsl_vector_int *ijump, satdata_mag *data);
static int jump_apply_offset(const size_t comp, const double offset, const size_t start_idx, const size_t end_idx, satdata_mag * data);

/*
jump_alloc()
  Allocate a workspace for detecting jumps in time series

Inputs: n - maximum length of time series to process at one time
        K - window length for impulse filter
*/

jump_workspace *
jump_alloc(const size_t n, const size_t K)
{
  jump_workspace *w;
  size_t i;

  w = calloc(1, sizeof(jump_workspace));
  if (!w)
    {
      fprintf(stderr, "jump_alloc: calloc failed: %s\n",
              strerror(errno));
      return NULL;
    }

  w->peak_workspace_p = peak_alloc(n, 7, 35);
  w->impulse_workspace_p = gsl_filter_impulse_alloc(K);
  w->gaussian_workspace_p = gsl_filter_gaussian_alloc(K);

  for (i = 0; i < 3; ++i)
    {
      w->input[i] = gsl_vector_alloc(n);
      w->output_impulse[i] = gsl_vector_alloc(n);
      w->ioutlier[i] = gsl_vector_int_alloc(n);
      w->ijump[i] = gsl_vector_int_alloc(n);

      w->cumjump[i] = 0.0;
    }

  w->n = n;

  return w;
}

void
jump_free(jump_workspace * w)
{
  size_t i;

  if (w->peak_workspace_p)
    peak_free(w->peak_workspace_p);

  if (w->impulse_workspace_p)
    gsl_filter_impulse_free(w->impulse_workspace_p);

  if (w->gaussian_workspace_p)
    gsl_filter_gaussian_free(w->gaussian_workspace_p);

  for (i = 0; i < 3; ++i)
    {
      if (w->input[i])
        gsl_vector_free(w->input[i]);

      if (w->output_impulse[i])
        gsl_vector_free(w->output_impulse[i]);

      if (w->ioutlier[i])
        gsl_vector_int_free(w->ioutlier[i]);

      if (w->ijump[i])
        gsl_vector_int_free(w->ijump[i]);
    }

  free(w);
}

/*
jump_detect()
  Detect jumps (level shifts) in time series data

Inputs: min_jump_threshold - minimum vertical distance between adjacent window
                             medians to declare a jump
        x                  - input signal with jumps, length n
        stats              - (output) array of window statistics (median etc), length n
        ijump              - (output) boolean array of jump positions, length n
        w                  - workspace

Return: number of jumps detected
*/

size_t
jump_detect(const double min_jump_threshold, const gsl_vector * x,
            medtest_stats_type * stats, gsl_vector_int * ijump, jump_workspace * w)
{
  const int n = x->size;
  const int H = 10;
  const double alpha = 1.0e-6;
  const double nsigma = 10.0;
  size_t njump = 0;
  double window1[100], window2[100];
  int i;
  int nforward = 1;

  gsl_vector_int_set_zero(ijump);

  /* first pass: compute window statistics for all samples */
  for (i = 0; i < n - 1; ++i)
    {
      size_t wsize1 = gsl_movstat_fill(GSL_MOVSTAT_END_TRUNCATE, x, i, H, 0, window1);
      size_t wsize2 = gsl_movstat_fill(GSL_MOVSTAT_END_TRUNCATE, x, i + 1, 0, H, window2);

#if 0
      medtest(alpha, window1, 1, wsize1, window2, 1, wsize2, &stats[i]);
#else
      medtest2(alpha, window1, 1, wsize1, window2, 1, wsize2, &stats[i]);
#endif
    }

  /* second pass: look for jumps; start at sample H to ensure a full window at edges */
  for (i = H; i < n - H - 1; i += nforward)
    {
      double height = fabs(stats[i].median1 - stats[i].median2);

      /* ensure sigma is above a certain threshold */
      if (stats[i].se > 0.1 && height > min_jump_threshold && height > nsigma * stats[i].se)
        {
          /* jump detected */
          int j;

          for (j = i + 1; j < GSL_MIN(n - 1, i + H - 1); ++j)
            {
              double xj = gsl_vector_get(x, j);
              double d1 = fabs(xj - stats[i].median1);
              double d2 = fabs(xj - stats[i].median2);

              if (d2 < d1)
                break;
            }

          /* most likely jump location is j */
          gsl_vector_int_set(ijump, j, 1);

          /* skip forward to avoid detecting the same jump multiple times */
          nforward = 2*H;

          ++njump;
        }
      else
        {
          nforward = 1;
        }
    }

  return njump;
}

#if 1

/*
jump_proc()
  Search/correct square wave jumps in time series using
t-test based on medians

Inputs: header    - write file header (1 or 0)
        fp        - file pointer
        start_idx - process data from [start_idx,end_idx]
        end_idx   - process data from [start_idx,end_idx]
        njump     - (output) number of corrected jumps in each component
        data      - satellite data
        w         - jump workspace
*/

int
jump_proc(const int header, FILE *fp, const size_t start_idx, const size_t end_idx, size_t njump[3],
          satdata_mag * data, jump_workspace * w)
{
  const size_t n = end_idx - start_idx + 1;
  gsl_vector_view input[3];
  gsl_vector_int_view ioutlier[3];
  gsl_vector_int_view ijump[3];
  medtest_stats_type *stats[3];
  size_t i;

  if (header && fp)
    {
      i = 1;
      fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: original X VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original Y VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original Z VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: test statistic X (nT)\n", i++);
      fprintf(fp, "# Field %zu: test statistic Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: test statistic Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: X outlier detected (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: Y outlier detected (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: Z outlier detected (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: X jump detected (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: Y jump detected (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: Z jump detected (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: corrected X VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected Y VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected Z VFM residual (nT)\n", i++);
      return 0;
    }

  if (n > w->n)
    {
      fprintf(stderr, "jump_proc: error: data length too large: %zu [%zu]\n", n, w->n);
      return -1;
    }

  /* set up vector views */
  for (i = 0; i < 3; ++i)
    {
      stats[i] = malloc(n * sizeof(medtest_stats_type));
      input[i] = gsl_vector_subvector(w->input[i], 0, n);
      ioutlier[i] = gsl_vector_int_subvector(w->ioutlier[i], 0, n);
      ijump[i] = gsl_vector_int_subvector(w->ijump[i], 0, n);
    }

  /* compute vector residuals in VFM frame with a-priori model */
  for (i = 0; i < n; ++i)
    {
      size_t didx = i + start_idx;
      double *q = &(data->q[4 * didx]);
      double *B_VFM = &(data->B_VFM[3 * didx]);
      double B_model[4], B_model_VFM[3];

      satdata_mag_model(didx, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      /* compute residuals in VFM frame */
      gsl_vector_set(&input[0].vector, i, B_VFM[0] - B_model_VFM[0]);
      gsl_vector_set(&input[1].vector, i, B_VFM[1] - B_model_VFM[1]);
      gsl_vector_set(&input[2].vector, i, B_VFM[2] - B_model_VFM[2]);
    }

  for (i = 0; i < 3; ++i)
    {
      /* search for sudden jumps in this component */
      jump_detect(5.0, &input[i].vector, stats[i], &ioutlier[i].vector, w);

      /* correct jumps in this component */
      njump[i] += jump_correct(i, start_idx, end_idx, &input[i].vector, stats[i], &ioutlier[i].vector, &ijump[i].vector, data);
    }

  /* print output file */
  if (fp)
    {
      for (i = 0; i < n; ++i)
        {
          size_t didx = i + start_idx;
          double *q = &(data->q[4 * didx]);
          double *B_VFM = &(data->B_VFM[3 * didx]);
          double B_model[4], B_model_VFM[3];

          /* separate tracks with newlines */
          if (didx > 0 && data->flags[didx] & SATDATA_FLG_TRACK_START)
            fprintf(fp, "\n\n");

          satdata_mag_model(didx, B_model, data);

          /* rotate B_model into VFM frame */
          quat_apply_inverse(q, B_model, B_model_VFM);

          fprintf(fp, "%ld %8.4f %8.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %d %d %d %d %d %d %10.4f %10.4f %10.4f\n",
                  satdata_epoch2timet(data->t[didx]),
                  data->latitude[didx],
                  data->qdlat[didx],
                  gsl_vector_get(&input[0].vector, i),
                  gsl_vector_get(&input[1].vector, i),
                  gsl_vector_get(&input[2].vector, i),
                  stats[0][i].tstat,
                  stats[1][i].tstat,
                  stats[2][i].tstat,
                  gsl_vector_int_get(&ioutlier[0].vector, i),
                  gsl_vector_int_get(&ioutlier[1].vector, i),
                  gsl_vector_int_get(&ioutlier[2].vector, i),
                  gsl_vector_int_get(&ijump[0].vector, i),
                  gsl_vector_int_get(&ijump[1].vector, i),
                  gsl_vector_int_get(&ijump[2].vector, i),
                  B_VFM[0] - B_model_VFM[0],
                  B_VFM[1] - B_model_VFM[1],
                  B_VFM[2] - B_model_VFM[2]);
        }
    }

  for (i = 0; i < 3; ++i)
    free(stats[i]);

  return 0;
}

#elif 0

/* based on Gaussian smoothing filter */
int
jump_proc(const int header, FILE *fp, const size_t start_idx, const size_t end_idx,
          satdata_mag * data, jump_workspace * w)
{
  const size_t n = end_idx - start_idx + 1;
  const double sigma = 0.25;
  const double nsigma = 4.0;
  size_t njump[3];
  gsl_vector_view input[3], output_smooth[3];
  gsl_vector_int_view ioutlier[3];
  size_t i;

  if (header && fp)
    {
      i = 1;
      fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: original X VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original Y VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original Z VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: smoothed and differenced X (nT)\n", i++);
      fprintf(fp, "# Field %zu: smoothed and differenced Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: smoothed and differenced Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: X outlier detected (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: Y outlier detected (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: Z outlier detected (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: corrected X VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected Y VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected Z VFM residual (nT)\n", i++);
      return 0;
    }

  if (n > w->n)
    {
      fprintf(stderr, "jump_proc: error: data length too large: %zu [%zu]\n", n, w->n);
      return -1;
    }

  /* set up vector views */
  for (i = 0; i < 3; ++i)
    {
      input[i] = gsl_vector_subvector(w->input[i], 0, n);
      output_smooth[i] = gsl_vector_subvector(w->output_impulse[i], 0, n);
      ioutlier[i] = gsl_vector_int_subvector(w->ioutlier[i], 0, n);
    }

  /* compute vector residuals in VFM frame with a-priori model */
  for (i = 0; i < n; ++i)
    {
      size_t didx = i + start_idx;
      double *q = &(data->q[4 * didx]);
      double *B_VFM = &(data->B_VFM[3 * didx]);
      double B_model[4], B_model_VFM[3];

      satdata_mag_model(didx, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      /* compute residuals in VFM frame */
      gsl_vector_set(&input[0].vector, i, B_VFM[0] - B_model_VFM[0]);
      gsl_vector_set(&input[1].vector, i, B_VFM[1] - B_model_VFM[1]);
      gsl_vector_set(&input[2].vector, i, B_VFM[2] - B_model_VFM[2]);
    }

  for (i = 0; i < 3; ++i)
    {
      gsl_filter_gaussian(sigma, 1, &input[i].vector, &output_smooth[i].vector, w->gaussian_workspace_p);
    }

  {
    size_t npeak;

    peak_find(1.0, nsigma, &output_smooth[2].vector, &npeak, &ioutlier[2].vector, w->peak_workspace_p);
    fprintf(stderr, "npeak = %zu\n", npeak);
  }

  /* print output file */
  if (fp)
    {
      for (i = 0; i < n; ++i)
        {
          size_t didx = i + start_idx;
          double *q = &(data->q[4 * didx]);
          double *B_VFM = &(data->B_VFM[3 * didx]);
          double B_model[4], B_model_VFM[3];

          satdata_mag_model(didx, B_model, data);

          /* rotate B_model into VFM frame */
          quat_apply_inverse(q, B_model, B_model_VFM);

          fprintf(fp, "%ld %8.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %d %d %d %10.4f %10.4f %10.4f\n",
                  satdata_epoch2timet(data->t[didx]),
                  data->qdlat[didx],
                  gsl_vector_get(&input[0].vector, i),
                  gsl_vector_get(&input[1].vector, i),
                  gsl_vector_get(&input[2].vector, i),
                  gsl_vector_get(&output_smooth[0].vector, i),
                  gsl_vector_get(&output_smooth[1].vector, i),
                  gsl_vector_get(&output_smooth[2].vector, i),
                  gsl_vector_int_get(&ioutlier[0].vector, i),
                  gsl_vector_int_get(&ioutlier[1].vector, i),
                  gsl_vector_int_get(&ioutlier[2].vector, i),
                  B_VFM[0] - B_model_VFM[0],
                  B_VFM[1] - B_model_VFM[1],
                  B_VFM[2] - B_model_VFM[2]);
        }

      fprintf(fp, "\n\n");
    }

  return 0;
}

#endif

/*
jump_correct()
  Correct jumps which have already been detected

Inputs: comp      - component (0 = X, 1 = Y, 2 = Z)
        start_idx - starting index into data
        end_idx   - end index into data
        signal    - signal with jumps
        stats     - statistics for each data point
        ioutlier  - previously detected jump positions
        ijump     - locations of corrected jumps
        data      - data

Return: number of jumps corrected
*/

static size_t
jump_correct(const size_t comp, const size_t start_idx, const size_t end_idx, const gsl_vector *signal,
             const medtest_stats_type * stats, const gsl_vector_int *ioutlier, gsl_vector_int *ijump, satdata_mag *data)
{
  const double max_dt = 10.0 * 60.0;     /* maximum dt between jumps in seconds */
  const double eta_threshold = 0.16;     /* threshold for comparing similar jump sizes */
  const int n = (int) end_idx - (int) start_idx + 1;
  int i;
  size_t njump = 0;
  deque *deque_p = deque_alloc(1000);
  jump_data cur_jump;
  gsl_vector *x = gsl_vector_alloc(n);

  assert(n == (int) signal->size);
  putenv("TZ=GMT");

  gsl_vector_int_set_zero(ijump);

  /* make a copy of signal */
  gsl_vector_memcpy(x, signal);

  for (i = 0; i < n - 1; ++i)
    {
      int outlier = gsl_vector_int_get(ioutlier, i);
      int didx = start_idx + i;

      /* correction isn't reliable at high latitudes for X/Y */
      if (comp != 2 && fabs(data->qdlat[didx]) > 55.0)
        continue;

      if (outlier)
        {
          double t = data->t[didx];
          double delta = stats[i].median2 - stats[i].median1;
          double dt;
          double sum = 0.0;
          int nqueue, j, njumps;
          int jump_detected = 0;
          jump_data prev_jump;

          cur_jump.idx = didx;
          cur_jump.t = t;
          cur_jump.delta = delta;

#if JUMP_DEBUG
          if (delta >= 0.0)
            {
              fprintf(stderr, "Up %f\n", delta);
            }
          else
            {
              fprintf(stderr, "Down %f\n", delta);
            }
#endif

          /* remove past detected jumps which are outside allowed time window */
          while (!deque_is_empty(deque_p))
            {
              deque_peek_back(&prev_jump, deque_p);

              dt = (t - prev_jump.t) * 1.0e-3;
              if (dt > max_dt)
                deque_pop_back(deque_p);
              else
                break;
            }

          nqueue = deque_n(deque_p);
          for (j = 0; j < nqueue; ++j)
            {
              int idx = (j + deque_p->head) % deque_p->size;
              double eta;

              prev_jump = deque_p->array[idx];
              sum += prev_jump.delta;

              eta = fabs((fabs(delta) - fabs(sum)) / fabs(delta));

              if ((GSL_SIGN(sum) != GSL_SIGN(delta)) && (eta < eta_threshold || fabs(delta + sum) < 2.0))
                {
                  jump_detected = 1;
                  njumps = j + 1;
#if JUMP_DEBUG
                  {
                    time_t unix_time = satdata_epoch2timet(0.5 * (prev_jump.t + t));
                    fprintf(stderr, "FOUND JUMP: last %d samples match current sample: %s", njumps, ctime(&unix_time));
                  }
#endif
                  break;
                }
            }

          if (jump_detected)
            {
              int idx = didx;

              /*
               * a full sequence of jumps has been detected leading back to the original signal baseline;
               * now correct the jump samples
               *
               * njumps : number of elements at front of queue in the jump sequence
               */
              sum = delta;
              for (j = 0; j < njumps; ++j)
                {
                  deque_peek_front(&prev_jump, deque_p);

                  jump_apply_offset(comp, sum, prev_jump.idx + 1, idx, data);

                  gsl_vector_int_set(ijump, prev_jump.idx + 1 - start_idx, 1);
                  gsl_vector_int_set(ijump, idx - start_idx, 1);

                  sum += prev_jump.delta;
                  idx = prev_jump.idx;

                  deque_pop_front(deque_p);

                  ++njump;
                }
            }
          else
            {
              /* add this jump to front of queue */
              deque_push_front(cur_jump, deque_p);
            }
        }
    }

  deque_free(deque_p);
  gsl_vector_free(x);

  return njump;
}

/* add offset to component between [start_idx,end_idx] */
static int
jump_apply_offset(const size_t comp, const double offset, const size_t start_idx, const size_t end_idx, satdata_mag * data)
{
  size_t i;

  for (i = start_idx; i <= end_idx; ++i)
    {
      double *B_VFM = &(data->B_VFM[3 * i]);
      B_VFM[comp] += offset;
    }

  return 0;
}
