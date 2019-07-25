/*
 * stage2a.c
 *
 * 1. Read DMSP file(s)
 * 2. Correct jumps
 * 3. Correct spikes
 * 4. Select quiet tracks
 * 5. Calculate scalar calibration parameters (scale factors, offsets, non-orthogonality angles)
 *
 * Usage: ./stage2a <-i dmsp_index_file> <-o dmsp_output_file> [-p parameter_file]
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
#include <gsl/gsl_blas.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics.h>

#include <satdata/satdata.h>
#include <indices/indices.h>
#include <common/common.h>
#include <common/quat.h>
#include <msynth/msynth.h>
#include <track/track.h>
#include <att/att_calc.h>

#include "eph.h"
#include "jump.h"
#include "magcal.h"

#include "stage2_align.c"
#include "stage2_calibrate.c"
#include "stage2_euler.c"
#include "stage2_filter.c"
#include "stage2_spikes.c"

#include "attitude.c"

#define WRITE_JUMP_DATA                   0

size_t
stage2_flag_time(const double tmin, const double tmax, satdata_mag *data)
{
  size_t i;
  size_t nflagged = 0;

  for (i = 0; i < data->n; ++i)
    {
      if (data->t[i] >= tmin && data->t[i] <= tmax)
        continue;

      data->flags[i] |= SATDATA_FLG_TIME;
      ++nflagged;
    }

  return nflagged;
}

int
stage2_unflag_time(satdata_mag *data)
{
  size_t i;

  for (i = 0; i < data->n; ++i)
    data->flags[i] &= ~SATDATA_FLG_TIME;

  return 0;
}

double
stage2_scalar_rms(const satdata_mag * data)
{
  size_t i;
  double rms = 0.0;
  double n = 0.0;

  for (i = 0; i < data->n; ++i)
    {
      double B_model[4], ri;

      if (fabs(data->qdlat[i]) > 50.0)
        continue;

      satdata_mag_model(i, B_model, data);

      ri = data->F[i] - B_model[3];

      rms += ri * ri;
      n += 1.0;
    }

  if (n > 0.0)
    rms = sqrt(rms / n);

  return rms;
}

/*
stage2_scalar_calibrate()
  Perform scalar calibration

Inputs: data    - satellite data
        track_p - track workspace
        c       - (output) scalar calibration parameters
        rms     - (output) scalar residual rms after calibration (nT)

Return: success/error
*/

int
stage2_scalar_calibrate(const char *res_file, satdata_mag * data, track_workspace * track_p,
                        gsl_vector * c, double *rms)
{
  int s = 0;
  size_t nflagged = satdata_nflagged(data);
  size_t n = data->n - nflagged;
  size_t i, j;
  magcal_workspace *magcal_p = magcal_alloc(n);
  struct timeval tv0, tv1;

  /* add unflagged data to magcal workspace */
  fprintf(stderr, "main: adding data for scalar calibration...");
  gettimeofday(&tv0, NULL);

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;

      if (tptr->flags)
        continue;

      for (j = start_idx; j <= end_idx; ++j)
        {
          double B_VFM[3], B_model[4];

          if (data->flags[j])
            continue;

          if (fabs(data->qdlat[j]) > 55.0)
            continue;

          B_VFM[0] = SATDATA_VEC_X(data->B_VFM, j);
          B_VFM[1] = SATDATA_VEC_Y(data->B_VFM, j);
          B_VFM[2] = SATDATA_VEC_Z(data->B_VFM, j);

          satdata_mag_model(j, B_model, data);

          magcal_add_datum(data->t[j], B_VFM, B_model[3], magcal_p);
        }
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, %zu data added)\n", time_diff(tv0, tv1), magcal_p->n);

  /* set initial values of calibration parameters */
  gsl_vector_set(c, MAGCAL_IDX_SX, 1.0);
  gsl_vector_set(c, MAGCAL_IDX_SY, 1.0);
  gsl_vector_set(c, MAGCAL_IDX_SZ, 1.0);
  gsl_vector_set(c, MAGCAL_IDX_OX, 0.0);
  gsl_vector_set(c, MAGCAL_IDX_OY, 0.0);
  gsl_vector_set(c, MAGCAL_IDX_OZ, 0.0);
  gsl_vector_set(c, MAGCAL_IDX_AXY, M_PI / 2.0);
  gsl_vector_set(c, MAGCAL_IDX_AXZ, M_PI / 2.0);
  gsl_vector_set(c, MAGCAL_IDX_AYZ, M_PI / 2.0);

  fprintf(stderr, "main: performing scalar calibration...");
  gettimeofday(&tv0, NULL);
  magcal_proc(c, magcal_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  *rms = magcal_rms(magcal_p);

  if (res_file)
    {
      fprintf(stderr, "main: writing scalar calibration residuals to %s...", res_file);
      magcal_print_residuals(res_file, c, magcal_p);
      fprintf(stderr, "done\n");
    }

  magcal_free(magcal_p);

  return s;
}

int
print_parameters(FILE *fp, const int header, const gsl_vector *c, const time_t t, const double rms)
{
  int s = 0;

  if (header)
    {
      size_t i;

      i = 1;
      fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      fprintf(fp, "# Field %zu: rms misfit (nT)\n", i++);
      fprintf(fp, "# Field %zu: scale factor X\n", i++);
      fprintf(fp, "# Field %zu: scale factor Y\n", i++);
      fprintf(fp, "# Field %zu: scale factor Z\n", i++);
      fprintf(fp, "# Field %zu: offset X\n", i++);
      fprintf(fp, "# Field %zu: offset Y\n", i++);
      fprintf(fp, "# Field %zu: offset Z\n", i++);
      fprintf(fp, "# Field %zu: angle U1 (degrees)\n", i++);
      fprintf(fp, "# Field %zu: angle U2 (degrees)\n", i++);
      fprintf(fp, "# Field %zu: angle U3 (degrees)\n", i++);
      return s;
    }

  fprintf(fp, "%ld %f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
          t,
          rms,
          gsl_vector_get(c, FLUXCAL_IDX_SX),
          gsl_vector_get(c, FLUXCAL_IDX_SY),
          gsl_vector_get(c, FLUXCAL_IDX_SZ),
          gsl_vector_get(c, FLUXCAL_IDX_OX),
          gsl_vector_get(c, FLUXCAL_IDX_OY),
          gsl_vector_get(c, FLUXCAL_IDX_OZ),
          gsl_vector_get(c, FLUXCAL_IDX_U1) * 180.0 / M_PI,
          gsl_vector_get(c, FLUXCAL_IDX_U2) * 180.0 / M_PI,
          gsl_vector_get(c, FLUXCAL_IDX_U3) * 180.0 / M_PI);

  fflush(fp);

  return s;
}

static int
print_data(const char *filename, const satdata_mag *data, const track_workspace *w)
{
  int s = 0;
  const size_t downsample = 120;
  FILE *fp;
  size_t i, j;
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);

  fp = fopen(filename, "w");

  i = 1;
  fprintf(fp, "# Field %zu: timestamp\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: scalar measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: VFM B_1 (nT)\n", i++);
  fprintf(fp, "# Field %zu: VFM B_2 (nT)\n", i++);
  fprintf(fp, "# Field %zu: VFM B_3 (nT)\n", i++);
  fprintf(fp, "# Field %zu: scalar model (nT)\n", i++);
  fprintf(fp, "# Field %zu: modeled VFM B_1 (nT)\n", i++);
  fprintf(fp, "# Field %zu: modeled VFM B_2 (nT)\n", i++);
  fprintf(fp, "# Field %zu: modeled VFM B_3 (nT)\n", i++);
  fprintf(fp, "# Field %zu: satellite direction (+1 north, -1 south)\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      size_t offset = (size_t) (gsl_rng_uniform(rng_p) * downsample);

      for (j = start_idx + offset; j <= end_idx; j += downsample)
        {
          double *q = &(data->q[4*j]);
          double B_model[4], B_model_VFM[3];

          /* compute model vector */
          satdata_mag_model(j, B_model, data);

          /* rotate model vector to VFM frame */
          quat_apply_inverse(q, B_model, B_model_VFM);

          fprintf(fp, "%ld %10.4f %10.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %d\n",
                  satdata_epoch2timet(data->t[j]),
                  data->qdlat[j],
                  data->F[j],
                  SATDATA_VEC_X(data->B_VFM, j),
                  SATDATA_VEC_Y(data->B_VFM, j),
                  SATDATA_VEC_Z(data->B_VFM, j),
                  B_model[3],
                  B_model_VFM[0],
                  B_model_VFM[1],
                  B_model_VFM[2],
                  satdata_satdir(j, data->n, data->latitude));
        }

      fprintf(fp, "\n\n");
    }

  fclose(fp);

  gsl_rng_free(rng_p);

  return s;
}

/*XXX*/
static int
print_data2(const char *filename, const satdata_mag *data, const track_workspace *w)
{
  int s = 0;
  const size_t downsample = 10;
  FILE *fp;
  size_t i, j;
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);

  fp = fopen(filename, "w");

  i = 1;
  fprintf(fp, "# Field %zu: timestamp\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: scalar measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC B_X (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC B_Y (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC B_Z (nT)\n", i++);
  fprintf(fp, "# Field %zu: scalar model (nT)\n", i++);
  fprintf(fp, "# Field %zu: model NEC B_X (nT)\n", i++);
  fprintf(fp, "# Field %zu: model NEC B_Y (nT)\n", i++);
  fprintf(fp, "# Field %zu: model NEC B_Z (nT)\n", i++);
  fprintf(fp, "# Field %zu: satellite direction (+1 north, -1 south)\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      size_t offset = (size_t) (gsl_rng_uniform(rng_p) * downsample);

      for (j = start_idx + offset; j <= end_idx; j += downsample)
        {
          double *q = &(data->q[4*j]);
          double *B_VFM = &(data->B_VFM[3*j]);
          double B_model[4], B_nec[3];

          /* compute model vector */
          satdata_mag_model(j, B_model, data);

          /* rotate model vector to VFM frame */
          quat_apply(q, B_VFM, B_nec);

          fprintf(fp, "%ld %10.4f %10.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %d\n",
                  satdata_epoch2timet(data->t[j]),
                  data->qdlat[j],
                  data->F[j],
                  B_nec[0],
                  B_nec[1],
                  B_nec[2],
                  B_model[3],
                  B_model[0],
                  B_model[1],
                  B_model[2],
                  satdata_satdir(j, data->n, data->latitude));
        }

      fprintf(fp, "\n\n");
    }

  fclose(fp);

  gsl_rng_free(rng_p);

  return s;
}

int
print_residuals(const char *filename, satdata_mag *data, track_workspace *w)
{
  int s = 0;
  FILE *fp;
  size_t i, j;

  fp = fopen(filename, "w");

  i = 1;
  fprintf(fp, "# Field %zu: timestamp\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: scalar measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: scalar model (nT)\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;

      if (tptr->flags)
        continue;

      for (j = start_idx; j <= end_idx; ++j)
        {
          double B_model[4];

          if (data->flags[j])
            continue;

          if (fabs(data->qdlat[j]) > 55.0)
            continue;

          satdata_mag_model(j, B_model, data);

          fprintf(fp, "%ld %.4f %f %f\n",
                  satdata_epoch2timet(data->t[j]),
                  data->qdlat[j],
                  data->F[j],
                  B_model[3]);
        }
    }

  fclose(fp);

  return s;
}

/*
preclean_jumps()
  Detect and remove jumps due to spacecraft fields in all 3 components in VFM frame
*/

int
preclean_jumps(const char *jump_file, satdata_mag * data, track_workspace *track_p)
{
  int s = 0;
  FILE *fp_jump = NULL;
  const time_t gap_threshold = 5;  /* maximum allowed gap in seconds */
  const size_t jump_K = 11; /* filter window size for jump detection */
  jump_workspace *jump_p;
  size_t njump[3] = { 0, 0, 0 };
  size_t start_idx = track_p->tracks[0].start_idx;
  size_t ngap = 0;            /* number of data gaps larger than gap_threshold */
  time_t max_gap = 0;         /* maximum data gap found */
  time_t tcur = satdata_epoch2timet(data->t[0]);
  struct timeval tv0, tv1;
  size_t i;

  fprintf(stderr, "\n");

  fprintf(stderr, "\t preclean_jumps: allocating workspace...");
  jump_p = jump_alloc(86400 * 90, jump_K);
  fprintf(stderr, "done\n");

  if (jump_file)
    {
      fp_jump = fopen(jump_file, "w");
      jump_proc(1, fp_jump, 0, 0, NULL, NULL, NULL);
    }

  /* flag positions of track starts in data */
  fprintf(stderr, "\t preclean_jumps: flagging track start positions...");

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      data->flags[tptr->start_idx] |= SATDATA_FLG_TRACK_START;
    }

  fprintf(stderr, "done\n");

  fprintf(stderr, "\t preclean_jumps: finding and correcting jumps...");
  gettimeofday(&tv0, NULL);

  for (i = start_idx; i < data->n - 1; ++i)
    {
      time_t tnext = satdata_epoch2timet(data->t[i + 1]);
      time_t tdiff = tnext - tcur;

      if (tdiff > gap_threshold)
        {
          jump_proc(0, fp_jump, start_idx, i, njump, data, jump_p);
          start_idx = i + 1;

          ++ngap;
          if (tdiff > max_gap)
            max_gap = tdiff;
        }

      tcur = tnext;
    }

  if (track_p->tracks[track_p->n - 1].end_idx > start_idx)
    jump_proc(0, fp_jump, start_idx, track_p->tracks[track_p->n - 1].end_idx, njump, data, jump_p);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "\t preclean_jumps: total data gaps found: %zu (maximum gap: %ld [sec])\n", ngap, max_gap);

  fprintf(stderr, "\t preclean_jumps: total jumps found and corrected: %zu (X), %zu (Y), %zu (Z)\n",
          njump[0], njump[1], njump[2]);

  if (fp_jump)
    fclose(fp_jump);

  jump_free(jump_p);

  return s;
}

/*
preclean_spikes()
  Detect and remove spikes due to spacecraft fields in all 3 components in VFM frame
*/

int
preclean_spikes(const char *spike_file, satdata_mag * data, track_workspace *track_p)
{
  int s = 0;
  FILE *fp_spike = NULL;
  const time_t gap_threshold = 5;  /* maximum allowed gap in seconds */
  const size_t spike_K = 15;  /* window size for impulse detection filter */
  const double nsigma[3] = { 5.0, 5.0, 5.0 }; /* tuning parameters for impulse rejection filter in each component */
  size_t nspikes[3] = { 0, 0, 0 };
  double min_spike[3] = { 1.0e6, 1.0e6, 1.0e6 };
  double max_spike[3] = { -1.0e6, -1.0e6, -1.0e6 };
  size_t start_idx = track_p->tracks[0].start_idx;
  size_t ngap = 0;            /* number of data gaps larger than gap_threshold */
  time_t max_gap = 0;         /* maximum data gap found */
  time_t tcur = satdata_epoch2timet(data->t[0]);
  struct timeval tv0, tv1;
  size_t i;

  fprintf(stderr, "\n");

  if (spike_file)
    {
      fp_spike = fopen(spike_file, "w");
      stage2_correct_spikes(1, fp_spike, 0, NULL, 0, 0, NULL, NULL, NULL, NULL);
    }

  /* flag positions of track starts in data */
  fprintf(stderr, "\t preclean_spikes: flagging track start positions...");

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      data->flags[tptr->start_idx] |= SATDATA_FLG_TRACK_START;
    }

  fprintf(stderr, "done\n");

  fprintf(stderr, "\t preclean_spikes: finding and correcting impulse spikes...");
  gettimeofday(&tv0, NULL);

  for (i = start_idx; i < data->n - 1; ++i)
    {
      time_t tnext = satdata_epoch2timet(data->t[i + 1]);
      time_t tdiff = tnext - tcur;

      if (tdiff > gap_threshold)
        {
          stage2_correct_spikes(0, fp_spike, spike_K, nsigma, start_idx, i, nspikes, min_spike, max_spike, data);
          start_idx = i + 1;

          ++ngap;
          if (tdiff > max_gap)
            max_gap = tdiff;
        }

      tcur = tnext;
    }

  if (track_p->tracks[track_p->n - 1].end_idx > start_idx)
    stage2_correct_spikes(0, fp_spike, spike_K, nsigma, start_idx, track_p->tracks[track_p->n - 1].end_idx, nspikes, min_spike, max_spike, data);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "\t preclean_spikes: total data gaps found: %zu (maximum gap: %ld [sec])\n", ngap, max_gap);

  fprintf(stderr, "\t preclean_spikes: total spikes found: %zu (X), %zu (Y), %zu (Z)\n",
          nspikes[0], nspikes[1], nspikes[2]);

  fprintf(stderr, "\t preclean_spikes: minimum spike amplitudes: %.2f [nT] X, %.2f [nT] Y, %.2f [nT] Z\n",
          min_spike[0], min_spike[1], min_spike[2]);

  fprintf(stderr, "\t preclean_spikes: maximum spike amplitudes: %.2f [nT] X, %.2f [nT] Y, %.2f [nT] Z\n",
          max_spike[0], max_spike[1], max_spike[2]);

  if (fp_spike)
    fclose(fp_spike);

  return s;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s <-i dmsp_index_file> <-b bowman_ephemeris_file> [-o output_file] [-p param_file] [-r residual_file] [-d data_file] [-t period]\n",
          argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -i dmsp_index_file     - Input DMSP index file\n");
  fprintf(stderr, "  -o output_file         - Output calibrated CDF file\n");
  fprintf(stderr, "  -p param_file          - Output ASCII file containing time series of calibration parameters\n");
  fprintf(stderr, "  -r residual_file       - Output ASCII file of residuals\n");
  fprintf(stderr, "  -d data_file           - Output ASCII file of data\n");
  fprintf(stderr, "  -t period              - Period in days of calibration parameter bins\n");
  fprintf(stderr, "  -h swarm_shc_file      - Replacement core field model in Swarm SHC format\n");
}

int
main(int argc, char *argv[])
{
#if WRITE_JUMP_DATA
  const char *scal_file = "stage2_scal.dat";
  const char *spike_file = "stage2_spikes.dat";
  const char *jump_file = "stage2_jumps.dat";
#else
  const char *scal_file = NULL;
  const char *spike_file = NULL;
  const char *jump_file = NULL;
#endif
  char *outfile = NULL;
  char *param_file = NULL;
  char *res_file = NULL;
  char *data_file = NULL;
  satdata_mag *data = NULL;
  satdata_mag *data2 = NULL;
  eph_data *eph = NULL;
  double rms0, rms1;    /* initial and final scalar rms for whole dataset (at low-latitudes) */
  double period = -1.0; /* period in days for fitting scalar calibration parameters */
  double period_ms = -1.0; /* period in ms for fitting scalar calibration parameters */
  size_t nbins = 1;     /* number of time bins for fitting calibration parameters */
  double *t;            /* array of timestamps of bin centers, size nbins */
  track_workspace *track_p;
  msynth_workspace *msynth_p = NULL;
  gsl_vector *coef = gsl_vector_alloc(FLUXCAL_P);
  att_calc_workspace * att_calc_p;
  gsl_vector *coef_align;
  FILE *fp_param = NULL;
  struct timeval tv0, tv1;
  double rms;
  size_t i;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "d:h:i:o:b:p:r:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_dmsp_read_idx(optarg, 1);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'd':
            data_file = optarg;
            break;

          case 'b':
            fprintf(stderr, "main: reading Bowman ephemerides from %s...", optarg);
            gettimeofday(&tv0, NULL);
            eph = eph_data_read_bowman(optarg);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu read, %g seconds)\n", eph->n, time_diff(tv0, tv1));

          case 'p':
            param_file = optarg;
            break;

          case 'r':
            res_file = optarg;
            break;

          case 't':
            period = atof(optarg);
            break;

          case 'h':
            msynth_p = msynth_swarm_read(optarg);
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  att_calc_p = att_calc_alloc(data->n, att_mrp);
  coef_align = gsl_vector_alloc(att_calc_p->p);

  if (outfile)
    fprintf(stderr, "output file = %s\n", outfile);

  /* determine number of bins for fitting calibration parameters */
  if (period > 0.0)
    {
      const double tmin = data->t[0];
      const double tmax = data->t[data->n - 1];
      const double dt = (tmax - tmin) / 8.64e7; /* convert to days */

      nbins = (size_t) (dt / period) + 1;
    }

  /*
   * this could be different than original period if time interval is not an exact multiple of
   * the original period
   */
  period_ms = (data->t[data->n - 1] - data->t[0]) / (double) nbins;

  /* build array of timestamps for each calibration bin */
  t = malloc(nbins * sizeof(double));

  if (nbins == 1)
    {
      t[0] = 0.5 * (data->t[data->n - 1] + data->t[0]);
    }
  else
    {
      for (i = 0; i < nbins; ++i)
        t[i] = data->t[0] + period_ms * (i + 0.5);
    }

  fprintf(stderr, "main: period for calibration:  %.1f [days]\n", period_ms / 8.64e7);
  fprintf(stderr, "main: number of temporal bins: %zu\n", nbins);

  if (msynth_p != NULL)
    {
      fprintf(stderr, "main: recalculating main field...");
      gettimeofday(&tv0, NULL);
      track_synth_core(data, msynth_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
    }

  /* calculate scalar rms of unmodified dataset */
  rms0 = stage2_scalar_rms(data);

  track_p = track_alloc();

  fprintf(stderr, "main: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

#if 1
  fprintf(stderr, "main: precleaning data for jumps...");
  gettimeofday(&tv0, NULL);
  preclean_jumps(jump_file, data, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, output jump file = %s)\n",
          time_diff(tv0, tv1), jump_file);

  fprintf(stderr, "main: precleaning data for impulse spikes...");
  gettimeofday(&tv0, NULL);
  preclean_spikes(spike_file, data, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, output spike file = %s)\n",
          time_diff(tv0, tv1), spike_file);
#endif

#if 1
  /* discard bad tracks according to rms test */
  fprintf(stderr, "main: filtering tracks with rms test...");
  gettimeofday(&tv0, NULL);
  stage2_filter(track_p, data);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (param_file)
    {
      fp_param = fopen(param_file, "w");
      print_parameters(fp_param, 1, NULL, 0, 0.0);
    }

  /* now do a scalar calibration separately for each time bin */
  if (nbins == 1)
    {
      stage2_calibrate(scal_file, data, track_p, coef, &rms);

      fprintf(stderr, "main: applying calibration parameters to data...");
      fluxcal_apply(coef, data);
      fprintf(stderr, "done\n");

      if (fp_param)
        {
          time_t unix_time = satdata_epoch2timet(t[0]);
          print_parameters(fp_param, 0, coef, unix_time, rms);
        }

      /* calculate alignment */
      stage2_align(data, track_p, coef_align, att_calc_p);

      /* apply alignment */
      stage2_align_vfm2nec(coef_align, data, att_calc_p);
    }
  else
    {
      for (i = 0; i < nbins; ++i)
        {
          double tmin = (i == 0) ? data->t[0] : (t[i] - 0.5 * period_ms);
          double tmax = (i == nbins - 1) ? data->t[data->n - 1] : (t[i] + 0.5 * period_ms);

          /* flag points outside of our time window */
          stage2_flag_time(tmin, tmax, data);

          /* calibrate points inside the time window [tmin,tmax] */
          stage2_calibrate(scal_file, data, track_p, coef, &rms);

          /* apply calibration to data inside time window */
          fluxcal_apply(coef, data);

          /* calculate alignment */
          stage2_align(data, track_p, coef_align, att_calc_p);

          /* apply alignment */
          stage2_align_vfm2nec(coef_align, data, att_calc_p);

          /* remove time flag */
          stage2_unflag_time(data);

          if (fp_param)
            {
              time_t unix_time = satdata_epoch2timet(t[i]);
              print_parameters(fp_param, 0, coef, unix_time, rms);
            }
        }
    }

  /* calculate scalar rms of final dataset */
  rms1 = stage2_scalar_rms(data);

  fprintf(stderr, "main: INITIAL SCALAR RMS = %.2f [nT]\n", rms0);
  fprintf(stderr, "main: FINAL SCALAR RMS   = %.2f [nT]\n", rms1);
#endif

#if 0

  fprintf(stderr, "main: performing attitude correction...");
  attitude_correct("attitude.txt", data, track_p);
  fprintf(stderr, "done\n");

#endif

  if (data_file)
    {
      fprintf(stderr, "main: printing data to %s...", data_file);
      print_data(data_file, data, track_p);
      fprintf(stderr, "done\n");
    }

  if (res_file)
    {
      fprintf(stderr, "main: printing residuals to %s...", res_file);
      print_residuals(res_file, data, track_p);
      fprintf(stderr, "done\n");
    }

  if (outfile)
    {
      fprintf(stderr, "main: writing data to %s...", outfile);
      satdata_dmsp_write(1, outfile, data);
      fprintf(stderr, "done\n");
    }

  gsl_vector_free(coef);
  gsl_vector_free(coef_align);
  track_free(track_p);
  satdata_mag_free(data);
  eph_data_free(eph);
  att_calc_free(att_calc_p);
  free(t);

  if (msynth_p)
    msynth_free(msynth_p);

  if (data2)
    satdata_mag_free(data2);

  if (fp_param)
    fclose(fp_param);

  return 0;
}
