/*
 * stage2a.c
 *
 * 1. Read DMSP file(s)
 * 2. Correct single point spikes with median filter
 * 3. Correct jumps
 * 4. Select quiet-time tracks (stage2_filter)
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
#include <gsl/gsl_deriv.h>
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

#include "eph.h"
#include "magcal.h"
#include "track.h"

#include "stage2_filter.c"
#include "stage2_jumps.c"
#include "stage2_quaternions.c"
#include "stage2_spikes.c"

#define WRITE_JUMP_DATA                   1

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
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

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
      fprintf(fp, "# Field %zu: angle AXY (degrees)\n", i++);
      fprintf(fp, "# Field %zu: angle AXZ (degrees)\n", i++);
      fprintf(fp, "# Field %zu: angle AYZ (degrees)\n", i++);
      return s;
    }

  fprintf(fp, "%ld %f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
          t,
          rms,
          gsl_vector_get(c, MAGCAL_IDX_SX),
          gsl_vector_get(c, MAGCAL_IDX_SY),
          gsl_vector_get(c, MAGCAL_IDX_SZ),
          gsl_vector_get(c, MAGCAL_IDX_OX),
          gsl_vector_get(c, MAGCAL_IDX_OY),
          gsl_vector_get(c, MAGCAL_IDX_OZ),
          gsl_vector_get(c, MAGCAL_IDX_AXY) * 180.0 / M_PI,
          gsl_vector_get(c, MAGCAL_IDX_AXZ) * 180.0 / M_PI,
          gsl_vector_get(c, MAGCAL_IDX_AYZ) * 180.0 / M_PI);

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
          double theta = M_PI / 2.0 - data->latitude[j] * M_PI / 180.0;
          double phi = data->longitude[j] * M_PI / 180.0;
          double *q = &(data->q[4*j]);
          double B_model[4], B_model_VFM[3], B_model_ell[3], r_ECEF[3];

          sph2ecef(data->r[j], theta, phi, r_ECEF);

          /* compute model vector */
          satdata_mag_model(j, B_model, data);

          /* rotate model vector to VFM frame */
          quat_apply_inverse(q, B_model, B_model_VFM);

          /* convert B_model into ellipsoid NEC */
          ellipsoid_nec2ell(r_ECEF, B_model, B_model_ell);

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

static int
stage2_vfm2nec(satdata_mag *data)
{
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      double *B_VFM = &(data->B_VFM[3*i]);
      double *B_NEC = &(data->B[3*i]);
      double *q = &(data->q[4*i]);

      quat_apply(q, B_VFM, B_NEC);
    }

  return 0;
}

int
preclean_proc(FILE *fp_spike, const size_t start_idx, const size_t end_idx, size_t nspikes[3],
              double min_spike[3], double max_spike[3], satdata_mag * data)
{
  int s = 0;
  const size_t nsamp = end_idx - start_idx + 1;
  const size_t spike_K = 11;                    /* window size for spike Hampel filter */
  const size_t min_samp = GSL_MAX(10, spike_K); /* minimum samples needed for spike/jump processing */
  const double nsigma[3] = { 8.0, 8.0, 5.0 };   /* tuning parameters for Hampel filter in each component */

  if (nsamp < min_samp)
    {
      /* insufficient number of samples for processing */
      return -1;
    }

  stage2_correct_spikes(0, fp_spike, spike_K, nsigma, start_idx, end_idx, nspikes, min_spike, max_spike, data);

  return s;
}

/*
preclean()
  Clean data prior to scalar calibration by detecting and
removing spikes and jumps. Also check for data gaps and
feed only continuous data streams into the spike/jump correction
routines
*/

int
preclean(const char *spike_file, satdata_mag * data)
{
  int s = 0;
  const time_t max_gap = 2;   /* maximum allowed gap in seconds */
  FILE *fp_spike = NULL;
  time_t tcur = satdata_epoch2timet(data->t[0]);
  size_t nspikes[3] = { 0, 0, 0 };
  double min_spike[3] = { 1.0e6, 1.0e6, 1.0e6 };
  double max_spike[3] = { -1.0e6, -1.0e6, -1.0e6 };
  size_t start_idx = 0;
  size_t i;

  fprintf(stderr, "\n");

  if (spike_file)
    {
      fp_spike = fopen(spike_file, "w");
      stage2_correct_spikes(1, fp_spike, 0, NULL, 0, 0, NULL, NULL, NULL, NULL);
    }

  for (i = 0; i < data->n - 1; ++i)
    {
      time_t tnext = satdata_epoch2timet(data->t[i + 1]);
      time_t tdiff = tnext - tcur;

      if (tdiff > max_gap)
        {
          preclean_proc(fp_spike, start_idx, i, nspikes, min_spike, max_spike, data);
          start_idx = i + 1;
        }

      tcur = tnext;
    }

  /* process end of data */
  preclean_proc(fp_spike, start_idx, data->n - 1, nspikes, min_spike, max_spike, data);

  fprintf(stderr, "\t preclean: minimum spike amplitudes: %.2f [nT] X, %.2f [nT] Y, %.2f [nT] Z\n",
          min_spike[0], min_spike[1], min_spike[2]);

  fprintf(stderr, "\t preclean: maximum spike amplitudes: %.2f [nT] X, %.2f [nT] Y, %.2f [nT] Z\n",
          max_spike[0], max_spike[1], max_spike[2]);

  fprintf(stderr, "\t preclean: total spikes found: %zu (X), %zu (Y), %zu (Z)\n",
          nspikes[0], nspikes[1], nspikes[2]);

  if (fp_spike)
    fclose(fp_spike);

  return s;
}

int
main(int argc, char *argv[])
{
  const char *track_file = "track_data.dat";
  const char *quat_file = "stage2_quat.dat";
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
  eph_data *eph = NULL;
  double period = -1.0; /* period in days for fitting scalar calibration parameters */
  double period_ms = -1.0; /* period in ms for fitting scalar calibration parameters */
  size_t nbins = 1;     /* number of time bins for fitting calibration parameters */
  double *t;            /* array of timestamps, size nbins */
  track_workspace *track_p;
  gsl_vector *coef = gsl_vector_alloc(MAGCAL_P);
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

      c = getopt_long(argc, argv, "d:i:o:b:p:r:t:", long_options, &option_index);
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

          default:
            break;
        }
    }

  if (!data)
    {
      fprintf(stderr, "Usage: %s <-i dmsp_index_file> <-b bowman_ephemeris_file> [-o output_file] [-p param_file] [-r residual_file] [-d data_file] [-t period]\n",
              argv[0]);
      exit(1);
    }

  if (outfile)
    fprintf(stderr, "output file = %s\n", outfile);

  /* determine number of bins for fitting calibration parameters */
  if (period > 0.0)
    {
      const double tmin = data->t[0];
      const double tmax = data->t[data->n - 1];
      const double dt = (tmax - tmin) / 8.64e7; /* convert to days */

      nbins = dt / period;
    }

  /* convert to ms */
  period_ms = period * 8.64e7;

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

  fprintf(stderr, "main: period for calibration:  %.1f [days]\n", period);
  fprintf(stderr, "main: number of temporal bins: %zu\n", nbins);

  track_p = track_alloc();

  fprintf(stderr, "main: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: precleaning data for spikes and jumps...");
  gettimeofday(&tv0, NULL);
  preclean(spike_file, data);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: fixing track jumps...");
  stage2_correct_jumps(jump_file, data);
  fprintf(stderr, "done\n");

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
      stage2_scalar_calibrate(scal_file, data, track_p, coef, &rms);

      fprintf(stderr, "main: applying calibration parameters to data...");
      magcal_apply(coef, data);
      fprintf(stderr, "done\n");

      if (fp_param)
        {
          time_t unix_time = satdata_epoch2timet(t[0]);
          print_parameters(fp_param, 0, coef, unix_time, rms);
        }
    }
  else
    {
      for (i = 0; i < nbins; ++i)
        {
          double tmin = t[i] - 0.5 * period_ms;
          double tmax = t[i] + 0.5 * period_ms;

          /* flag points outside of our time window */
          stage2_flag_time(tmin, tmax, data);

          /* calibrate points inside the time window [tmin,tmax] */
          stage2_scalar_calibrate(scal_file, data, track_p, coef, &rms);

          /* apply calibration to data inside time window */
          magcal_apply(coef, data);

          /* remove time flag */
          stage2_unflag_time(data);

          if (fp_param)
            {
              time_t unix_time = satdata_epoch2timet(t[i]);
              print_parameters(fp_param, 0, coef, unix_time, rms);
            }
        }
    }

#if 1
  fprintf(stderr, "main: correcting quaternions for satellite drift...");
  stage2_correct_quaternions(quat_file, data, track_p);
  fprintf(stderr, "done (data written to %s)\n", quat_file);
#endif

#if 0
  {
    print_data2("data2.txt", data, track_p);
    exit(1);
  }
#endif

  fprintf(stderr, "main: computing B_NEC...");
  stage2_vfm2nec(data);
  fprintf(stderr, "done\n");

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
  track_free(track_p);
  satdata_mag_free(data);
  eph_data_free(eph);
  free(t);

  if (fp_param)
    fclose(fp_param);

  return 0;
}
