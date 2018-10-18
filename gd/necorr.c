/*
 * necorr.c
 *
 * Correlate gravity magnetic signature with electron density peaks
 *
 * Pre-processing steps are:
 * 1. Instrument flags (recommended CHAMP flags except 1 star camera allowed)
 * 2. Track rms test
 * 3. Downsample by factor 15
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <flow/flow.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>

#include <common/bin3d.h>
#include <common/common.h>
#include <track/track.h>

#include "peak.h"

/* maximum number of peaks in Ne data */
#define MAX_PEAKS                  10

#define MIN_KP                     (0.0)
#define MAX_KP                     (2.0)

typedef struct
{
  int all;            /* print all tracks */
  double lt_min;      /* local time minimum (hours) */
  double lt_max;      /* local time maximum (hours) */
  double alt_min;     /* altitude minimum (km) */
  double alt_max;     /* altitude maximum (km) */
} preprocess_parameters;

satdata_mag *
read_swarm(const char *filename)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;

  fprintf(stderr, "read_swarm: reading %s...", optarg);
  gettimeofday(&tv0, NULL);
  data = satdata_swarm_read_idx(optarg, 0);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */

  fprintf(stderr, "read_swarm: filtering for instrument flags...");
  nflag = satdata_swarm_filter_instrument(1, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
} /* read_swarm() */

int
print_track_stats(const satdata_mag *data, const track_workspace *track_p)
{
  size_t nflagged = satdata_nflagged(data);
  size_t nleft = data->n - nflagged;
  size_t nflagged_track = track_nflagged(track_p);
  size_t nleft_track = track_p->n - nflagged_track;

  fprintf(stderr, "preprocess_data: total flagged data: %zu/%zu (%.1f%%)\n",
          nflagged, data->n, (double)nflagged / (double)data->n * 100.0);
  fprintf(stderr, "preprocess_data: total remaining data: %zu/%zu (%.1f%%)\n",
          nleft, data->n, (double)nleft / (double)data->n * 100.0);

  fprintf(stderr, "preprocess_data: total flagged tracks: %zu/%zu (%.1f%%)\n",
          nflagged_track, track_p->n, (double)nflagged_track / (double)track_p->n * 100.0);
  fprintf(stderr, "preprocess_data: total remaining tracks: %zu/%zu (%.1f%%)\n",
          nleft_track, track_p->n, (double)nleft_track / (double)track_p->n * 100.0);

  return 0;
} /* print_track_stats() */

/*
preprocess_data()

Inputs: params - preprocess parameters
          lt_min     - minimum local time
          lt_max     - maximum local time
        data   - satellite data

Return: pointer to sorted track workspace (should be freed by caller)
*/

track_workspace *
preprocess_data(const preprocess_parameters *params, satdata_mag *data)
{
  struct timeval tv0, tv1;
  track_workspace *track_p = track_alloc();

  fprintf(stderr, "preprocess_data: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (params->all)
    {
      print_track_stats(data, track_p);
      return track_p; /* no further filtering */
    }

  /* flag local time */
  if (params->lt_min >= 0.0 && params->lt_max >= 0.0)
    {
      size_t nlt = track_flag_lt(params->lt_min, params->lt_max, NULL, data, track_p);

      fprintf(stderr, "preprocess_data: flagged data outside LT window [%g,%g]: %zu/%zu (%.1f%%) tracks flagged)\n",
              params->lt_min, params->lt_max,
              nlt, track_p->n, (double)nlt / (double)track_p->n * 100.0);
    }

  /* flag altitude */
  {
    size_t nalt = track_flag_meanalt(params->alt_min, params->alt_max, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data due to altitude: %zu/%zu (%.1f%%) data flagged)\n",
            nalt, data->n, (double)nalt / (double)data->n * 100.0);
  }

  /* flag high kp data */
  {
    const double kp_min = MIN_KP;
    const double kp_max = MAX_KP;
    size_t nkp = track_flag_kp(kp_min, kp_max, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data outside kp window [%g,%g]: %zu/%zu (%.1f%%) data flagged)\n",
            kp_min, kp_max,
            nkp, data->n, (double)nkp / (double)data->n * 100.0);
  }

  /* print track statistics */
  {
    char *jstat_file = "track_jump_stats.dat";

    fprintf(stderr, "preprocess_data: printing jump track statistics to %s...", jstat_file);
    track_print_stats_flag(jstat_file, TRACK_FLG_JUMP, track_p);
    fprintf(stderr, "done\n");
  }

  print_track_stats(data, track_p);

  return track_p;
} /* preprocess_data() */

size_t
find_ne_peaks(const double qd_min, const double qd_max, const size_t width_points,
              const gsl_vector *qdlat, const gsl_vector *ne, size_t idx[],
              peak_workspace *peak_p)
{
  int s;
  const size_t smooth_window = width_points / 2;                  /* smoothing window size */
  const double minslope = 0.4 * pow((double) width_points, -2.0); /* minimum slope of first derivative */
  const double minheight = 8.0e5;                                 /* minimum peak height */
  size_t npeak = 0;

  /* initialize peak search */
  peak_init(smooth_window, qdlat, ne, peak_p);

  do
    {
      s = peak_find(1, minslope, minheight, qdlat, ne, peak_p);

      if (s == GSL_CONTINUE)
        {
          size_t pidx = peak_p->pidx; /* index into arrays of peak location */
          double center = gsl_vector_get(qdlat, pidx);

          /* test if peak center is in [qd_min,qd_max] */
          if (fabs(center) >= qd_min && fabs(center) <= qd_max)
            {
              /* found a peak, store Gaussian parameters */
              idx[npeak++] = pidx;

              if (npeak >= MAX_PEAKS)
                break;
            }
        }
    }
  while (s == GSL_CONTINUE);

  return npeak;
}

size_t
find_F_peaks(const double qd_min, const double qd_max, const size_t width_points,
             const gsl_vector *qdlat, const gsl_vector *F, size_t idx[],
             peak_workspace *peak_p)
{
  int s;
  const size_t smooth_window = width_points / 2;                  /* smoothing window size */
  const double minslope = 1.0 * pow((double) width_points, -2.0); /* minimum slope of first derivative */
  const double minheight = -1.0e6;                                /* minimum peak height */
  size_t npeak = 0;

  /* initialize peak search */
  peak_init(smooth_window, qdlat, F, peak_p);

  do
    {
      /*s = peak_find(-1, minslope, minheight, qdlat, F, peak_p);*/
      s = peak_find(0, minslope, minheight, qdlat, F, peak_p);

      if (s == GSL_CONTINUE)
        {
          size_t pidx = peak_p->pidx; /* index into arrays of peak location */
          double center = gsl_vector_get(qdlat, pidx);

          /* test if peak center is in [qd_min,qd_max] */
          if (center >= qd_min && center <= qd_max)
            {
              /* found a peak, store Gaussian parameters */
              idx[npeak++] = pidx;

              if (npeak >= MAX_PEAKS)
                break;
            }
        }
    }
  while (s == GSL_CONTINUE);

  return npeak;
}

int
correlate(const char *peak_file, const char *corr_file, const satdata_mag *data, const track_workspace *w)
{
  int s = 0;
  size_t i, j;
  const size_t width_points = 100;                                /* number of points in half-width of peak */
  const size_t fit_width = width_points * 1.5;                    /* gaussian fit window size */
  const size_t half_width = fit_width / 2;
  FILE *fp_peak, *fp_corr;

  size_t ne_idx[MAX_PEAKS];                                       /* indices for centers of Ne peaks */
  size_t F_idx[MAX_PEAKS];                                        /* indices for centers of F peaks */
  gsl_vector *ne_coef[MAX_PEAKS];                                 /* save gauss parameters for each peak */
  const size_t p = 3;

  const double qd_min = 4.0;                                      /* search for peaks only in [qd_min,qd_max] */
  const double qd_max = 30.0;

  fp_peak = fopen(peak_file, "w");
  fp_corr = fopen(corr_file, "w");

  i = 1;
  fprintf(fp_peak, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp_peak, "# Field %zu: electron density (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: scalar field residual (nT)\n", i++);
  fprintf(fp_peak, "# Field %zu: N_e first derivative (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: smoothed N_e first derivative (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: F first derivative (nT)\n", i++);
  fprintf(fp_peak, "# Field %zu: smoothed F first derivative (nT)\n", i++);
  fprintf(fp_peak, "# Field %zu: N_e data used for 1st peak fit (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: N_e data used for 2nd peak fit (cm^{-3})\n", i++);
  fprintf(fp_peak, "# Field %zu: scalar field peak 1 used for correlation (nT)\n", i++);
  fprintf(fp_peak, "# Field %zu: scalar field peak 2 used for correlation (nT)\n", i++);

  i = 1;
  fprintf(fp_corr, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp_corr, "# Field %zu: time (decimal year)\n", i++);
  fprintf(fp_corr, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp_corr, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp_corr, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp_corr, "# Field %zu: season (day of year)\n", i++);
  fprintf(fp_corr, "# Field %zu: N_e peak center in northern hemisphere (degrees QD latitude)\n", i++);
  fprintf(fp_corr, "# Field %zu: N_e peak height in northern hemisphere (cm^{-3})\n", i++);
  fprintf(fp_corr, "# Field %zu: N_e peak center in southern hemisphere (degrees QD latitude)\n", i++);
  fprintf(fp_corr, "# Field %zu: N_e peak height in southern hemisphere (cm^{-3})\n", i++);
  fprintf(fp_corr, "# Field %zu: F peak center in northern hemisphere (degrees QD latitude)\n", i++);
  fprintf(fp_corr, "# Field %zu: F peak height in northern hemisphere (nT)\n", i++);
  fprintf(fp_corr, "# Field %zu: F peak center in southern hemisphere (degrees QD latitude)\n", i++);
  fprintf(fp_corr, "# Field %zu: F peak height in southern hemisphere (nT)\n", i++);
  fprintf(fp_corr, "# Field %zu: corr(F,N_e) in northern hemisphere\n", i++);
  fprintf(fp_corr, "# Field %zu: corr(F,N_e) in southern hemisphere\n", i++);

  for (i = 0; i < MAX_PEAKS; ++i)
    ne_coef[i] = gsl_vector_alloc(p);

  for (i = 0; i < w->n; ++i)
    {
      time_t unix_time;
      track_data *tptr = &(w->tracks[i]);
      size_t sidx = tptr->start_idx;
      size_t n = tptr->n;
      peak_workspace *peak_ne = peak_alloc(n);
      peak_workspace *peak_F = peak_alloc(n);
      gsl_vector_view qd = gsl_vector_view_array(&(data->qdlat[sidx]), n);
      gsl_vector_view ne = gsl_vector_view_array(&(data->ne[sidx]), n);
      gsl_vector_view F = gsl_vector_view_array(tptr->Bf, n);
      size_t npeak_ne, npeak_F;
      double corr1, corr2;
      size_t ne_idx0[2], ne_idx1[2];
      size_t F_idx0[2], F_idx1[2];
      double ne_center1, ne_center2;

      /* discard flagged tracks */
      if (tptr->flags)
        continue;

      npeak_ne = find_ne_peaks(qd_min, qd_max, width_points, &qd.vector, &ne.vector, ne_idx, peak_ne);

      if (npeak_ne != 2)
        continue;

      /* store QD latitudes of peak centers */
      ne_center1 = gsl_vector_get(&qd.vector, ne_idx[0]);
      ne_center2 = gsl_vector_get(&qd.vector, ne_idx[1]);

      if (ne_center1 * ne_center2 > 0.0)
        {
          fprintf(stderr, "correlate: detected two peaks in the same hemisphere\n");
          continue;
        }

      /* search for corresponding peaks in F data */
      {
        const double qd_width = 5.0;

        /* search for an F peak within ne_center +/- qd_width */
        npeak_F = find_F_peaks(ne_center1 - qd_width, ne_center1 + qd_width, width_points, &qd.vector, &F.vector, F_idx, peak_F);
        if (npeak_F != 1)
          continue; /* no corresponding F peak found */

        npeak_F = find_F_peaks(ne_center2 - qd_width, ne_center2 + qd_width, width_points, &qd.vector, &F.vector, &(F_idx[1]), peak_F);
        if (npeak_F != 1)
          continue; /* no corresponding F peak found */
      }

      ne_idx0[0] = ne_idx[0] - half_width;
      ne_idx1[0] = ne_idx[0] + half_width;
      ne_idx0[1] = ne_idx[1] - half_width;
      ne_idx1[1] = ne_idx[1] + half_width;

      F_idx0[0] = F_idx[0] - half_width;
      F_idx1[0] = F_idx[0] + half_width;
      F_idx0[1] = F_idx[1] - half_width;
      F_idx1[1] = F_idx[1] + half_width;

      fprintf(stderr, "correlate: track %zu: found %zu N_e peaks\n", i + 1, npeak_ne);

      {
        size_t n1 = ne_idx1[0] - ne_idx0[0] + 1; /* number of points for correlation of peak 1 */
        size_t n2 = ne_idx1[1] - ne_idx0[1] + 1; /* number of points for correlation of peak 2 */

        corr1 = gsl_stats_correlation(&ne.vector.data[ne_idx0[0]], 1, &tptr->Bf[F_idx0[0]], 1, n1);
        corr2 = gsl_stats_correlation(&ne.vector.data[ne_idx0[1]], 1, &tptr->Bf[F_idx0[1]], 1, n2);

        /* make sure corr1 is always the peak in the northern hemisphere */
        if (ne_center1 < 0.0)
          {
            SWAP(corr1, corr2);
          }
      }

      for (j = 0; j < n; ++j)
        {
          double qdj = gsl_vector_get(&qd.vector, j);
          double ne1, ne2, F1, F2;

          if (j >= ne_idx0[0] && j <= ne_idx1[0])
            ne1 = gsl_vector_get(&ne.vector, j);
          else
            ne1 = 1.0/0.0;

          if (j >= ne_idx0[1] && j <= ne_idx1[1])
            ne2 = gsl_vector_get(&ne.vector, j);
          else
            ne2 = 1.0/0.0;

          if (j >= F_idx0[0] && j <= F_idx1[0])
            F1 = tptr->Bf[j];
          else
            F1 = 1.0/0.0;

          if (j >= F_idx0[1] && j <= F_idx1[1])
            F2 = tptr->Bf[j];
          else
            F2 = 1.0/0.0;

          fprintf(fp_peak, "%f %.12e %f %.12e %.12e %f %f %.12e %.12e %f %f\n",
                  qdj,
                  gsl_vector_get(&ne.vector, j),
                  tptr->Bf[j],
                  peak_deriv(j, peak_ne),
                  peak_sderiv(j, peak_ne),
                  peak_deriv(j, peak_F),
                  peak_sderiv(j, peak_F),
                  ne1,
                  ne2,
                  F1,
                  F2);
        }

      fprintf(fp_peak, "\n\n");

      unix_time = satdata_epoch2timet(tptr->t_eq);
      fprintf(fp_corr, "%ld %f %f %f %f %f %f %e %f %e %f %f %f %f %f %f\n",
              unix_time,
              satdata_epoch2year(tptr->t_eq),
              tptr->lon_eq,
              tptr->lt_eq,
              tptr->meanalt,
              get_season(unix_time),
              gsl_vector_get(&qd.vector, ne_idx[0]),
              gsl_vector_get(&ne.vector, ne_idx[0]),
              gsl_vector_get(&qd.vector, ne_idx[1]),
              gsl_vector_get(&ne.vector, ne_idx[1]),
              gsl_vector_get(&qd.vector, F_idx[0]),
              gsl_vector_get(&F.vector, F_idx[0]),
              gsl_vector_get(&qd.vector, F_idx[1]),
              gsl_vector_get(&F.vector, F_idx[1]),
              corr1,
              corr2);

      /*exit(1);*/

      peak_free(peak_ne);
      peak_free(peak_F);
    }

  for (i = 0; i < MAX_PEAKS; ++i)
    gsl_vector_free(ne_coef[i]);

  fclose(fp_peak);
  fclose(fp_corr);

  return s;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --swarm_file | -s swarm_index_file          - Swarm index file\n");
  fprintf(stderr, "\t --champ_file | -c champ_index_file          - CHAMP index file\n");
  fprintf(stderr, "\t --champ_plp_file | -b champ_plp_index_file  - CHAMP PLP index file\n");
  fprintf(stderr, "\t --swarm_lp_file | -f swarm_lp_index_file    - Swarm LP index file\n");
  fprintf(stderr, "\t --all | -a                                  - print all tracks (no filtering)\n");
  fprintf(stderr, "\t --lt_min | -j lt_min                        - local time minimum\n");
  fprintf(stderr, "\t --lt_max | -k lt_max                        - local time maximum\n");
  fprintf(stderr, "\t --alt_min | -l alt_min                      - altitude minimum\n");
  fprintf(stderr, "\t --alt_max | -m alt_max                      - altitude maximum\n");
}

int
main(int argc, char *argv[])
{
  const char *peak_file = "peak.dat";
  const char *corr_file = "corr.dat";
  const char *track_file = "track.dat";
  satdata_mag *data = NULL;
  satdata_lp *lp_data = NULL;
  struct timeval tv0, tv1;
  track_workspace *track_p;
  preprocess_parameters params;

  /* defaults */
  params.all = 0;
  params.lt_min = -1.0;
  params.lt_max = -1.0;
  params.alt_min = 0.0;
  params.alt_max = 1000.0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "all", no_argument, NULL, 'a' },
          { "champ_file", required_argument, NULL, 'c' },
          { "swarm_file", required_argument, NULL, 's' },
          { "champ_plp_file", required_argument, NULL, 'b' },
          { "swarm_plp_file", required_argument, NULL, 'f' },
          { "lt_min", required_argument, NULL, 'j' },
          { "lt_max", required_argument, NULL, 'k' },
          { "alt_min", required_argument, NULL, 'l' },
          { "alt_max", required_argument, NULL, 'm' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "ab:c:f:j:k:l:m:s:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            params.all = 1;
            break;

          case 's':
            data = read_swarm(optarg);
            break;

          case 'c':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_champ_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu data read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));

            /* check for instrument flags since we use Stage1 data */
            {
              size_t nflag;
              size_t champ_flags = 0;

              fprintf(stderr, "main: filtering for instrument flags...");
              nflag = satdata_champ_filter_instrument(1, champ_flags, data);
              fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
                      nflag, data->n, (double)nflag / (double)data->n * 100.0);
            }

            break;

          case 'b':
            fprintf(stderr, "main: reading CHAMP PLP data from %s...", optarg);
            lp_data = satdata_champ_plp_read_idx(optarg, 0);
            fprintf(stderr, "done\n");
            break;

          case 'f':
            fprintf(stderr, "main: reading Swarm LP data from %s...", optarg);
            lp_data = satdata_swarm_lp_read_idx(optarg, 0);
            fprintf(stderr, "done\n");
            break;

          case 'j':
            params.lt_min = atof(optarg);
            break;

          case 'k':
            params.lt_max = atof(optarg);
            break;

          case 'l':
            params.alt_min = atof(optarg);
            break;

          case 'm':
            params.alt_max = atof(optarg);
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

  fprintf(stderr, "main: LT minimum       = %.1f\n", params.lt_min);
  fprintf(stderr, "main: LT maximum       = %.1f\n", params.lt_max);
  fprintf(stderr, "main: altitude minimum = %.1f\n", params.alt_min);
  fprintf(stderr, "main: altitude maximum = %.1f\n", params.alt_max);

  if (lp_data)
    {
      fprintf(stderr, "main: adding electron densities to magnetic data...");
      satdata_mag_fill_ne(data, lp_data);
      fprintf(stderr, "done\n");
    }

  track_p = preprocess_data(&params, data);

  fprintf(stderr, "main: writing track data to %s...", track_file);
  track_print_stats_flag(track_file, 0, track_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: computing correlations...");
  correlate(peak_file, corr_file, data, track_p);
  fprintf(stderr, "done (outputs: %s, %s)\n", peak_file, corr_file);

  satdata_mag_free(data);
  track_free(track_p);

  if (lp_data)
    satdata_lp_free(lp_data);

  return 0;
} /* main() */
