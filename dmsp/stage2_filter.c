
static size_t filter_kp(const double kp_min, const double kp_max, satdata_mag *data, track_workspace *w);
static size_t filter_dRC(const double dRC_max, track_workspace *track_p, satdata_mag *data);
static size_t filter_LT(const double min_LT, const double max_LT, track_workspace *track_p, satdata_mag *data);

/*
filter_kp()
  Flag any tracks with kp outside of [kp_min,kp_max].
  
4 kp values are compared:

1. beginning of track
2. equator crossing
3. end of track
4. 2 hours prior to beginning of track

Inputs: kp_min - minimum kp
        kp_max - maximum kp
        data   - satellite data
        w      - track workspace

Return: number of tracks flagged
*/

static size_t
filter_kp(const double kp_min, const double kp_max, satdata_mag *data, track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;        /* number of points flagged */
  size_t ntrack_flagged = 0;  /* number of entire tracks flagged for UT */
  kp_workspace *kp_p;
  int s;

  if (data->n == 0)
    return 0;

  kp_p = kp_alloc(KP_IDX_FILE);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      time_t t1 = satdata_epoch2timet(data->t[start_idx]);
      time_t t2 = satdata_epoch2timet(tptr->t_eq);
      time_t t3 = satdata_epoch2timet(data->t[end_idx]);
      double kp1, kp2, kp3, kp4;

      s = kp_get(t1, &kp1, kp_p);
      s += kp_get(t2, &kp2, kp_p);
      s += kp_get(t3, &kp3, kp_p);
      s += kp_get(t1 - 2*3600, &kp4, kp_p); /* 2 hours before track */
      if (s)
        {
          fprintf(stderr, "filter_kp: error: kp not available for track %zu\n", i);
          continue;
        }

      if ((kp1 < kp_min || kp1 > kp_max) ||
          (kp2 < kp_min || kp2 > kp_max) ||
          (kp3 < kp_min || kp3 > kp_max) ||
          (kp4 < kp_min || kp4 > kp_max))
        {
          nflagged += track_flag_track(i, TRACK_FLG_KP, data, w);
          ++ntrack_flagged;
        }
    }

  kp_free(kp_p);

  return ntrack_flagged;
}

/*
filter_dRC()
  Flag tracks for dRC/dt criteria

Reject track if:

1. |dRC/dt| > dRC_max at any time along the track
2. |dRC/dt| > dRC_max within 3 hours prior to track

Inputs: dRC_max - maximum allowed dRC/dt index (nT/hour)
        track_p - track workspace
        data    - data

Return: number of tracks flagged
*/

static size_t
filter_dRC(const double dRC_max, track_workspace *track_p, satdata_mag *data)
{
  size_t ntrack_flagged = 0; /* number of tracks flagged */
  rc_workspace *rc_p;
  size_t i;
  int s;

  if (data->n == 0)
    return 0;

  rc_p = rc_alloc(RC_IDX_FILE);

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      time_t t1 = satdata_epoch2timet(data->t[start_idx]);
      time_t t2 = satdata_epoch2timet(tptr->t_eq);
      time_t t3 = satdata_epoch2timet(data->t[end_idx]);
      double dRC1, dRC2, dRC3;
      int flag_track = 0;

      s = rc_deriv_get(t1, &dRC1, rc_p);
      s += rc_deriv_get(t2, &dRC2, rc_p);
      s += rc_deriv_get(t3, &dRC3, rc_p);
      if (s)
        {
          fprintf(stderr, "preprocess_flag_RC: error: dRC/dt not available for track %zu\n", i);
          flag_track = 1;
        }

      if ((fabs(dRC1) > dRC_max) || (fabs(dRC2) > dRC_max) || (fabs(dRC3) > dRC_max))
        {
          flag_track = 1;
        }

      /* now check |dRC/dt| for 3 hours prior to track */
      s = rc_deriv_get(t1 - 3600, &dRC1, rc_p);
      s += rc_deriv_get(t1 - 2*3600, &dRC2, rc_p);
      s += rc_deriv_get(t1 - 3*3600, &dRC3, rc_p);
      if (s)
        {
          fprintf(stderr, "preprocess_flag_RC: error: dRC/dt not available for track %zu\n", i);
          flag_track = 1;
        }

      if ((fabs(dRC1) > dRC_max) || (fabs(dRC2) > dRC_max) || (fabs(dRC3) > dRC_max))
        {
          flag_track = 1;
        }

      if (flag_track)
        {
          track_flag_track(i, TRACK_FLG_RC, data, track_p);
          ++ntrack_flagged;
        }
    }

  rc_free(rc_p);

  return ntrack_flagged;
}

/*
filter_LT()
  Flag points for local time criteria. Points are analyzed individually
since different criteria are used for mid and high-latitudes

Reject point if:

mid-latitudes:
  ! LT_eq \in [min_LT,max_LT] and ! LT_eq \in [euler_min_LT,euler_max_LT]

high-latitudes:
  zenith < min_zenith

Inputs: params        - parameters
        track_p       - track workspace
        data          - data

Return: number of points flagged

Notes:
1) The purpose of this routine is to flag points which don't fit
   our LT/zenith criteria so they won't be copied into the magdata
   structure. However, once the magdata structure is built, we need
   to repeat this process in order to flag points for MF and/or Euler
   angle fitting.
*/

static size_t
filter_LT(const double min_LT, const double max_LT,
          track_workspace *track_p, satdata_mag *data)
{
  size_t ntrack_flagged = 0; /* number of points flagged */
  size_t i;

  if (data->n == 0)
    return 0;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      double LT = tptr->lt_eq;
      int good_LT = check_LT(LT, min_LT, max_LT);

      if (!good_LT)
        {
          track_flag_track(i, TRACK_FLG_RC, data, track_p);
          ++ntrack_flagged;
        }
    }

  return ntrack_flagged;
}

/*
stage2_filter()
  Select data for geomagnetically quiet periods
*/

static int
stage2_filter(const size_t downsample, const double thresh[4], track_workspace *track_p, satdata_mag *data)
{
  int s = 0;
  const char *rms_file = "satrms.dat";
  const double max_kp = 3.0;
  const double max_dRC = 5.0; /* nT/hour */
  const double min_LT = 17.0;
  const double max_LT = 5.0;
  size_t nflagged_rms,
         nflagged_kp,
         nflagged_dRC,
         nflagged_LT,
         nflagged_lon,
         nflagged_tot;
  size_t i;

  fprintf(stderr, "\n");

  /*
   * recompute along-track residuals in case previous processing steps have modified data;
   * otherwise rms test could fail
   */
  fprintf(stderr, "\t stage2_filter: recomputing track residuals...");

  for (i = 0; i < track_p->n; ++i)
    track_calc_residuals(&(track_p->tracks[i]), data);

  fprintf(stderr, "done\n");

  nflagged_rms = track_flag_rms(rms_file, thresh, NULL, data, track_p);
  fprintf(stderr, "\t stage2_filter: flagged %zu/%zu (%.1f%%) tracks due to scalar rms [%.1f nT]\n",
          nflagged_rms, track_p->n, (double) nflagged_rms / (double) track_p->n * 100.0, thresh[3]);

#if 0
  nflagged_lon = track_flag_lon(-100.0, 150.0, NULL, data, track_p);
  fprintf(stderr, "\t stage2_filter: flagged %zu/%zu (%.1f%%) tracks due to longitude\n",
          nflagged_lon, track_p->n, (double) nflagged_lon / (double) track_p->n * 100.0);
#endif

#if 0
  nflagged_kp = filter_kp(0.0, max_kp, data, track_p);
  fprintf(stderr, "\t stage2_filter: flagged %zu/%zu (%.1f%%) tracks due to kp [%.1f]\n",
          nflagged_kp, track_p->n, (double) nflagged_kp / (double) track_p->n * 100.0, max_kp);

  nflagged_dRC = filter_dRC(max_dRC, track_p, data);
  fprintf(stderr, "\t stage2_filter: flagged %zu/%zu (%.1f%%) tracks due to dRC/dt [%.1f nT/hour]\n",
          nflagged_dRC, track_p->n, (double) nflagged_dRC / (double) track_p->n * 100.0, max_dRC);

  nflagged_LT = filter_LT(min_LT, max_LT, track_p, data);
  fprintf(stderr, "\t stage2_filter: flagged %zu/%zu (%.1f%%) tracks due to LT [%g,%g]\n",
          nflagged_LT, track_p->n, (double) nflagged_LT / (double) track_p->n * 100.0, min_LT, max_LT);
#endif

  nflagged_tot = track_nflagged(track_p);
  fprintf(stderr, "\t stage2_filter: flagged %zu/%zu (%.1f%%) total tracks\n",
          nflagged_tot, track_p->n, (double) nflagged_tot / (double) track_p->n * 100.0);

  /* downsample data */
  {
    size_t i;

    fprintf(stderr, "stage2_filter: downsampling data by factor %zu...", downsample);

    for (i = 0; i < data->n; ++i)
      {
        if (i % downsample != 0)
          data->flags[i] |= SATDATA_FLG_DOWNSAMPLE;
      }

    fprintf(stderr, "done\n");
  }

  return s;
}
