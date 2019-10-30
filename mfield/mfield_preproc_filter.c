/*
 * mfield_preproc_filter.c
 *
 * This module contains routines for filtering data according
 * to various indices requirements
 */

static size_t mfield_flag_kp(const double kp_min, const double kp_max, satdata_mag *data, track_workspace *w);
static size_t mfield_flag_RC(const double RC_max, track_workspace *track_p, satdata_mag *data);
static size_t mfield_flag_dRC(const double dRC_max, track_workspace *track_p, satdata_mag *data);
static size_t mfield_flag_LT(const size_t magdata_flags, const preprocess_parameters * params,
                             track_workspace *track_p, satdata_mag *data);
static size_t mfield_flag_EEJ(const size_t magdata_flags, const preprocess_parameters * params,
                              track_workspace *track_p, satdata_mag *data);
static size_t mfield_flag_zenith(const preprocess_parameters * params, track_workspace *track_p, satdata_mag *data);
static int mfield_preprocess_filter(const size_t magdata_flags, const preprocess_parameters *params,
                                    track_workspace *track_p, satdata_mag *data);
static int mfield_calc_residual(const size_t idx, double B[4], satdata_mag * data);

/*
mfield_flag_kp()
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
mfield_flag_kp(const double kp_min, const double kp_max, satdata_mag *data, track_workspace *w)
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
          fprintf(stderr, "mfield_flag_kp: error: kp not available for track %zu\n", i);
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
mfield_flag_RC()
  Flag tracks for RC criteria

Inputs: RC_max  - maximum allowed RC index (nT)
        track_p - track workspace
        data    - data

Return: number of tracks flagged
*/

static size_t
mfield_flag_RC(const double RC_max, track_workspace *track_p, satdata_mag *data)
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
      double RC1, RC2, RC3;

      s = rc_get_RC(t1, &RC1, rc_p);
      s += rc_get_RC(t2, &RC2, rc_p);
      s += rc_get_RC(t3, &RC3, rc_p);
      if (s)
        {
          fprintf(stderr, "preprocess_flag_RC: error: RC not available for track %zu\n", i);
          continue;
        }

      if ((fabs(RC1) > RC_max) || (fabs(RC2) > RC_max) || (fabs(RC3) > RC_max))
        {
          track_flag_track(i, TRACK_FLG_RC, data, track_p);
          ++ntrack_flagged;
        }
    }

  rc_free(rc_p);

  return ntrack_flagged;
}

/*
mfield_flag_dRC()
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
mfield_flag_dRC(const double dRC_max, track_workspace *track_p, satdata_mag *data)
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
mfield_flag_LT()
  Flag points for local time criteria. Points are analyzed individually
since different criteria are used for mid and high-latitudes

Reject point if:

mid-latitudes:
  ! LT_eq \in [min_LT,max_LT] and ! LT_eq \in [align_min_LT,align_max_LT]

high-latitudes:
  zenith < min_zenith

Inputs: magdata_flags - MAGDATA_GLOBFLG_xxx
        params        - parameters
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
mfield_flag_LT(const size_t magdata_flags, const preprocess_parameters * params,
               track_workspace *track_p, satdata_mag *data)
{
  size_t nflagged = 0; /* number of points flagged */
  size_t i, j;

  if (data->n == 0)
    return 0;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      double LT = tptr->lt_eq;
      int good_MF = mfield_check_LT(LT, params->min_LT, params->max_LT);
      int good_align = mfield_check_LT(LT, params->align_min_LT, params->align_max_LT) &&
                       (magdata_flags & MAGDATA_GLOBFLG_ALIGN);
      int flag_LT = !good_MF && !good_align;

      if (!flag_LT)
        continue;

      for (j = start_idx; j <= end_idx; ++j)
        {
          if (fabs(data->qdlat[j]) <= params->qdlat_preproc_cutoff)
            {
              /* we are at mid-latitudes and this point won't be used for field modeling
               * or Euler angle fitting, so flag the point */
              data->flags[j] |= SATDATA_FLG_OUTLIER;
              ++nflagged;
            }
        }
    }

  return nflagged;
}

/*
mfield_flag_EEJ()
  Flag points in the EEJ region. This is used for Cryosat when we
use both ascending and descending tracks.

Reject point if:

LT_eq \in [6,17] and |qdlat| < 20 deg

Inputs: magdata_flags - MAGDATA_GLOBFLG_xxx
        params        - parameters
        track_p       - track workspace
        data          - data

Return: number of points flagged
*/

static size_t
mfield_flag_EEJ(const size_t magdata_flags, const preprocess_parameters * params,
                track_workspace *track_p, satdata_mag *data)
{
  size_t nflagged = 0; /* number of points flagged */
  size_t i, j;

  if (data->n == 0)
    return 0;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      double LT = tptr->lt_eq;

      if (LT < 6.0 || LT > 17.0)
        continue; /* don't flag night-side data */

      for (j = start_idx; j <= end_idx; ++j)
        {
          if (fabs(data->qdlat[j]) <= 20.0)
            {
              /* we are at low-latitudes on the dayside, so flag point */
              data->flags[j] |= SATDATA_FLG_OUTLIER;
              ++nflagged;
            }
        }
    }

  return nflagged;
}

/*
mfield_flag_zenith()
  Flag high-latitude points for zenith angle criteria.

Reject point if:

1. qdlat > qdlat_preproc_cutoff
2. zenith < min_zenith

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
mfield_flag_zenith(const preprocess_parameters * params, track_workspace *track_p, satdata_mag *data)
{
  size_t nflagged = 0; /* number of points flagged */
  size_t i, j;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;

      for (j = start_idx; j <= end_idx; ++j)
        {
          int flag_point = 0;

          if (fabs(data->qdlat[j]) > params->qdlat_preproc_cutoff)
            {
              time_t unix_time = satdata_epoch2timet(data->t[j]);
              double lat_rad = data->latitude[j] * M_PI / 180.0;
              double phi = data->longitude[j] * M_PI / 180.0;
              double zenith;

              solarpos_calc_zenith(unix_time, lat_rad, phi, &zenith, solarpos_workspace_p);
              zenith *= 180.0 / M_PI;
              assert(zenith >= 0.0);

              /* small zenith angle means sunlit */
              if (zenith < params->min_zenith)
                flag_point = 1;
            }

          if (flag_point)
            {
              data->flags[j] |= SATDATA_FLG_OUTLIER;
              ++nflagged;
            }
        }
    }

  return nflagged;
}

/*
mfield_flag_IMF()
  Flag high-latitude points for IMF criteria

Reject point if:

1. qdlat > qdlat_preproc_cutoff
2. IMF_Bz < 0

Inputs: params        - parameters
        track_p       - track workspace
        data          - data

Return: number of points flagged
*/

static size_t
mfield_flag_IMF(const preprocess_parameters * params, track_workspace *track_p, satdata_mag *data)
{
  ace_workspace *ace_p = ace_alloc(ACE_IDX_FILE);
  size_t nflagged = 0; /* number of points flagged */
  size_t i, j;
  int s;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      time_t t1 = satdata_epoch2timet(data->t[start_idx]);
      time_t t2 = satdata_epoch2timet(tptr->t_eq);
      time_t t3 = satdata_epoch2timet(data->t[end_idx]);
      double IMF_B1[3], IMF_B2[3], IMF_B3[3], SW_vel;
      int flag_IMF = 0;

      s = ace_get(t1, IMF_B1, &SW_vel, ace_p);
      s += ace_get(t2, IMF_B2, &SW_vel, ace_p);
      s += ace_get(t3, IMF_B3, &SW_vel, ace_p);
      if (s)
        {
          /* missing data - flag track since we don't know what IMF is */
          fprintf(stderr, "mfield_flag_IMF: IMF data not available for track %zu [t = %ld]\n", i, t2);
          flag_IMF = 1;
        }
      else
        {
          /* flag if IMF_Bz is outside of [0,6] nT, or if IMF By is > 8 nT */
          if ((IMF_B1[2] < 0.0) || (IMF_B1[2] > 6.0) ||
              (IMF_B2[2] < 0.0) || (IMF_B2[2] > 6.0) ||
              (IMF_B3[2] < 0.0) || (IMF_B3[2] > 6.0) ||
              (IMF_B1[1] > 8.0) || (IMF_B2[1] > 8.0) || (IMF_B3[1] > 8.0))
            {
              flag_IMF = 1;
            }
        }

      if (!flag_IMF)
        continue;

      /* flag all high-latitude points in this track */

      for (j = start_idx; j <= end_idx; ++j)
        {
          if (fabs(data->qdlat[j]) > params->qdlat_preproc_cutoff)
            {
              data->flags[j] |= SATDATA_FLG_OUTLIER;
              ++nflagged;
            }
        }
    }

  ace_free(ace_p);

  return nflagged;
}

/*
mfield_flag_pb()
  Attempt to detect plasma bubbles and flag any possible tracks.
We compute residuals by removing a core, crust and external field.
Then we calcualte N/S gradients between the residuals at the current
point and the along-track point. Then we search if any low-latitude
N/S gradients exceed a given threshold.

Reject track if:

|dB_ns| > eps and qdlat \in [-qdmax,qdmax] (low-latitudes)

where:

1. dB_ns_{ij} = dB_i - dB_j (i is measurement point and j is along-track point)
2. dB = B - B_core - B_crust - B_ext
3. eps is a user-defined threshold

Inputs: params  - parameters
                  pb_qdmax               - search for PBs with qdlat \in [-qdmax,qdmax]
                  pb_threshold_{X,Y,Z,F} - thresholds for N/S gradients (nT)
        track_p - track workspace
        data    - data

Return: number of tracks flagged

Notes:
1) This function should be called after LT selection to avoid
comparing gradient data during daytime
*/

static size_t
mfield_flag_pb(const char *pb_datafile, const preprocess_parameters * params, track_workspace *track_p, satdata_mag *data)
{
  const size_t grad_idx = (size_t) params->gradient_ns;
  const double grad_dt_min = params->gradient_ns - 2.0;
  const double grad_dt_max = params->gradient_ns + 2.0;
  FILE *fp;
  size_t ntrack_flagged = 0; /* number of tracks flagged */
  size_t i, j;

  if (data->n == 0)
    return 0;

  fp = fopen(pb_datafile, "w");
  track_print_track(1, fp, NULL, NULL);

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      int flag_track = 0;

      if (tptr->flags)
        continue;

      for (j = start_idx; j <= end_idx; ++j)
        {
          size_t k = GSL_MIN(j + grad_idx, data->n - 1);
          double dt = fabs(data->t[k] - data->t[j]) / 1000.0; /* in s */
          double B[4], B_ns[4];

          if (!SATDATA_AvailableData(data->flags[j]) || !SATDATA_AvailableData(data->flags[k]))
            continue;
          else if (fabs(data->qdlat[j]) > params->pb_qdmax)
            continue;
          else if (dt < grad_dt_min || dt > grad_dt_max)
            continue;

          mfield_calc_residual(j, B, data);
          mfield_calc_residual(k, B_ns, data);

          if (fabs(B[0] - B_ns[0]) > params->pb_thresh[0] ||
              fabs(B[1] - B_ns[1]) > params->pb_thresh[1] ||
              fabs(B[2] - B_ns[2]) > params->pb_thresh[2])
            {
              flag_track = 1;
              break;
            }
        }

      if (flag_track)
        {
          track_flag_track(i, TRACK_FLG_PB, data, track_p);
          track_print_track(0, fp, tptr, data);
          ++ntrack_flagged;
        }
    }

  return ntrack_flagged;
}

/*
mfield_preprocess_filter()
  Select data for geomagnetically quiet periods
*/

static int
mfield_preprocess_filter(const size_t magdata_flags, const preprocess_parameters *params,
                         track_workspace *track_p, satdata_mag *data)
{
  int s = 0;
  size_t nflagged_kp,
         nflagged_RC,
         nflagged_dRC,
         nflagged_pb,
         nflagged_LT,
         nflagged_EEJ,
         nflagged_zenith,
         nflagged_IMF;

  fprintf(stderr, "\n");

  nflagged_kp = mfield_flag_kp(0.0, params->max_kp, data, track_p);
  fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) tracks due to kp [%.1f]\n",
          nflagged_kp, track_p->n, (double) nflagged_kp / (double) track_p->n * 100.0, params->max_kp);

  /*nflagged_RC = mfield_flag_RC(RC_max, track_p, data);*/
  /*fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) tracks due to RC\n",
          nflagged_RC, track_p->n, (double) nflagged_RC / (double) track_p->n * 100.0);*/

  nflagged_dRC = mfield_flag_dRC(params->max_dRC, track_p, data);
  fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) tracks due to dRC/dt [%.1f nT/hour]\n",
          nflagged_dRC, track_p->n, (double) nflagged_dRC / (double) track_p->n * 100.0, params->max_dRC);

  nflagged_LT = mfield_flag_LT(magdata_flags, params, track_p, data);
  fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) mid-latitude points due to LT [cutoff: %.1f deg]\n",
          nflagged_LT, data->n, (double) nflagged_LT / (double) data->n * 100.0, params->qdlat_preproc_cutoff);

  nflagged_zenith = mfield_flag_zenith(params, track_p, data);
  fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) high-latitude points due to zenith angle [cutoff: %.1f deg]\n",
          nflagged_zenith, data->n, (double) nflagged_zenith / (double) data->n * 100.0, params->qdlat_preproc_cutoff);

  if (magdata_flags & MAGDATA_GLOBFLG_CRYOSAT)
    {
      /* we use ascending and descending tracks for Cryosat, so exclude the EEJ region */
      nflagged_EEJ = mfield_flag_EEJ(magdata_flags, params, track_p, data);
      fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) points due to EEJ region\n",
              nflagged_EEJ, data->n, (double) nflagged_EEJ / (double) data->n * 100.0);
    }

#if 0
  nflagged_IMF = mfield_flag_IMF(params, track_p, data);
  fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) high-latitude points due to IMF [cutoff: %.1f deg]\n",
          nflagged_IMF, data->n, (double) nflagged_IMF / (double) data->n * 100.0, params->qdlat_preproc_cutoff);
#endif

  /* look for plasma bubbles last, after LT selection */
  if (params->pb_flag)
    {
      const char *pb_datafile = "PB_tracks.dat";

      nflagged_pb = mfield_flag_pb(pb_datafile, params, track_p, data);
      fprintf(stderr, "\t mfield_preprocess_filter: flagged %zu/%zu (%.1f%%) tracks due to plasma bubbles (tracks written to %s)\n",
              nflagged_pb, track_p->n, (double) nflagged_pb / (double) track_p->n * 100.0, pb_datafile);
    }

  return s;
}

static int
mfield_calc_residual(const size_t idx, double B[4], satdata_mag * data)
{
  int s = 0;
  double B_model[3] = { 0.0, 0.0, 0.0 };

  B_model[0] = SATDATA_VEC_X(data->B_main, idx) +
               SATDATA_VEC_X(data->B_crust, idx) +
               SATDATA_VEC_X(data->B_ext, idx);

  B_model[1] = SATDATA_VEC_Y(data->B_main, idx) +
               SATDATA_VEC_Y(data->B_crust, idx) +
               SATDATA_VEC_Y(data->B_ext, idx);

  B_model[2] = SATDATA_VEC_Z(data->B_main, idx) +
               SATDATA_VEC_Z(data->B_crust, idx) +
               SATDATA_VEC_Z(data->B_ext, idx);

  B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

  B[0] = SATDATA_VEC_X(data->B, idx) - B_model[0];
  B[1] = SATDATA_VEC_Y(data->B, idx) - B_model[1];
  B[2] = SATDATA_VEC_Z(data->B, idx) - B_model[2];
  B[3] = data->F[idx] - B_model[3];

  return s;
}
