/*
 * track.h
 */

#ifndef INCLUDED_track_h
#define INCLUDED_track_h

#include <satdata/satdata.h>

#include <msynth/msynth.h>

/* maximum number of tracks to store */
#define TRACK_MAX            150000

#define TRACK_DEBUG          0

/* along-track rms is computed from [-QDLAT_MAX,QDLAT_MAX] */
#define TRACK_QDLAT_MAX      (55.0)

/* minimum number of along-track data points for rms calculation */
#define TRACK_RMS_NDATA      200

#define TRACK_FLG_RMS        (1 << 0)  /* flagged due to high rms */
#define TRACK_FLG_JUMP       (1 << 1)  /* flagged due to data jump */
#define TRACK_FLG_SATDIR     (1 << 2)  /* flagged due to satellite direction */
#define TRACK_FLG_UT         (1 << 3)  /* flagged due to UT */
#define TRACK_FLG_LT         (1 << 4)  /* flagged due to local time */
#define TRACK_FLG_KP         (1 << 5)  /* flagged due to kp */
#define TRACK_FLG_RC         (1 << 6)  /* flagged due to RC */
#define TRACK_FLG_DOY        (1 << 7)  /* flagged due to doy / season */
#define TRACK_FLG_ALTITUDE   (1 << 8)  /* flagged due to altitude */
#define TRACK_FLG_LONGITUDE  (1 << 9)  /* flagged due to longitude */
#define TRACK_FLG_INCOMPLETE (1 << 10) /* flagged due to missing data */
#define TRACK_FLG_PB         (1 << 11) /* flagged due to plasma bubble detected */
#define TRACK_FLG_TIME       (1 << 12) /* flagged due to timestamp */
#define TRACK_FLG_IMF        (1 << 13) /* flagged due to IMF conditions */

typedef struct
{
  size_t start_idx; /* starting index of track in 'data' */
  size_t end_idx;   /* ending index of track in 'data' */
  double *Bx;       /* residual X component (nT) */
  double *By;       /* residual Y component (nT) */
  double *Bz;       /* residual Z component (nT) */
  double *Bf;       /* residual F component (nT) */
  size_t n;         /* number of data points in this track */
  double t_eq;      /* timestamp of equator crossing (CDF_EPOCH) */
  double lon_eq;    /* longitude of equator crossing (deg) */
  double lat_eq;    /* geocentric latitude of equator crossing (deg) */
  double lt_eq;     /* local time of equator crossing (hours) */
  int satdir;       /* satellite direction: +1 north, -1 south */
  double rms[4];    /* along-track rms (nT) */
  size_t nrms_scal; /* number of data used for scalar rms calculation */
  size_t nrms_vec;  /* number of data used for vector rms calculation */
  double rmspol[4]; /* along-track rms at polar latitudes (nT) */
  size_t nrmspol_scal; /* number of data used for scalar polar rms calculation */
  size_t nrmspol_vec;  /* number of data used for vector polar rms calculation */
  double meanalt;   /* mean track altitude (km) */
  size_t flags;     /* TRACK_FLG_xxx */
  double k_ext;     /* external field correction coefficient (nT) */
} track_data;

typedef struct
{
  track_data *tracks; /* track data */
  size_t n;           /* number of tracks stored */
  size_t ntot;        /* number of tracks allocated */
  satdata_mag *data;  /* satellite data */
  msynth_workspace *msynth_workspace_p;
} track_workspace;

/*
 * Prototypes
 */

track_workspace *track_alloc(void);
void track_free(track_workspace *w);
size_t track_init(satdata_mag *data, msynth_workspace *msynth_p, track_workspace *w);
int track_smooth(const double alpha, satdata_mag *data, track_workspace *w);
int track_residual(const size_t track_idx, const size_t data_idx, double B[3],
                   const satdata_mag *data, track_workspace *w);

size_t track_flag_data(size_t sidx, size_t eidx, satdata_mag *data);
size_t track_flag_track(const size_t track_idx, const size_t flags,
                        satdata_mag *data, track_workspace *w);
size_t track_nflagged(const track_workspace *w);
size_t track_data_nflagged(const track_data *tptr, const satdata_mag *data);
int track_find(const double t_eq, const double phi_eq, const double dt_min,
               const double dphi, size_t *idx, const track_workspace *w);
int track_find_t(const double t, size_t *idx, const track_workspace *w);
int track_print(const char *filename, const size_t flags,
                const satdata_mag *data, track_workspace *w);
int track_print_track(const int header, FILE *fp, const track_data *tptr,
                      const satdata_mag *data);
int track_print_stats_flag(const char *filename, const size_t flag,
                           track_workspace *w);
int track_print_stats(const char *filename, track_workspace *w);
int track_calc_residuals(track_data *tptr, const satdata_mag *data);

/* track_filter.c */
int track_filter(const char *filename, track_workspace *w);

/* track_flag.c */
size_t track_flag_rms(const char *outfile, const double thresh[4],
                      size_t *ndata_flagged, satdata_mag *data, track_workspace *w);
size_t track_flag_jumps(const double thresh, satdata_mag *data, track_workspace *w);
size_t track_flag_satdir(const int satdir, satdata_mag *data, track_workspace *w);
size_t track_flag_time(const double t_min, const double t_max, satdata_mag *data, track_workspace *w);
size_t track_flag_ut(const double ut_min, const double ut_max, satdata_mag *data,
                     track_workspace *w);
size_t track_flag_lt(const double lt_min, const double lt_max, size_t *ndata_flagged,
                     satdata_mag *data, track_workspace *w);
size_t track_flag_kp(const double kp_min, const double kp_max, satdata_mag *data,
                     track_workspace *w);
size_t track_flag_IMF(const double Bz_min, const double Bz_max, satdata_mag *data, track_workspace *w);
size_t track_flag_season(int (*callback)(const double doy, const void *params),
                         const void *params, satdata_mag *data, track_workspace *w);
size_t track_flag_lon(const double lon_min, const double lon_max,
                      size_t *ndata_flagged, satdata_mag *data, track_workspace *w);
size_t track_flag_meanalt(const double alt_min, const double alt_max, satdata_mag *data,
                          track_workspace *w);
size_t track_flag_incomplete(const double qd_min, const double qd_max, satdata_mag *data,
                             track_workspace *w);
size_t track_flag_n(const size_t nmin, satdata_mag *data, track_workspace *w);

/* track_offsets.c */
int track_fix_offsets(const satdata_mag *data, track_workspace *w);

/* track_synth.c */
int track_synth_int(satdata_mag *data, msynth_workspace *msynth_core_p, msynth_workspace *msynth_crust_p);
int track_synth_core(satdata_mag *data, msynth_workspace *msynth_p);
int track_synth_QD(satdata_mag *data);
int track_synth_MLT(satdata_mag *data);
int track_synth_tilt(satdata_mag *data);
int track_synth_pomme(satdata_mag *data);

#endif /* INCLUDED_track_h */
