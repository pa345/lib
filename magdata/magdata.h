/*
 * magdata.h
 */

#ifndef INCLUDED_magdata_h
#define INCLUDED_magdata_h

#include <satdata/satdata.h>
#include <obsdata/obsdata.h>

#include "track.h"
#include "track_weight.h"

#define MAGDATA_NAME_LENGTH           128

/* data flags */
#define MAGDATA_FLG_X                 (1 << 0)  /* X measurement available (both VFM and NEC) */
#define MAGDATA_FLG_Y                 (1 << 1)  /* Y measurement available (both VFM and NEC) */
#define MAGDATA_FLG_Z                 (1 << 2)  /* Z measurement available (both VFM and NEC) */
#define MAGDATA_FLG_F                 (1 << 3)  /* F measurement available */
#define MAGDATA_FLG_DXDT              (1 << 4)  /* dX/dt measurement available (NEC) */
#define MAGDATA_FLG_DYDT              (1 << 5)  /* dY/dt measurement available (NEC) */
#define MAGDATA_FLG_DZDT              (1 << 6)  /* dZ/dt measurement available (NEC) */
#define MAGDATA_FLG_DX_NS             (1 << 7)  /* along-track DX measurement available (both VFM and NEC) */
#define MAGDATA_FLG_DY_NS             (1 << 8)  /* along-track DY measurement available (both VFM and NEC) */
#define MAGDATA_FLG_DZ_NS             (1 << 9)  /* along-track DZ measurement available (both VFM and NEC) */
#define MAGDATA_FLG_DF_NS             (1 << 10) /* along-track DF measurement available */
#define MAGDATA_FLG_DX_EW             (1 << 11) /* east-west DX measurement available (both VFM and NEC) */
#define MAGDATA_FLG_DY_EW             (1 << 12) /* east-west DY measurement available (both VFM and NEC)*/
#define MAGDATA_FLG_DZ_EW             (1 << 13) /* east-west DZ measurement available (both VFM and NEC) */
#define MAGDATA_FLG_DF_EW             (1 << 14) /* east-west DF measurement available */
#define MAGDATA_FLG_TRACK_START       (1 << 15) /* start of track in data structure */
#define MAGDATA_FLG_FIT_MF            (1 << 16) /* fit internal main field to this data point */
#define MAGDATA_FLG_FIT_EULER         (1 << 17) /* fit Euler angles to this data point */
#define MAGDATA_FLG_DISCARD           (1 << 18) /* discard this data point */

/* global flags */
#define MAGDATA_GLOBFLG_EULER         (1 << 0)  /* fit Euler angles to this dataset */
#define MAGDATA_GLOBFLG_SCALAR_GRID   (1 << 1)  /* dataset is a scalar-only grid like EMAG2 */
#define MAGDATA_GLOBFLG_OBSERVATORY   (1 << 2)  /* dataset is from ground observatory */
#define MAGDATA_GLOBFLG_SATELLITE     (1 << 3)  /* dataset is from satellite */

/* k_b * mu_0 in units of: nT^2 cm^3 / K */
#define MAGDATA_KB_MU0                (1.73497445090703e-05)

/* check if data point is discarded */
#define MAGDATA_Discarded(x)          ( (x) & MAGDATA_FLG_DISCARD )

/* vector X/Y/Z measurements available */
#define MAGDATA_ExistX(x)             (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_X))
#define MAGDATA_ExistY(x)             (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_Y))
#define MAGDATA_ExistZ(x)             (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_Z))

/* total vector measurement available */
#define MAGDATA_ExistVector(x)        (!MAGDATA_Discarded(x) && \
                                       ((x) & MAGDATA_FLG_X) && \
                                       ((x) & MAGDATA_FLG_Y) && \
                                       ((x) & MAGDATA_FLG_Z))

/* scalar measurement available */
#define MAGDATA_ExistScalar(x)        (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_F))

/* vector dX/dt,dY/dt,dZ/dt measurements available */
#define MAGDATA_ExistDXDT(x)          (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DXDT))
#define MAGDATA_ExistDYDT(x)          (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DYDT))
#define MAGDATA_ExistDZDT(x)          (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DZDT))

#define MAGDATA_ExistVectorNS(x)      (!MAGDATA_Discarded(x) && \
                                       ((x) & MAGDATA_FLG_DX_NS) && \
                                       ((x) & MAGDATA_FLG_DY_NS) && \
                                       ((x) & MAGDATA_FLG_DZ_NS))

/* vector north-south gradient X/Y/Z measurements available */
#define MAGDATA_ExistDX_NS(x)         (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DX_NS))
#define MAGDATA_ExistDY_NS(x)         (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DY_NS))
#define MAGDATA_ExistDZ_NS(x)         (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DZ_NS))

/* scalar north-south gradient measurement available */
#define MAGDATA_ExistDF_NS(x)         (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DF_NS))

/* vector east-west gradient X/Y/Z measurements available */
#define MAGDATA_ExistDX_EW(x)         (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DX_EW))
#define MAGDATA_ExistDY_EW(x)         (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DY_EW))
#define MAGDATA_ExistDZ_EW(x)         (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DZ_EW))

/* scalar east-west gradient measurement available */
#define MAGDATA_ExistDF_EW(x)         (!MAGDATA_Discarded(x) && ((x) & MAGDATA_FLG_DF_EW))

/* check if fitting Euler angles to this data point */
#define MAGDATA_FitEuler(x)           ((MAGDATA_ExistVector(x) || MAGDATA_ExistVectorNS(x)) && \
                                       ((x) & MAGDATA_FLG_FIT_EULER))

/* check if fitting MF model to this data point */
#define MAGDATA_FitMF(x)              ((x) & MAGDATA_FLG_FIT_MF)

typedef struct
{
  char name[MAGDATA_NAME_LENGTH]; /* name of dataset */
  double t;              /* timestamp (CDF_EPOCH) */
  double r;              /* radius (km) */
  double theta;          /* geocentric colatitude (radians) */
  double phi;            /* longitude (radians) in [-pi,pi] */
  double qdlat;          /* QD latitude (degrees) */
  double MLT;            /* magnetic local time (hours) */
  double B_nec[3];       /* magnetic field vector in NEC frame (nT) */
  double B_vfm[3];       /* magnetic field vector in VFM frame (nT) */
  double B_model[3];     /* main/crustal/external field model in NEC frame (nT) */
  double F;              /* magnetic field scalar measurement (nT) */
  double dBdt_nec[3];    /* time derivative of magnetic field vector in NEC frame (nT/year) */
  double q[4];           /* quaternions */
  double ne;             /* electron density cm^{-3} */
  int satdir;            /* +1 north, -1 south */
  double lt;             /* local time in [0,24] */
  double lt_eq;          /* local time of equator crossing */

  double t_ns;           /* timestamp for along-track point (CDF_EPOCH) */
  double r_ns;           /* radius for along-track point (km) */
  double theta_ns;       /* geocentric (or QD) colatitude for along-track point (radians) */
  double phi_ns;         /* longitude for along-track point (radians) */
  double qdlat_ns;       /* QD latitude for along-track point (degrees) */
  double B_nec_ns[3];    /* along-track magnetic field vector in NEC frame (nT) */
  double B_vfm_ns[3];    /* along-track magnetic field vector in VFM frame (nT) */
  double B_model_ns[3];  /* along-track main/crustal/external field model in NEC frame (nT) */
  double F_ns;           /* along-track magnetic field scalar measurement (nT) */
  double q_ns[4];        /* along-track quaternions */
  double lt_ns;          /* local time of gradient point in [0,24] */
  double lt_eq_ns;       /* local time of equator crossing for gradient point */

  size_t flags;          /* MAGDATA_FLG_xxx flags */
} magdata_datum;

typedef struct
{
  char name[MAGDATA_NAME_LENGTH]; /* name of dataset */
  double *t;           /* timestamp (CDF_EPOCH) */
  double *ts;          /* scaled timestamp for MF modeling (dimensionless), not written to output file */
  double *r;           /* radius (km) */
  double *theta;       /* geocentric colatitude (radians) */
  double *phi;         /* longitude (radians) */
  double *qdlat;       /* QD latitude (degrees) */
  double *MLT;         /* magnetic local time (hours) */
  double *Bx_nec;      /* NEC X measurement (nT) */
  double *By_nec;      /* NEC Y measurement (nT) */
  double *Bz_nec;      /* NEC Z measurement (nT) */
  double *Bx_vfm;      /* VFM X measurement (nT) */
  double *By_vfm;      /* VFM Y measurement (nT) */
  double *Bz_vfm;      /* VFM Z measurement (nT) */
  double *Bx_model;    /* NEC X main/crustal/external field (nT) */
  double *By_model;    /* NEC Y main/crustal/external field (nT) */
  double *Bz_model;    /* NEC Z main/crustal/external field (nT) */
  double *F;           /* scalar measurement (nT) */
  double *dXdt_nec;    /* NEC dX/dt measurement (nT/year) */
  double *dYdt_nec;    /* NEC dY/dt measurement (nT/year) */
  double *dZdt_nec;    /* NEC dZ/dt measurement (nT/year) */
  double *q;           /* quaternions */
  double *weights;     /* spatial weights */
  int *satdir;         /* +1 north, -1 south */
  double *lt;          /* local time in [0,24] */
  double *lt_eq;       /* local time of equator crossing */

  double *ne;          /* electron density cm^{-3} */

  double *t_ns;        /* timestamp for along-track point (CDF_EPOCH) */
  double *ts_ns;       /* along-track scaled timestamp for MF modeling (dimensionless), not written to output file */
  double *r_ns;        /* radius for along-track point (km) */
  double *theta_ns;    /* geocentric colatitude for along-track point (radians) */
  double *phi_ns;      /* longitude for along-track point (radians) */
  double *qdlat_ns;    /* QD latitude for along-track point (degrees) */
  double *Bx_nec_ns;   /* NEC along-track X measurement (nT) */
  double *By_nec_ns;   /* NEC along-track Y measurement (nT) */
  double *Bz_nec_ns;   /* NEC along-track Z measurement (nT) */
  double *Bx_vfm_ns;   /* VFM along-track X measurement (nT) */
  double *By_vfm_ns;   /* VFM along-track Y measurement (nT) */
  double *Bz_vfm_ns;   /* VFM along-track Z measurement (nT) */
  double *Bx_model_ns; /* NEC along-track X main/crustal/external field (nT) */
  double *By_model_ns; /* NEC along-track Y main/crustal/external field (nT) */
  double *Bz_model_ns; /* NEC along-track Z main/crustal/external field (nT) */
  double *F_ns;        /* along-track scalar measurement (nT) */
  double *q_ns;        /* along-track quaternions */
  double *lt_ns;       /* local time of gradient point in [0,24] */
  double *lt_eq_ns;    /* along-track local time of equator crossing */

  size_t *flags;       /* MAGDATA_FLG_xxx flags */

  size_t *index;       /* indexing of residuals, skipping flagged data */

  size_t n;            /* total number of data */
  size_t ntot;         /* total array allocation */
  double R;            /* Earth radius (km) */
  double rmin;         /* minimum radius stored (km) */
  double rmax;         /* maximum radius stored (km) */
  size_t nx;           /* number of X measurements */
  size_t ny;           /* number of Y measurements */
  size_t nz;           /* number of Z measurements */
  size_t nf;           /* number of F measurements */
  size_t ndXdt;        /* number of dX/dt measurements */
  size_t ndYdt;        /* number of dY/dt measurements */
  size_t ndZdt;        /* number of dZ/dt measurements */
  size_t ndx;          /* number of DX measurements */
  size_t ndy;          /* number of DY measurements */
  size_t ndz;          /* number of DZ measurements */
  size_t ndf;          /* number of DF measurements */
  size_t nvec;         /* number of vector measurements */
  size_t ngrad;        /* number of along-track difference measurements */
  size_t nres;         /* number of total residuals */

  size_t euler_flags;  /* EULER_FLG_xxx flags for Euler angle convention */
  size_t global_flags; /* MAGDATA_GLOBFLG_xxx flags applying to all data */

  track_weight_workspace *weight_workspace_p;
} magdata;

/* parameters for copying tracks into magdata structure */
typedef struct
{
  double grad_dt_ns;    /* time interval for along-track differences (seconds) */
  double grad_dt_ew;    /* maximum time interval for E/W differences (seconds) */
  double grad_dphi_max; /* maximum allowed longitudinal separation for E/W gradients (degrees) */
  double grad_dlat_max; /* maximum allowed latitudinal separation for E/W gradients (degrees) */
  int model_main;       /* include main field in B_model */
  int model_crust;      /* include crustal field in B_model */
  int model_ext;        /* include external field in B_model */
} magdata_params;

/* parameters for preprocessing data */
typedef struct
{
  size_t downsample;      /* downsampling factor */
  double min_LT;          /* minimum local time for field modeling */
  double max_LT;          /* maximum local time for field modeling */
  double euler_min_LT;    /* minimum local time for Euler angles */
  double euler_max_LT;    /* maximum local time for Euler angles */
  double rms_thresh[4];   /* rms thresholds (X,Y,Z,F) (nT) */
  double qdlat_preproc_cutoff; /* QD latitude cutoff for high-latitudes */
  double min_zenith;      /* minimum zenith angle for high-latitude data selection */
  size_t gradient_ns;     /* number of seconds between N/S gradient samples */
  int fit_track_RC;       /* fit track-by-track RC field */
  int flag_IMF;           /* flag high-latitude points according to IMF criteria */

  double season_min;      /* season minimum [0,366] */
  double season_max;      /* season maximum [0,366] */
  double season_min2;     /* season minimum [0,366] */
  double season_max2;     /* season maximum [0,366] */

  double rmin;            /* minimum radius (km) */
  double rmax;            /* maximum radius (km) */

  double gradew_dphi_max; /* maximum longitude distance for east-west gradients (degrees) */
  double gradew_dlat_max; /* maximum latitude distance for east-west gradients (degrees) */
  double gradew_dt_max;   /* maximum time difference for east-west gradients (seconds) */

  int subtract_B_main;    /* subtract a-priori main field from data */
  int subtract_B_crust;   /* subtract a-priori crustal field from data */
  int subtract_B_ext;     /* subtract a-priori external field from data */

  double max_kp;          /* maximum kp */
  double max_dRC;         /* maximum dRC/dt (nT/hour) */

  int polar_damping;      /* apply cosine damping window to polar data */
  double polar_qdlat;     /* QD latitude cutoff for polar damping window */

  int pb_flag;            /* flag tracks with plasma bubble signatures */
  double pb_qdmax;        /* QD latitude range for PB search */
  double pb_thresh[4];    /* threshold values for N/S gradients (X,Y,Z,F) (nT) */
} magdata_preprocess_parameters;

/*
 * Prototypes
 */

magdata *magdata_alloc(const size_t n, const double R);
magdata *magdata_realloc(const size_t n, const double R, magdata *data);
void magdata_free(magdata *data);
int magdata_set_euler(const size_t flags, magdata *data);
int magdata_datum_init(magdata_datum *datum);
int magdata_add(const magdata_datum *datum, magdata *data);
int magdata_init(magdata *data);
int magdata_unit_weights(magdata *data);
int magdata_calc(magdata *data);
int magdata_print(const char *prefix, const magdata *data);
int magdata_map(const char *prefix, const magdata *data);
int magdata_residual(const size_t idx, double B[4], const magdata *data);
int magdata_residual_ns(const size_t idx, double B[4], const magdata *data);
int magdata_residual_dB_ns(const size_t idx, double dB[4], const magdata *data);
int magdata_prior(const size_t idx, double B_prior[4], const magdata *data);
int magdata_prior_grad(const size_t idx, double B_prior[4], const magdata *data);
int magdata_t(double *t0, double *t1, const magdata *data);
int magdata_write(const char *filename, magdata *data);
magdata *magdata_read(const char *filename, magdata *data);
size_t magdata_ndiscard(const magdata *data);
size_t magdata_neuler(const magdata *data);
int magdata_clear(magdata *data);
int magdata_flag_t(const double t0, const double t1, magdata *data);
int magdata_flag_scalar(magdata *data);
int magdata_copy_track(const magdata_params *params, const size_t track_idx,
                       const satdata_mag *data, const track_workspace *track_p,
                       magdata *mdata, size_t ntype[6]);
int magdata_copy_track_EW(const magdata_params *params, const size_t track_idx,
                          const satdata_mag *data, const track_workspace *track_p,
                          const satdata_mag *data2, const track_workspace *track_p2,
                          magdata *mdata, size_t ntype[6]);
int magdata_copy_station(const magdata_params *params, const obsdata_station * station,
                         magdata *mdata, size_t ntype[6]);
satdata_mag *magdata_mag2sat(const magdata *mdata);
int magdata_replace_phi_LT(const double lt0, magdata *data);

/* preproc.c */
magdata_preprocess_parameters magdata_preprocess_default_parameters(void);
int magdata_preprocess_parse(const char *filename, magdata_preprocess_parameters *params);
int magdata_preprocess_check(magdata_preprocess_parameters * params);
track_workspace *magdata_preprocess(const magdata_preprocess_parameters *params, const size_t magdata_flags,
                                    satdata_mag *data);
magdata *magdata_preprocess_fill(const size_t magdata_flags, const satdata_mag *data, const track_workspace *track_p,
                                 const size_t magdata_flags2, const satdata_mag *data2, const track_workspace *track_p2,
                                 magdata_preprocess_parameters * preproc_params);

#endif /* INCLUDED_magdata_h */
