/*
 * mag.h
 */

#ifndef INCLUDED_mag_h
#define INCLUDED_mag_h

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_apex.h>
#include <mainlib/ml_track.h>
#include <mainlib/ml_pomme.h>
#include <mainlib/ml_mageq.h>
#include <mainlib/ml_green.h>

#include "estist_calc.h"
#include "log.h"
#include "inverteef.h"
#include "magfit.h"
#include "pde.h"

/* maximum data in one track */
#define MAG_MAX_TRACK          10000

typedef struct
{
  int year;                       /* year of measurements for QD transforms */
  char *log_dir;                  /* log directory */
  char *output_file;              /* output file */
  double kp_max;                  /* maximum allowed kp */
  double lt_min;                  /* minimum local time (hours) */
  double lt_max;                  /* maximum local time (hours) */
  double lon_min;                 /* minimum longitude (degrees) */
  double lon_max;                 /* maximum longitude (degrees) */
  double season_min;              /* minimum season (doy) */
  double season_max;              /* maximum season (doy) */
  double season_min2;             /* minimum season 2 (doy) */
  double season_max2;             /* maximum season 2 (doy) */
  double track_qdmax;             /* maximum QD latitude for processing */
  size_t sq_nmax_int;             /* Sq filter internal nmax */
  size_t sq_mmax_int;             /* Sq filter internal mmax */
  size_t sq_nmax_ext;             /* Sq filter external nmax */
  size_t sq_mmax_ext;             /* Sq filter external mmax */
  double sq_qdmin;                /* minimum QD latitude for Sq fitting */
  double sq_qdmax;                /* maximum QD latitude for Sq fitting */
  double curr_altitude;           /* altitude of line currents (km) */
  size_t ncurr;                   /* number of line currents */
  double qdlat_max;               /* maximum QD latitude for line currents (deg) */
  double dlat_max;                /* maximum latitude data gap allowed (deg) */
  double dtime_max;               /* maximum time data gap allowed (minutes) */
  int profiles_only;              /* compute profiles only (no EEF) */
  int use_vector;                 /* use vector data instead of scalar */
  double r_earth;                 /* Earth radius (km) */
  int calc_field_models;          /* compute along-track field models? */
  int main_nmax_int;              /* spherical harmonic nmax for core field model */
  int crust_nmax_int;             /* spherical harmonic nmax for crustal field model */

  char *prev_day_file;            /* previous day MAGx_LR file */
  char *curr_day_file;            /* current day MAGx_LR file */
  char *kp_file;                  /* KP data file */
  char *f107_file;                /* F10.7 data file */
  char *dst_file;                 /* DST data file */
  char *estist_file;              /* Est/Ist data file */

  char *core_file;                /* Core field file */
  char *lith_file;                /* lithospheric field file */
} mag_params;

/* store 1 satellite track */
typedef struct
{
  double t[MAG_MAX_TRACK];        /* timestamp (CDF_EPOCH) */
  double theta[MAG_MAX_TRACK];    /* colatitude in radians */
  double phi[MAG_MAX_TRACK];      /* longitude in radians */
  double r[MAG_MAX_TRACK];        /* radius in km */
  double thetaq[MAG_MAX_TRACK];   /* QD colatitude in radians */
  double qdlat[MAG_MAX_TRACK];    /* QD latitude in degrees */
  double lat_deg[MAG_MAX_TRACK];  /* geocentric latitude in degrees */
  double F[MAG_MAX_TRACK];        /* scalar field measurement in nT */
  double X[MAG_MAX_TRACK];        /* X field measurement in nT */
  double Y[MAG_MAX_TRACK];        /* Y field measurement in nT */
  double Z[MAG_MAX_TRACK];        /* Z field measurement in nT */
  double Bx_int[MAG_MAX_TRACK];   /* Bx_internal NEC (nT) */
  double By_int[MAG_MAX_TRACK];   /* By_internal NEC (nT) */
  double Bz_int[MAG_MAX_TRACK];   /* Bz_internal NEC (nT) */
  double F_int[MAG_MAX_TRACK];    /* F_internal (nT) */
  double dF_ext[MAG_MAX_TRACK];   /* b_int . B_ext (nT) */
  double F1[MAG_MAX_TRACK];       /* F_sat - F_int - dF_ext (nT) */
  double F2[MAG_MAX_TRACK];       /* F^(1) - b . (M + K) (nT) */
  double X1[MAG_MAX_TRACK];       /* X_sat - X_int - X_ext (nT) */
  double Y1[MAG_MAX_TRACK];       /* Y_sat - Y_int - Y_ext (nT) */
  double Z1[MAG_MAX_TRACK];       /* Z_sat - Z_int - Z_ext (nT) */
  double X2[MAG_MAX_TRACK];       /* X^(1) - (M + K)_X (nT) */
  double Y2[MAG_MAX_TRACK];       /* Y^(1) - (M + K)_Y (nT) */
  double Z2[MAG_MAX_TRACK];       /* Z^(1) - (M + K)_Z (nT) */

  double t_grad[MAG_MAX_TRACK];        /* timestamp for gradient point (CDF_EPOCH) */
  double theta_grad[MAG_MAX_TRACK];    /* colatitude for gradient point in radians */
  double phi_grad[MAG_MAX_TRACK];      /* longitude for gradient point in radians */
  double r_grad[MAG_MAX_TRACK];        /* radius for gradient point in km */
  double thetaq_grad[MAG_MAX_TRACK];   /* QD colatitude for gradient point in radians */
  double qdlat_grad[MAG_MAX_TRACK];    /* QD latitude for gradient point in degrees */
  double lat_deg_grad[MAG_MAX_TRACK];  /* geocentric latitude for gradient point in degrees */
  double F_grad[MAG_MAX_TRACK];        /* scalar field measurement for gradient point in nT */
  double X_grad[MAG_MAX_TRACK];        /* X field measurement for gradient point in nT */
  double Y_grad[MAG_MAX_TRACK];        /* Y field measurement for gradient point in nT */
  double Z_grad[MAG_MAX_TRACK];        /* Z field measurement for gradient point in nT */
  double Bx_int_grad[MAG_MAX_TRACK];   /* Bx_internal NEC for gradient point (nT) */
  double By_int_grad[MAG_MAX_TRACK];   /* By_internal NEC for gradient point (nT) */
  double Bz_int_grad[MAG_MAX_TRACK];   /* Bz_internal NEC for gradient point (nT) */
  double F_int_grad[MAG_MAX_TRACK];    /* F_internal for gradient point (nT) */
  double dF_ext_grad[MAG_MAX_TRACK];   /* b_int . B_ext for gradient point (nT) */
  double F1_grad[MAG_MAX_TRACK];       /* F_sat - F_int - dF_ext for gradient point (nT) */
  double F2_grad[MAG_MAX_TRACK];       /* F^(1) - b . (M + K) for gradient point (nT) */
  double X1_grad[MAG_MAX_TRACK];       /* X_sat - X_int - X_ext for gradient point (nT) */
  double Y1_grad[MAG_MAX_TRACK];       /* Y_sat - Y_int - Y_ext for gradient point (nT) */
  double Z1_grad[MAG_MAX_TRACK];       /* Z_sat - Z_int - Z_ext for gradient point (nT) */
  double X2_grad[MAG_MAX_TRACK];       /* X^(1) - (M + K)_X for gradient point (nT) */
  double Y2_grad[MAG_MAX_TRACK];       /* Y^(1) - (M + K)_Y for gradient point (nT) */
  double Z2_grad[MAG_MAX_TRACK];       /* Z^(1) - (M + K)_Z for gradient point (nT) */

  double Sq_int[MAG_MAX_TRACK];   /* internal Sq model b . M (nT) */
  double Sq_ext[MAG_MAX_TRACK];   /* external Sq model b . K (nT) */
  double X_Sq_int[MAG_MAX_TRACK]; /* internal Sq model M_x (nT) */
  double Y_Sq_int[MAG_MAX_TRACK]; /* internal Sq model M_y (nT) */
  double Z_Sq_int[MAG_MAX_TRACK]; /* internal Sq model M_z (nT) */
  double X_Sq_ext[MAG_MAX_TRACK]; /* external Sq model K_x (nT) */
  double Y_Sq_ext[MAG_MAX_TRACK]; /* external Sq model K_y (nT) */
  double Z_Sq_ext[MAG_MAX_TRACK]; /* external Sq model K_z (nT) */
  double F2_fit[MAG_MAX_TRACK];   /* fit to F^(2) from line current model */
  double X2_fit[MAG_MAX_TRACK];   /* fit to X^(2) from current model (nT) */
  double Y2_fit[MAG_MAX_TRACK];   /* fit to Y^(2) from current model (nT) */
  double Z2_fit[MAG_MAX_TRACK];   /* fit to Z^(2) from current model (nT) */
  double J_fit[MAG_MAX_TRACK];    /* current model fit to F^(2) data (A/km) */
  double phi_eq;                  /* longitude of equator crossing (rad) */
  double theta_eq;                /* geocentric colatitude of equator crossing (rad) */
  double t_eq;                    /* time of equator crossing (CDF_EPOCH) */
  size_t n;                       /* number of data in this track */
} mag_track;

typedef struct
{
  double qd_alt;    /* altitude at which QD map is made (km) */
  size_t nlon;      /* number of longitude segments per current arc */
  size_t p;         /* number of model parameters (number of line currents) */

  double qdlat_max; /* maximum QD latitude for arc currents */
  double dqdlat;    /* distance in deg between arc currents */
  double curr_dist_km; /* approx. distance in km between arc currents */

  gsl_matrix *Jx;   /* current density greens functions */
  gsl_matrix *Jy;
  gsl_matrix *Jz;

  /* ECEF Cartesian positions of all current segment midpoints (km) */
  gsl_matrix *mid_pos_x;
  gsl_matrix *mid_pos_y;
  gsl_matrix *mid_pos_z;

  double *mid_lon; /* approx. geocentric longitudes of segment midpts (deg) */

  gsl_matrix *X;     /* least squares matrix */
  gsl_matrix *cov;   /* covariance matrix */
  gsl_vector *S;     /* current strength coefficients (S_j) */
  gsl_vector *rhs;   /* right hand side vector (F^(2)) */
  gsl_matrix *L;     /* regularization matrix */
  gsl_vector *Ltau;  /* Householder vector for L decomposition */
  gsl_matrix *Xs;    /* least squares matrix in standard form */
  gsl_vector *ys;    /* right hand side vector in standard form */
  gsl_vector *cs;    /* solution vector in standard form */
  gsl_matrix *M;     /* matrix for standard to general form conversion */
  double rnorm;      /* residual norm || y - X c || */
  double snorm;      /* solution norm || L c || */
  gsl_multifit_linear_workspace *multifit_p;

  /* parameters for L-curve analysis */
  size_t nreg;       /* number of regularization parameters on L-curve */
  gsl_vector *rho;   /* vector of residual norms */
  gsl_vector *eta;   /* vector of solution norms */
  gsl_vector *reg_param; /* vector of regularization parameters */
  size_t reg_idx;    /* index of optimal regularization parameter */

  apex_workspace *apex_workspace_p;
  mageq_workspace *mageq_workspace_p;
} mag_eej_workspace;

typedef struct
{
  gsl_matrix *X;     /* least squares matrix */
  gsl_vector *c;     /* coefficient vector */
  gsl_vector *rhs;   /* right hand side vector */
  gsl_vector *L;     /* regularization matrix L */
  size_t p;          /* number of coefficients for Sq model */
  size_t n;          /* number of data for LS fit */
  size_t ntot;       /* total data size allocated */
  double rnorm;      /* residual norm || y - X c || */
  double snorm;      /* solution norm || L c || */

  size_t nmax_int;   /* maximum internal spherical harmonic degree */
  size_t mmax_int;   /* maximum internal spherical harmonic order */
  size_t nmax_ext;   /* maximum external spherical harmonic degree */
  size_t mmax_ext;   /* maximum external spherical harmonic order */
  size_t p_int;      /* number of internal coefficients */
  size_t p_ext;      /* number of external coefficients */
  size_t int_offset; /* offset in 'c' of internal coefficients */
  size_t ext_offset; /* offset in 'c' of external coefficients */

  /*
   * indexing for coefficients, base[n] gives the offset in c for
   * all coefficients of degree n
   */
  size_t *base_int;
  size_t *base_ext;

  gsl_vector *work_p; /* workspace of size p */

  size_t nreg;       /* number of regularization parameters for L-curve analysis */
  gsl_vector *rho;   /* vector of residual norms */
  gsl_vector *eta;   /* vector of solution norms */
  gsl_vector *reg_param; /* vector of regularization parameters */
  size_t reg_idx;    /* index of optimal regularization parameter */
  double Rsq;        /* coefficient of determination of fit */

  gsl_multifit_linear_workspace *multifit_workspace_p;
} mag_sqfilt_scalar_workspace;

typedef struct
{
  mag_params *params;

  mag_track track;           /* stores current track for processing */

  double *EEJ;               /* EEJ current density (A/m) */
  size_t ncurr;              /* number of line currents */

  double EEF;                /* final EEF value (mV/m) */
  double RelErr;             /* relative error output parameter */

  log_workspace *log_general;
  log_workspace *log_profile;
  log_workspace *log_F2;
  log_workspace *log_B2;
  log_workspace *log_B2_grad;
  log_workspace *log_Sq_Lcurve;
  log_workspace *log_Sq_Lcorner;
  log_workspace *log_Sq_svd;
  log_workspace *log_LC;
  log_workspace *log_EEJ;
  log_workspace *log_EEJ_Lcurve;
  log_workspace *log_EEJ_Lcorner;
  log_workspace *log_EEJ_svd;
  log_workspace *log_PDE;
  log_workspace *log_model;
  log_workspace *log_EEF;
  log_workspace *log_grids;
  log_workspace *log_fields;

  mag_sqfilt_scalar_workspace *sqfilt_scalar_workspace_p;
  mag_eej_workspace *eej_workspace_p;
  pde_workspace *pde_workspace_p;
  inverteef_workspace *inverteef_workspace_p;

  const char *output_file; /* output file */
  FILE *fp_output;

  double *dX;              /* internal X green's functions */
  double *dY;              /* internal Y green's functions */
  double *dZ;              /* internal Z green's functions */
  double *dX_ext;          /* external X green's functions */
  double *dY_ext;          /* external Y green's functions */
  double *dZ_ext;          /* external Z green's functions */

  green_workspace *green_int_p;
  green_workspace *green_ext_p;
  kp_workspace *kp_workspace_p;
  pomme_workspace *pomme_workspace_p;
  estist_calc_workspace *estist_calc_workspace_p;
  msynth_workspace *core_workspace_p;
  msynth_workspace *lith_workspace_p;
  magfit_workspace *magfit_workspace_p;
  magfit_workspace *magfit_eej_workspace_p;
} mag_workspace;

/*
 * Prototypes
 */

mag_workspace *mag_alloc(mag_params *params);
void mag_free(mag_workspace *w);
int mag_preproc(const mag_params *params, track_workspace *track_p,
                satdata_mag *data, mag_workspace *w);
int mag_proc(const mag_params *params, track_workspace *track_p,
             satdata_mag *data, mag_workspace *w);

/* mag_log.c */
int mag_log_profile(const int header, const size_t ntrack,
                    const double kp, const int dir,
                    const mag_workspace *w);
int mag_log_F2(const int header, const mag_workspace *w);
int mag_log_B2(const int header, const mag_workspace *w);
int mag_log_Sq_Lcurve(const int header, const mag_workspace *w);
int mag_log_Sq_Lcorner(const int header, const mag_workspace *w);
int mag_log_Sq_svd(const int header, const mag_workspace *w);
int mag_log_LC(const int header, const mag_workspace *w);
int mag_log_EEJ(const int header, const mag_workspace *w);
int mag_log_EEJ_Lcurve(const int header, const mag_workspace *w);
int mag_log_EEJ_Lcorner(const int header, const mag_workspace *w);
int mag_log_EEJ_svd(const int header, const mag_workspace *w);
int mag_log_PDE(const int header, const mag_workspace *w);
int mag_log_model(const int header, const mag_workspace *w);
int mag_log_EEF(const int header, const time_t t, const double phi,
                const double kp, const mag_workspace *w);
int mag_log_grids(const int header, const mag_workspace *w);
int mag_log_fields(const int header, const mag_workspace *w);

/* mag_log_grad.c */
int mag_log_B2_grad(const int header, const mag_workspace *w);

/* mag_sqfilt_scalar.c */
mag_sqfilt_scalar_workspace *mag_sqfilt_scalar_alloc(const size_t nmax_int, const size_t mmax_int,
                                                     const size_t nmax_ext, const size_t mmax_ext);
void mag_sqfilt_scalar_free(mag_sqfilt_scalar_workspace *w);
int mag_sqfilt_scalar(mag_workspace *mag_p, mag_sqfilt_scalar_workspace *w);

/* mag_sqfilt_vector.c */
int mag_sqfilt_vector(mag_workspace *mag_p, mag_sqfilt_scalar_workspace *w);

/* mag_eej.c */
mag_eej_workspace *mag_eej_alloc(const int year, const size_t ncurr,
                                 const double altitude,
                                 const double qdlat_max);
void mag_eej_free(mag_eej_workspace *w);
int mag_eej_proc_scalar(mag_track *track, gsl_vector *J, magfit_workspace *w);
int mag_eej_proc_vector(mag_track *track, gsl_vector *J, magfit_workspace *w);
int mag_eej_vector_proc(mag_track *track, double *J, mag_eej_workspace *w);

#endif /* INCLUDED_mag_h */
