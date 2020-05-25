/*
 * invert_data.h
 */

#ifndef INCLUDED_invert_data_h
#define INCLUDED_invert_data_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

#include <mainlib/ml_magdata.h>

/* indices for each residual type */
#define INVERT_DATA_IDX_X              0
#define INVERT_DATA_IDX_Y              1
#define INVERT_DATA_IDX_Z              2
#define INVERT_DATA_IDX_F              3
#define INVERT_DATA_IDX_END            4

typedef struct
{
  double epoch;          /* model epoch in decimal years */
  double qdlat_fit_cutoff; /* QD latitude cutoff separating high-latitudes for fitting components (degrees) */

  /*
   * separate fitted components into mid/high latitudes; if 0, use the fit_xxx flags
   * below for all latitudes; if 1, the fit_xxx are for mid/low latitudes, and the fit_xxx_highlat
   * are for high latitudes
   */
  int fit_seplat;

  /* mid/low latitude components for fitting */

  int fit_X;             /* fit X vector component */
  int fit_Y;             /* fit Y vector component */
  int fit_Z;             /* fit Z vector component */
  int fit_F;             /* fit F scalar component */
  int fit_DXDT;          /* fit dX/dt vector component */
  int fit_DYDT;          /* fit dY/dt vector component */
  int fit_DZDT;          /* fit dZ/dt vector component */
  int fit_DX_NS;         /* fit DX N/S difference component */
  int fit_DY_NS;         /* fit DY N/S difference component */
  int fit_DZ_NS;         /* fit DZ N/S difference component */
  int fit_DF_NS;         /* fit DF N/S difference component */
  int fit_DX_EW;         /* fit DX E/W difference component */
  int fit_DY_EW;         /* fit DY E/W difference component */
  int fit_DZ_EW;         /* fit DZ E/W difference component */
  int fit_DF_EW;         /* fit DF E/W difference component */

  /* high latitude components for fitting */

  int fit_Z_highlat;     /* fit high-latitude Z vector component */
  int fit_F_highlat;     /* fit high-latitude F scalar component */
  int fit_DZ_NS_highlat; /* fit high-latitude DZ N/S difference component */
  int fit_DF_NS_highlat; /* fit high-latitude DF N/S difference component */
  int fit_DZ_EW_highlat; /* fit high-latitude DZ E/W difference component */
  int fit_DF_EW_highlat; /* fit high-latitude DF E/W difference component */
} invert_data_parameters;

typedef struct
{
  size_t nsources;   /* number of data sources (satellites) */
  magdata **mdata;

  double *t0;        /* array of size nsources for first time of each satellite (CDF_EPOCH) */
  double *t1;        /* array of size nsources for last time of each satellite (CDF_EPOCH) */

  double t0_data;    /* timestamp of first data point (CDF_EPOCH) */
  double t1_data;    /* timestamp of last data point (CDF_EPOCH) */

  invert_data_parameters params;

  gsl_rstat_workspace *rstat_workspace_p;
} invert_data_workspace;

/*
 * Prototypes
 */

invert_data_workspace *invert_data_alloc(const size_t nsources,
                                         const invert_data_parameters * params);
void invert_data_free(invert_data_workspace *w);
int invert_data_copy(const size_t sat_idx, satdata_mag *data,
                     const size_t flags, invert_data_workspace *w);
size_t invert_data_filter_time(const double tmin, const double tmax,
                               invert_data_workspace *w);
size_t invert_data_filter_comp(invert_data_workspace *w);
size_t invert_data_filter_observatory(invert_data_workspace *w);
int invert_data_compact(invert_data_workspace * w);
int invert_data_init(invert_data_workspace *w);
int invert_data_map(const char *dir_prefix, const invert_data_workspace *w);
int invert_data_print(const char *dir_prefix, const gsl_vector *wts_spatial, const invert_data_workspace *w);
magdata *invert_data_ptr(const size_t idx, const invert_data_workspace *w);
int invert_data_t(double *t0, double *t1, const magdata *data);
int invert_data_weights(gsl_vector * wts, const double weightfac[], const invert_data_workspace * w);

#endif /* INCLUDED_invert_data_h */
