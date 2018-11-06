/*
 * mfield_data.h
 */

#ifndef INCLUDED_mfield_data_h
#define INCLUDED_mfield_data_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

#include "magdata.h"

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
} mfield_data_parameters;

typedef struct
{
  size_t nsources;   /* number of data sources (satellites) */
  magdata **mdata;

  double *t0;        /* array of size nsources for first time of each satellite (CDF_EPOCH) */
  double *t1;        /* array of size nsources for last time of each satellite (CDF_EPOCH) */

  double t_mu;       /* mean of timestamps (years) */
  double t_sigma;    /* stddev of timestamps (years) */
  double t0_data;    /* timestamp of first data point (CDF_EPOCH) */
  double t1_data;    /* timestamp of last data point (CDF_EPOCH) */

  mfield_data_parameters params;

  gsl_rstat_workspace *rstat_workspace_p;
} mfield_data_workspace;

/*
 * Prototypes
 */

mfield_data_workspace *mfield_data_alloc(const size_t nsources,
                                         const mfield_data_parameters * params);
void mfield_data_free(mfield_data_workspace *w);
int mfield_data_copy(const size_t sat_idx, satdata_mag *data,
                     const size_t flags, mfield_data_workspace *w);
size_t mfield_data_filter_time(const double tmin, const double tmax,
                               mfield_data_workspace *w);
size_t mfield_data_filter_euler(mfield_data_workspace *w);
size_t mfield_data_filter_comp(mfield_data_workspace *w);
size_t mfield_data_filter_observatory(mfield_data_workspace *w);
int mfield_data_init(mfield_data_workspace *w);
double mfield_data_epoch(mfield_data_workspace *w);
int mfield_data_map(const char *dir_prefix, const mfield_data_workspace *w);
int mfield_data_print(const char *dir_prefix, const gsl_vector *wts_spatial, const mfield_data_workspace *w);
magdata *mfield_data_ptr(const size_t idx, const mfield_data_workspace *w);
int mfield_data_t(double *t0, double *t1, const magdata *data);
int mfield_data_add_noise(const double sigma, const double bias, mfield_data_workspace * w);

#endif /* INCLUDED_mfield_data_h */
