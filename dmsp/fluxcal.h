/*
 * fluxcal.h
 */

#ifndef INCLUDED_fluxcal_h
#define INCLUDED_fluxcal_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit_nlinear.h>

#include <mainlib/ml_satdata.h>

/* number of fit parameters */
#define FLUXCAL_P              9

#define FLUXCAL_IDX_SX         0
#define FLUXCAL_IDX_SY         1
#define FLUXCAL_IDX_SZ         2
#define FLUXCAL_IDX_OX         3
#define FLUXCAL_IDX_OY         4
#define FLUXCAL_IDX_OZ         5
#define FLUXCAL_IDX_U1         6
#define FLUXCAL_IDX_U2         7
#define FLUXCAL_IDX_U3         8

typedef struct
{
  size_t ntot;      /* total data allocated */
  size_t n;         /* number of data for current LS fit */
  size_t p;         /* number of calibration parameters */

  double *t;        /* timestamp array for interpolation (CDF_EPOCH) */
  double *Ex;       /* original X measurement (nT) */
  double *Ey;       /* original Y measurement (nT) */
  double *Ez;       /* original Z measurement (nT) */

  double *F;        /* main field values (nT) */

  gsl_vector *weights; /* robust weights */
  double lambda;    /* Tikhonov damping parameter */
  size_t max_iter;  /* maximum nls iterations */

  gsl_matrix *covar; /* parameter covariance matrix */

  gsl_multifit_nlinear_workspace *nlinear_workspace_p;
} fluxcal_workspace;

typedef struct
{
  fluxcal_workspace *w;
  double dt; /* time shift (ms) */
} fluxcal_params;

/*
 * Prototypes
 */

fluxcal_workspace *fluxcal_alloc(const size_t n);
void fluxcal_free(fluxcal_workspace *w);
int fluxcal_add_datum(const double t, const double B_VFM[3], const double F, fluxcal_workspace *w);
int fluxcal_proc(gsl_vector *c, fluxcal_workspace *w);
int fluxcal_nls(const gsl_vector * weights, gsl_vector * c, fluxcal_workspace *w);
double fluxcal_rms(const fluxcal_workspace * w);
time_t fluxcal_mean_time(const fluxcal_workspace * w);
int fluxcal_apply(const gsl_vector *m, satdata_mag *data);
int fluxcal_apply_datum(const gsl_vector *m, const double E[3], double B[4]);
int fluxcal_print_residuals(const char *filename, const gsl_vector * m, const fluxcal_workspace *w);

#endif /* INCLUDED_fluxcal_h */
