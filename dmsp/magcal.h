/*
 * magcal.h
 */

#ifndef INCLUDED_magcal_h
#define INCLUDED_magcal_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit_nlinear.h>

#include <mainlib/ml_satdata.h>

/* number of fit parameters */
#define MAGCAL_P              9

#define MAGCAL_IDX_SX         0
#define MAGCAL_IDX_SY         1
#define MAGCAL_IDX_SZ         2
#define MAGCAL_IDX_OX         3
#define MAGCAL_IDX_OY         4
#define MAGCAL_IDX_OZ         5
#define MAGCAL_IDX_AXY        6
#define MAGCAL_IDX_AXZ        7
#define MAGCAL_IDX_AYZ        8

typedef struct
{
  size_t ntot;      /* total data allocated */
  size_t n;         /* number of data for current LS fit */
  size_t p;         /* number of calibration parameters */

  double *t;        /* timestamp array for interpolation (ms) */
  double *Ex;       /* original X measurement (nT) */
  double *Ey;       /* original Y measurement (nT) */
  double *Ez;       /* original Z measurement (nT) */

  double *F;        /* main field values (nT) */

  double lambda;    /* Tikhonov damping parameter */

  gsl_matrix *covar; /* parameter covariance matrix */

  gsl_multifit_nlinear_workspace *nlinear_workspace_p;
  gsl_spline *spline_f;
  gsl_interp_accel *acc;
} magcal_workspace;

typedef struct
{
  magcal_workspace *w;
  double dt; /* time shift (ms) */
} magcal_params;

/*
 * Prototypes
 */

magcal_workspace *magcal_alloc(const size_t n);
void magcal_free(magcal_workspace *w);
int magcal_add_datum(const double t, const double B_VFM[3], const double F, magcal_workspace *w);
int magcal_proc(gsl_vector *c, magcal_workspace *w);
double magcal_rms(const magcal_workspace * w);
time_t magcal_mean_time(const magcal_workspace * w);
int magcal_apply(const gsl_vector *m, satdata_mag *data);
int magcal_apply_cal(const gsl_vector *m, const double E[3], double B[4]);
int magcal_print_residuals(const char *filename, const gsl_vector * m, const magcal_workspace *w);

#endif /* INCLUDED_magcal_h */
