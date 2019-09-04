/*
 * magcal.h
 */

#ifndef INCLUDED_magcal_h
#define INCLUDED_magcal_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit_nlin.h>

#include <mainlib/ml_satdata.h>

/* number of fit parameters */
#define MAGCAL_P              9

#define MAGCAL_IDX_DT         0   /* time shift parameter */

/* scalar calibration parameter indices */
#define MAGCAL_IDX_SX         1   /* X scale factor */
#define MAGCAL_IDX_SY         2   /* Y scale factor */
#define MAGCAL_IDX_SZ         3   /* Z scale factor */
#define MAGCAL_IDX_OX         4   /* X offset */
#define MAGCAL_IDX_OY         5   /* Y offset */
#define MAGCAL_IDX_OZ         6   /* Z offset */
#define MAGCAL_IDX_AXY        7   /* a_{xy} non-orthogonality angle */
#define MAGCAL_IDX_AXZ        8   /* a_{xz} non-orthogonality angle */
#define MAGCAL_IDX_AYZ        9   /* a_{yz} non-orthogonality angle */

typedef struct
{
  size_t n;         /* number of data for current LS fit */
  size_t p;         /* number of calibration parameters */
  size_t n_max;     /* maximum number of data (size of arrays below) */

  double *Ex;       /* original X measurement (nT) */
  double *Ey;       /* original Y measurement (nT) */
  double *Ez;       /* original Z measurement (nT) */

  double *F;        /* scalar field data (dimensionless) */

  double B_s;       /* B scale factor (nT) */
  double lambda;    /* Tikhonov damping parameter */

  gsl_matrix *covar; /* parameter covariance matrix */

  const gsl_multifit_fdfsolver_type *fdf_type;
  gsl_multifit_fdfsolver *fdf_s;
  gsl_multifit_fdfridge *fdf_ridge;
} magcal_workspace;

typedef struct
{
  magcal_workspace *w;
} magcal_params;

/*
 * Prototypes
 */

magcal_workspace *magcal_alloc(const size_t n);
void magcal_free(magcal_workspace *w);
int magcal_initcond(gsl_vector *m);
int magcal_proc(gsl_vector *m, const satdata_mag *data, magcal_workspace *w);
int magcal_apply(const gsl_vector *m, satdata_mag *data);
int magcal_apply_cal(const gsl_vector *m, const double E[3], double B[4]);

#endif /* INCLUDED_magcal_h */
