/*
 * inverteef.h
 * Patrick Alken
 */

#ifndef INCLUDED_inverteef_h
#define INCLUDED_inverteef_h

#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>

/* define to invert for electric field */
#define PDE_INVERT_E_FIELD

/* define to allow a DC shift in inversion */
#define PDE_ALLOW_DC_SHIFT

/*
 * define to constrain inversion to have solution match satellite
 * profile at 90 degrees
 */
#define PDE_USE_LSE

typedef struct
{
  size_t ntheta;    /* number of theta grid points for PDE solution */
  double theta_min; /* minimum theta for PDE solution */
  double theta_max; /* maximum theta for PDE solution */
  size_t ncurr;     /* number of line currents in satellite profile */
  double qdlat_max; /* QD latitude range for line currents */
} inverteef_parameters;

typedef struct
{
  size_t ntheta;
  double theta_min;
  double theta_max;
  double dtheta;
  size_t ncurr;      /* number of line currents in satellite profile */
  double qdlat_max;  /* QD latitude range for line currents (degrees) */
  double qdlat_step; /* QD latitude step size for line currents (degrees) */

  gsl_vector *J_pde_E; /* pde profile for E_0, interpolated to line current grid */
  gsl_vector *J_pde_u; /* corresponding pde profile for winds, interpolated to line current grid */
  gsl_vector *J_pde;   /* total PDE modeled profile */
  gsl_vector *J_rhs;   /* RHS for inversion */
  gsl_vector *J_diff;  /* J_pde - J_sat */
  gsl_matrix *X;       /* least squares design matrix */
  gsl_vector *coeffs;  /* least squares coefficients */
  gsl_matrix *cov;     /* covariance matrix */
  gsl_matrix *B;       /* LSE constraint matrix */
  gsl_vector *d;       /* LSE rhs */
  gsl_vector *r;       /* fit residuals */
  double chisq;
  double Rsq;
  double R;
  gsl_multifit_linear_workspace *multifit_p;
  double E_scale;     /* E field scaling factor */
  double J_DC;        /* DC current shift */
  double J_champ;     /* CHAMP peak current value */

  gsl_interp_accel *acc;
  gsl_spline *spline_E;
  gsl_spline *spline_u;

  double RelErr;      /* relative error between model and satellite profile */
} inverteef_workspace;

/*
 * Prototypes
 */

int invert_calc(gsl_vector *J_sat, const gsl_matrix * J_phi);

inverteef_workspace *inverteef_alloc(inverteef_parameters *params);
void inverteef_free(inverteef_workspace *w);
int inverteef_calc(gsl_vector *J_sat, const double *qdlat_pde, gsl_vector *J_lat_E,
                   gsl_vector *J_lat_u, inverteef_workspace *w);

#endif /* INCLUDED_inverteef_h */
