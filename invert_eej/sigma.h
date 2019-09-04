/*
 * sigma.h
 */

#ifndef INCLUDED_sigma_h
#define INCLUDED_sigma_h

#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#include <mainlib/ml_mageq.h>

#include "cond.h"

/* force latitude-symmetric conductivities */
#define SIGMA_SYMMETRIC

typedef struct
{
  size_t nr;         /* number of radial grid points */
  size_t ntheta;     /* number of theta grid points */

  double rmin;       /* lower boundary in m */
  double rmax;       /* upper boundary in m */
  double rstep;      /* radial step size in m */

  double theta_min;  /* minimum theta in radians */
  double theta_max;  /* maximum theta in radians */
  double theta_step; /* theta step size in radians */

  gsl_matrix *s0;    /* direct conductivity */
  gsl_matrix *s1;    /* Pedersen conductivity */
  gsl_matrix *s2;    /* Hall conductivity */

  cond_workspace *cond_workspace_p;
  mageq_workspace *mageq_workspace_p;
} sigma_workspace;

/*
 * Prototypes
 */

sigma_workspace *sigma_alloc(const char *f107_file, size_t nr, size_t ntheta, double rmin,
                             double rmax, double theta_min,
                             double theta_max);
void sigma_free(sigma_workspace *w);
int sigma_calc(time_t t, double longitude, sigma_workspace *w);
int sigma_result(size_t i, size_t j, double *s0, double *s1, double *s2,
                 sigma_workspace *w);
int sigma_max(double *s0_max, double *s1_max, double *s2_max,
              sigma_workspace *w);

#endif /* INCLUDED_sigma_h */
