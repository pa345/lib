/*
 * peak.h
 */

#ifndef INCLUDED_peak_h
#define INCLUDED_peak_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_filter.h>

typedef struct
{
  gsl_vector *deriv;        /* first derivative of input data */
  gsl_vector *median;       /* median of window */
  gsl_vector *sigma;        /* dispersion of window */
  gsl_vector *work;         /* workspace, size n */
  gsl_vector_int *ioutlier; /* outlier detected */
  size_t n;                 /* number of data points to process */
  gsl_filter_gaussian_workspace *gaussian_workspace_p;
  gsl_movstat_workspace *movstat_workspace_p;
  gsl_filter_impulse_workspace *impulse_workspace_p;
} peak_workspace;

/*
 * Prototypes
 */

peak_workspace *peak_alloc(const size_t n, const size_t K_gauss, const size_t K_scale);
void peak_free(peak_workspace *w);
int peak_find(const double minheight, const double nsigma, const gsl_vector * x, size_t * npeak,
              gsl_vector_int * ipeak, peak_workspace *w);

#endif /* INCLUDED_peak_h */
