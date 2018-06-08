/*
 * peak.c
 *
 * Routines for finding peaks in time series data. Calling sequence:
 *
 * 1. peak_alloc    - allocate peak workspace
 * 2. peak_init     - initialize by computing and smoothing first derivative
 * 3. peak_find     - find peaks by looking for zero crossings of first derivative
 * 4. peak_gaussian - fit gaussian to detected peak
 * 5. peak_free     - free memory
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_filter.h>

#include "peak.h"

/*
peak_alloc()
  Allocate peak workspace

Inputs: n       - maximum length of input vector
        K_gauss - size of window for Gaussian filter
        K_scale - size of window for scale estimate
*/

peak_workspace *
peak_alloc(const size_t n, const size_t K_gauss, const size_t K_scale)
{
  peak_workspace *w;

  w = calloc(1, sizeof(peak_workspace));
  if (!w)
    return 0;

  w->n = n;

  w->deriv = gsl_vector_alloc(n);
  w->median = gsl_vector_alloc(n);
  w->sigma = gsl_vector_alloc(n);
  w->work = gsl_vector_alloc(n);
  w->ioutlier = gsl_vector_int_alloc(n);
  w->gaussian_workspace_p = gsl_filter_gaussian_alloc(K_gauss);
  w->movstat_workspace_p = gsl_movstat_alloc(K_scale);
  w->impulse_workspace_p = gsl_filter_impulse_alloc(K_scale);

  return w;
}

void
peak_free(peak_workspace *w)
{
  if (w->deriv)
    gsl_vector_free(w->deriv);

  if (w->median)
    gsl_vector_free(w->median);

  if (w->sigma)
    gsl_vector_free(w->sigma);

  if (w->work)
    gsl_vector_free(w->work);

  if (w->ioutlier)
    gsl_vector_int_free(w->ioutlier);

  if (w->gaussian_workspace_p)
    gsl_filter_gaussian_free(w->gaussian_workspace_p);

  if (w->impulse_workspace_p)
    gsl_filter_impulse_free(w->impulse_workspace_p);

  if (w->movstat_workspace_p)
    gsl_movstat_free(w->movstat_workspace_p);

  free(w);
}

/*
peak_find()

Inputs: minheight - minimum height needed to detect peak; set to
                    0.0 to disable this test
        nsigma    - number of standard deviations a peak must be from its
                    local median
        x         - input data
        npeak     - (output) number of peaks found
        ipeak     - (output) boolean vector of peak locations
                    ipeak[j] = -1 if minimum peak found at location j
                             =  0 if no peak found at location j
                             = +1 if maximum peak found at location j
        w         - workspace

Return:
*/

int
peak_find(const double minheight, const double nsigma, const gsl_vector * x, size_t * npeak, gsl_vector_int * ipeak, peak_workspace *w)
{
  const size_t n = x->size;

  if (n > w->n)
    {
      GSL_ERROR("x vector does not match workspace", GSL_EBADLEN);
    }
  else if (ipeak->size != n)
    {
      GSL_ERROR("ipeak vector does not match x size", GSL_EBADLEN);
    }
  else
    {
#if 0
      const double sigma = 0.25;
      gsl_vector_view d = gsl_vector_subvector(w->deriv, 0, n);
      gsl_vector_view med = gsl_vector_subvector(w->median, 0, n);
      gsl_vector_view scale = gsl_vector_subvector(w->sigma, 0, n);
      gsl_vector_int_view iout = gsl_vector_int_subvector(w->ioutlier, 0, n);
      gsl_vector_view y = gsl_vector_subvector(w->work, 0, n);
      size_t nimpulse;
      size_t i;

      /* initialize */
      gsl_vector_int_set_zero(ipeak);
      *npeak = 0;

      /* apply smoothing (Gaussian) filter to input and differentiate at same time */
      gsl_filter_gaussian(sigma, 1, x, &d.vector, w->gaussian_workspace_p);

      /* apply impulse detection filter to signal to identify candidate peaks */
      gsl_filter_impulse(GSL_FILTER_END_PADVALUE, GSL_FILTER_SCALE_MAD, nsigma, x, &y.vector,
                         &med.vector, &scale.vector, &nimpulse, &iout.vector, w->impulse_workspace_p);

      for (i = 0; i < n - 1; ++i)
        {
          double di = gsl_vector_get(&d.vector, i);
          double dip1 = gsl_vector_get(&d.vector, i + 1);
          int iouti = gsl_vector_int_get(&iout.vector, i);

          /* look for an outlier (sample is nsigma stddevs above surrounding noise) and also
           * a zero crossing of first derivative indicating a maxima/minima */
          if (iouti && di * dip1 < 0.0)
            {
              double xi = gsl_vector_get(x, i);
              double xip1 = gsl_vector_get(x, i + 1);

              /* ensure peak height is above minimum threshold */
              if (fabs(xi) < minheight || fabs(xip1) < minheight)
                continue;

              if (di > 0.0)
                gsl_vector_int_set(ipeak, i, 1);  /* peak */
              else
                gsl_vector_int_set(ipeak, i, -1); /* valley */

              ++(*npeak);
            }

          printf("%zu %.12e %.12e %.12e %.12e %d\n",
                 i,
                 gsl_vector_get(x, i),
                 di,
                 gsl_vector_get(&med.vector, i),
                 gsl_vector_get(&scale.vector, i),
                 gsl_vector_int_get(ipeak, i));
        }
#endif

      return GSL_SUCCESS;
    }
}
