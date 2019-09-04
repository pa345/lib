/*
 * sigma.c
 *
 * This file contains routines related to computing the conductivity
 * tensor for each grid point
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <setjmp.h>
#include <signal.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_mageq.h>

#include "cond.h"

#include "sigma.h"

/*
 * Global
 */
static int cond_success = 0;
static jmp_buf exception_env;

static void
pde_sighandler(int sig)
{
  if (sig == SIGALRM)
    {
      if (cond_success == 0)
        {
          /*
           * The IRI call failed - jump back to main loop and continue
           * processing profiles
           */
          siglongjmp(exception_env, 1);
        }
    }
} /* pde_sighandler() */

sigma_workspace *
sigma_alloc(const char *f107_file, size_t nr, size_t ntheta,
            double rmin, double rmax, double theta_min, double theta_max)
{
  sigma_workspace *w;

  w = calloc(1, sizeof(sigma_workspace));
  if (!w)
    {
      fprintf(stderr, "sigma_alloc: calloc failed: %s\n", strerror(errno));
      return 0;
    }

  w->s0 = gsl_matrix_alloc(nr, ntheta);
  w->s1 = gsl_matrix_alloc(nr, ntheta);
  w->s2 = gsl_matrix_alloc(nr, ntheta);
  if (!w->s0 || !w->s1 || !w->s2)
    {
      sigma_free(w);
      return 0;
    }

  w->cond_workspace_p = cond_alloc(nr, f107_file);
  if (!w->cond_workspace_p)
    {
      sigma_free(w);
      return 0;
    }

  w->mageq_workspace_p = mageq_alloc();
  if (!w->mageq_workspace_p)
    {
      sigma_free(w);
      return 0;
    }

  cond_set_alpha(4.0, w->cond_workspace_p);

  /*cond_set_error_scale(1.0, 4.0, w->cond_workspace_p);*/

  w->nr = nr;
  w->ntheta = ntheta;

  w->rmin = rmin;
  w->rmax = rmax;
  w->rstep = (rmax - rmin) / (nr - 0);

  w->theta_min = theta_min;
  w->theta_max = theta_max;
  w->theta_step = (theta_max - theta_min) / (ntheta - 0);

  return w;
} /* sigma_alloc() */

void
sigma_free(sigma_workspace *w)
{
  if (w->s0)
    gsl_matrix_free(w->s0);

  if (w->s1)
    gsl_matrix_free(w->s1);

  if (w->s2)
    gsl_matrix_free(w->s2);

  if (w->cond_workspace_p)
    cond_free(w->cond_workspace_p);

  if (w->mageq_workspace_p)
    mageq_free(w->mageq_workspace_p);

  free(w);
}

/*
pde_sigma()
  Compute conductivity tensor for entire grid

Inputs: t         - timestamp
        longitude - longitude (radians)
        w         - workspace

Return: success or failure

Notes:

1) w->sigma is filled in with the conductivity tensor

2) sigma = (sigma_0 - sigma_1) b b^t + sigma_1 * I_3 + sigma_2 * A_2

where:

sigma_0 = parallel conductivity
sigma_1 = Pederson conductivity
sigma_2 = Hall conductivity
b = unit magnetic field vector
A_2 = (   0  -b_3  b_2 )
      (  b_3   0  -b_1 )
      ( -b_2  b_1   0  )
*/

int
sigma_calc(time_t t, double longitude, sigma_workspace *w)
{
  size_t i, j;
  int s;
  double lat_eq; /* latitude of magnetic equator for this longitude */
  const double r = 6371.2 + 108.0;
  const double tyr = get_year(t);

  lat_eq = mageq_calc(longitude, r, tyr, w->mageq_workspace_p);

  /* compute conductivities */
  for (j = 0; j < w->ntheta; ++j)
    {
      /* colatitude offset */
      double theta_off = w->theta_min + j * w->theta_step;
      double theta = theta_off - lat_eq;

      /*
       * set an alarm for 60 seconds so if IRI fails we can continue
       * processing new profiles
       */
      cond_success = 0;
      /*alarm(60);*/

      s = cond_calc(theta,
                    longitude,
                    t,
                    w->rmin * 1.0e-3 - R_EARTH_KM,
                    w->rstep * 1.0e-3,
                    w->nr,
                    w->cond_workspace_p);

      alarm(0);
      cond_success = 1;

      if (s)
        return GSL_FAILURE; /* error occurred */
        
      /* save conductivity results */
      for (i = 0; i < w->nr; ++i)
        {
          cond_result *result = cond_get_result(i, w->cond_workspace_p);

          if (!gsl_finite(result->sigma_0) ||
              !gsl_finite(result->sigma_p) ||
              !gsl_finite(result->sigma_h))
            return GSL_FAILURE;

          gsl_matrix_set(w->s0, i, j, result->sigma_0);
          gsl_matrix_set(w->s1, i, j, result->sigma_p);
          gsl_matrix_set(w->s2, i, j, result->sigma_h);
        } /* for (i = 0; i < w->nr; ++i) */
    } /* for (j = 0; j < w->ntheta; ++j) */

#ifdef SIGMA_SYMMETRIC

  for (j = 0; j < w->ntheta / 2; ++j)
    {
      for (i = 0; i < w->nr; ++i)
        {
          double tmp;

          tmp = 0.5 * (gsl_matrix_get(w->s0, i, j) +
                       gsl_matrix_get(w->s0, i, w->ntheta - j - 1));
          gsl_matrix_set(w->s0, i, j, tmp);
          gsl_matrix_set(w->s0, i, w->ntheta - j - 1, tmp);

          tmp = 0.5 * (gsl_matrix_get(w->s1, i, j) +
                       gsl_matrix_get(w->s1, i, w->ntheta - j - 1));
          gsl_matrix_set(w->s1, i, j, tmp);
          gsl_matrix_set(w->s1, i, w->ntheta - j - 1, tmp);

          tmp = 0.5 * (gsl_matrix_get(w->s2, i, j) +
                       gsl_matrix_get(w->s2, i, w->ntheta - j - 1));
          gsl_matrix_set(w->s2, i, j, tmp);
          gsl_matrix_set(w->s2, i, w->ntheta - j - 1, tmp);
        }
    }

#endif /* SIGMA_SYMMETRIC */

  return GSL_SUCCESS;
} /* sigma_calc() */

/*
sigma_result()
  Access the conductivities computed previously by pde_sigma()

Inputs: i  - radial grid point in [0,nr-1]
        j  - theta grid point in [0,ntheta-1]
        s0 - (output) where to store direct conductivity
        s1 - (output) where to store Pedersen conductivity
        s2 - (output) where to store Hall conductivity
        w  - workspace
*/

int
sigma_result(size_t i, size_t j, double *s0, double *s1, double *s2,
             sigma_workspace *w)
{
  int s = 0;

  *s0 = gsl_matrix_get(w->s0, i, j);
  *s1 = gsl_matrix_get(w->s1, i, j);
  *s2 = gsl_matrix_get(w->s2, i, j);

  return s;
} /* sigma_result() */

int
sigma_max(double *s0_max, double *s1_max, double *s2_max,
          sigma_workspace *w)
{
  int s = 0;

  *s0_max = gsl_matrix_max(w->s0);
  *s1_max = gsl_matrix_max(w->s1);
  *s2_max = gsl_matrix_max(w->s2);

  return s;
} /* sigma_max() */
