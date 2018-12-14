/*
 * invert.c
 *
 * Invert PDE solution vectors for EEF and wind perturbation field
 *
 * Steps:
 *
 * 1. Interpolate PDE solutions J_lat_E and J_lat_u to the same grid
 *    used for line currents
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_spline.h>

#include <indices/indices.h>
#include <common/common.h>

#include "inverteef.h"

#include "lse.c"

/*
invert_calc()
  Invert satellite profile

Inputs: J_sat     - satellite current vector profile (length ncurr)
        J_phi     - J_phi matrix (height-integrated), ntheta-by-nrhs

Notes:
1) On output, w->J_pde contains the modeled profile
2) On output, w->RelErr contains relative error between modeled and
satellite profile
*/

int
invert_calc(gsl_vector *J_sat, const gsl_matrix * J_phi)
{
  const size_t ncurr = J_sat->size;
  const size_t ntheta = J_phi->size1;
  const size_t nrhs = J_phi->size2;
  size_t i;

  for (i = 0; i < ncurr; ++i)
    printf("%.12e\n", gsl_vector_get(J_sat, i));

  printf("\n\n");
  for (i = 0; i < ntheta; ++i)
    printf("%.12e %.12e\n", gsl_matrix_get(J_phi, i, 0), gsl_matrix_get(J_phi, i, 1));

#if 0
  fprintf(stderr, "inverteef_calc: interpolating PDE solution...");

  inverteef_interp_pdesol(J_sat, qdlat_pde, J_lat_E, J_lat_u, w);

  fprintf(stderr, "done\n");
  fflush(stderr);
  fprintf(stderr, "inverteef_calc: constructing least squares matrix...");

  inverteef_design_matrix(w);

  fprintf(stderr, "done\n");
  fflush(stderr);
  fprintf(stderr, "inverteef_calc: inverting profile...");

  inverteef_invert_profile(J_sat, w);

  /* compute final modeled profile and store in w->J_pde */
  {
    size_t j;
    double x[2];
    gsl_vector_view xv = gsl_vector_view_array(x, 2);

    for (j = 0; j < w->ncurr; ++j)
      {
        double J, Jerr;

        x[0] = gsl_vector_get(w->J_pde_E, j);
        x[1] = -1.0;

        gsl_multifit_linear_est(&xv.vector, w->coeffs, w->cov,
                                &J, &Jerr);

        J += gsl_vector_get(w->J_pde_u, j);

        gsl_vector_set(w->J_pde, j, J);
        gsl_vector_set(w->J_diff, j, J - gsl_vector_get(J_sat, j));
      }
  }

  w->R = gsl_stats_correlation(J_sat->data, 1,
                               w->J_pde->data, 1,
                               w->ncurr);

  /* compute relative error between modeled and satellite profiles */
  w->RelErr = gsl_blas_dnrm2(w->J_diff) / gsl_blas_dnrm2(J_sat);

  fprintf(stderr, "done (chisq = %e, Rsq = %f, R = %f)\n",
          w->chisq, w->Rsq, w->R);
  fprintf(stderr, "inverteef_calc: scale factor = %f\n", w->E_scale);
  fflush(stderr);
#endif

  return GSL_SUCCESS;
}
