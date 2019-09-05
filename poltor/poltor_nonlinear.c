/*
 * poltor_nonlinear.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multilarge_nlinear.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_green.h>
#include <mainlib/ml_lapack_wrapper.h>

#include "poltor.h"

static int poltor_nonlinear_jac(const gsl_vector_complex *x, const gsl_vector *weights, poltor_workspace *w);
static int poltor_nonlinear_regularize(poltor_workspace *w);
static int poltor_calc_f(gsl_vector_complex * x, void * params, gsl_vector * f);
static int poltor_robust_print_stat(const char *str, const double sigma, const gsl_rstat_workspace *rstat_p);
static int poltor_robust_weights(const gsl_vector * f, gsl_vector * weights, poltor_workspace * w);
static double huber(const double x);
static double bisquare(const double x);

/*
poltor_calc_nonlinear()
  Perform one least squares iteration

Inputs: iter - iteration number
        c    - model coefficients; on input, c contains a best estimate for
               the coefficient vector, which is replaced on output by an
               updated estimate, size p
        w    - workspace
*/

int
poltor_calc_nonlinear(const size_t iter, gsl_vector_complex *c, poltor_workspace *w)
{
  int s = 0;
  const poltor_parameters *params = &(w->params);
  struct timeval tv0, tv1;
  double rcond;

  if (params->use_weights)
    {
      gsl_vector_memcpy(w->wts_final, w->wts_spatial);

      if (iter > 1)
        {
          /* calculate residuals with coefficients from previous iteration */
          fprintf(stderr, "poltor_calc_nonlinear: computing residuals with previous coefficients...");
          gettimeofday(&tv0, NULL);
          poltor_calc_f(c, w, w->f);
          gettimeofday(&tv1, NULL);
          fprintf(stderr, "done (%g seconds, ||f|| = %.12e)\n", time_diff(tv0, tv1), gsl_blas_dnrm2(w->f));

          poltor_robust_weights(w->f, w->wts_robust, w);

          /* final weights = spatial * robust */
          gsl_vector_mul(w->wts_final, w->wts_robust);
        }
    }
  else
    {
      gsl_vector_set_all(w->wts_final, 1.0);
    }

  if (w->lls_solution == 1)
    {
      /* no scalar residuals so it is a linear problem */

      /* compute J^H W J and J^H W b (set x = 0) */
      fprintf(stderr, "poltor_calc_nonlinear: computing J^H W J and J^H W f [LINEAR CASE]...");
      gettimeofday(&tv0, NULL);
      poltor_nonlinear_jac(NULL, w->wts_final, w);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      if (params->regularize)
        {
          fprintf(stderr, "poltor_calc_nonlinear: regularizing system...");
          poltor_nonlinear_regularize(w);
          fprintf(stderr, "done\n");
        }

      fprintf(stderr, "poltor_calc_nonlinear: solving linear system J^H W J = J^H W f...");
      s = lapack_complex_zposv(w->JHf, w->JHJ, c, &rcond);
      fprintf(stderr, "done (s = %d, cond = %.12e)\n", s, 1.0 / rcond);

      /* reverse sign since precompute calculates -J^H W b */
      gsl_vector_complex_scale(c, GSL_COMPLEX_NEGONE);
    }
  else
    {
      const gsl_complex step = gsl_complex_rect(-0.7, 0.0);
      gsl_vector_complex *delta = gsl_vector_complex_alloc(w->p);

      /* compute J^T W J and J^T W f, with f = model - data */
      fprintf(stderr, "poltor_calc_nonlinear: computing J^H W J and J^H W f [NONLINEAR CASE]...");
      gettimeofday(&tv0, NULL);
      poltor_nonlinear_jac(c, w->wts_final, w);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      if (params->regularize)
        {
          fprintf(stderr, "poltor_calc_nonlinear: regularizing system...");
          poltor_nonlinear_regularize(w);
          fprintf(stderr, "done\n");
        }

      fprintf(stderr, "poltor_calc_nonlinear: solving linear system J^H W J = J^H W f for step delta...");
      s = lapack_zposv(w->JHJ, w->JHf, delta, &rcond);
      fprintf(stderr, "done (s = %d, cond = %.12e)\n", s, 1.0 / rcond);

      /* c += step * delta */
      gsl_blas_zaxpy(step, delta, c);

      gsl_vector_complex_free(delta);
    }

  return s;
}

/*
poltor_nonlinear_jac()
  Precompute J_int^T W J_int for vector measurements, since
this submatrix is independent of the model parameters and
only needs to be computed once per iteration. This function
uses OpenMP to speed up the calculation

Inputs: x       - model parameter vector (may be NULL in which case only
                  vector measurements are used (which don't depend on x))
        weights - weight vector
        w       - workspace

Notes:
1) w->JTJ_vec is updated with J_int^T W J_int for vector
residuals
*/

static int
poltor_nonlinear_jac(const gsl_vector_complex *x, const gsl_vector *weights, poltor_workspace *w)
{
  int s = GSL_SUCCESS;
  const magdata_list *list = w->data;
  size_t i;

  gsl_matrix_complex_set_zero(w->JHJ);
  gsl_vector_complex_set_zero(w->JHf);

  /*
   * omp_rowidx[thread_id] contains the number of currently filled rows
   * of omp_J[thread_id]. When omp_J[thread_id] is filled, it is folded
   * into the JTJ_vec matrix and then omp_rowidx[thread_id] is reset to 0
   */
  for (i = 0; i < w->max_threads; ++i)
    {
      w->omp_nrows[i] = 0;
      w->omp_rowidx[i] = 0;
    }

  fprintf(stderr, "\n");

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);
      size_t j;

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index of this data point in [0,w->n-1] */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          gsl_complex fval;
          double B_prior[4], B_prior_grad[4]; /* a priori field models (main, crust, magnetosphere) */
          gsl_complex B_model[3];       /* B_model (ionosphere + prior) */
          gsl_complex B_model_grad[3];  /* B_model (ionosphere + prior) for gradient point */
          double F_model, F_model_grad; /* || B_model ||, || B_model_grad || */
          gsl_complex b_model[3];       /* B_model / || B_model || */
          gsl_complex b_model_grad[3];  /* B_model_grad / || B_model_grad || */
          gsl_vector_complex_view vb_model = gsl_vector_complex_view_array((double *) b_model, 3);
          gsl_vector_complex_view vb_model_grad = gsl_vector_complex_view_array((double *) b_model_grad, 3);
          double alpha = gsl_vector_get(w->solar_flux, ridx);           /* alpha factor = 1 + N*EUVAC (same for all residuals for this data point) */
          double alpha_grad = gsl_vector_get(w->solar_flux_grad, ridx); /* alpha for gradient point */
          size_t k;

          GSL_SET_IMAG(&fval, 0.0);

          gsl_vector_complex_view vx_int, vy_int, vz_int;
          gsl_vector_complex_view vx_ext, vy_ext, vz_ext;
          gsl_vector_complex_view vx_sh, vy_sh, vz_sh;
          gsl_vector_complex_view vx_tor, vy_tor, vz_tor;
          gsl_vector_complex_view vx_grad_int, vy_grad_int, vz_grad_int;
          gsl_vector_complex_view vx_grad_ext, vy_grad_ext, vz_grad_ext;
          gsl_vector_complex_view vx_grad_sh, vy_grad_sh, vz_grad_sh;
          gsl_vector_complex_view vx_grad_tor, vy_grad_tor, vz_grad_tor;
          gsl_vector_complex_view v;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute prior model (B_main + B_crust + B_ext) */
          magdata_prior(j, B_prior, mptr);

          /* compute Y_{nm} = P_{nm}(cos(theta)) exp(i m phi) and d/dtheta Y_{nm} functions,
           * stored in w->green_p[thread_id]->{Ynm,dYnm} */
          green_complex_Ynm_deriv(theta, phi, w->green_p[thread_id]);

          if (w->p_pint > 0)
            {
              vx_int = gsl_matrix_complex_subrow(w->omp_dX, thread_id, w->pint_offset, w->p_pint);
              vy_int = gsl_matrix_complex_subrow(w->omp_dY, thread_id, w->pint_offset, w->p_pint);
              vz_int = gsl_matrix_complex_subrow(w->omp_dZ, thread_id, w->pint_offset, w->p_pint);

              /* calculate internal Green's functions */
              green_complex_calc_int(w->nmax_int, w->mmax_int, w->R, r, theta,
                                     w->green_p[thread_id]->Ynm,
                                     w->green_p[thread_id]->dYnm,
                                     (complex double *) vx_int.vector.data,
                                     (complex double *) vy_int.vector.data,
                                     (complex double *) vz_int.vector.data);
            }

          if (w->p_pext > 0)
            {
              vx_ext = gsl_matrix_complex_subrow(w->omp_dX, thread_id, w->pext_offset, w->p_pext);
              vy_ext = gsl_matrix_complex_subrow(w->omp_dY, thread_id, w->pext_offset, w->p_pext);
              vz_ext = gsl_matrix_complex_subrow(w->omp_dZ, thread_id, w->pext_offset, w->p_pext);

              /* calculate external Green's functions */
              green_complex_calc_ext(w->nmax_ext, w->mmax_ext, w->R, r, theta,
                                     w->green_p[thread_id]->Ynm,
                                     w->green_p[thread_id]->dYnm,
                                     (complex double *) vx_ext.vector.data,
                                     (complex double *) vy_ext.vector.data,
                                     (complex double *) vz_ext.vector.data);
            }

          if (w->p_psh > 0)
            {
              vx_sh = gsl_matrix_complex_subrow(w->omp_dX, thread_id, w->psh_offset, w->p_psh);
              vy_sh = gsl_matrix_complex_subrow(w->omp_dY, thread_id, w->psh_offset, w->p_psh);
              vz_sh = gsl_matrix_complex_subrow(w->omp_dZ, thread_id, w->psh_offset, w->p_psh);

              /* calculate shell Green's functions */
              poltor_calc_psh(r, theta,
                              w->green_p[thread_id]->Ynm,
                              w->green_p[thread_id]->dYnm,
                              (complex double *) vx_sh.vector.data,
                              (complex double *) vy_sh.vector.data,
                              (complex double *) vz_sh.vector.data,
                              w);
            }

          if (w->p_tor > 0)
            {
              vx_tor = gsl_matrix_complex_subrow(w->omp_dX, thread_id, w->tor_offset, w->p_tor);
              vy_tor = gsl_matrix_complex_subrow(w->omp_dY, thread_id, w->tor_offset, w->p_tor);
              vz_tor = gsl_matrix_complex_subrow(w->omp_dZ, thread_id, w->tor_offset, w->p_tor);

              /* calculate shell Green's functions */
              poltor_calc_tor(r, theta,
                              w->green_p[thread_id]->Ynm,
                              w->green_p[thread_id]->dYnm,
                              (complex double *) vx_tor.vector.data,
                              (complex double *) vy_tor.vector.data,
                              (complex double *) vz_tor.vector.data,
                              w);
            }

          if (x != NULL)
            {
              /* compute ionsphere part of B_model = x . dB */
              v = gsl_matrix_complex_row(w->omp_dX, thread_id);
              gsl_blas_zdotu(&v.vector, x, &B_model[0]);

              v = gsl_matrix_complex_row(w->omp_dY, thread_id);
              gsl_blas_zdotu(&v.vector, x, &B_model[1]);

              v = gsl_matrix_complex_row(w->omp_dZ, thread_id);
              gsl_blas_zdotu(&v.vector, x, &B_model[2]);

              /* compute total model: B_model = alpha*B_iono + B_prior */
              for (k = 0; k < 3; ++k)
                {
                  GSL_REAL(B_model[k]) *= alpha;
                  GSL_REAL(B_model[k]) += B_prior[k];
                }

              F_model = gsl_hypot3(GSL_REAL(B_model[0]), GSL_REAL(B_model[1]), GSL_REAL(B_model[2]));

              /* compute unit vector */
              for (k = 0; k < 3; ++k)
                {
                  assert(fabs(GSL_IMAG(B_model[k])) < 1.0e-8); /* sanity check */
                  b_model[k] = gsl_complex_div_real(B_model[k], F_model);
                }
            }
          else
            {
              /* set B_iono = 0, so B_model = B_prior */
              GSL_SET_COMPLEX(&B_model[0], B_prior[0], 0.0);
              GSL_SET_COMPLEX(&B_model[1], B_prior[1], 0.0);
              GSL_SET_COMPLEX(&B_model[2], B_prior[2], 0.0);
              F_model = gsl_hypot3(GSL_REAL(B_model[0]), GSL_REAL(B_model[1]), GSL_REAL(B_model[2]));
            }

          /* calculate internal Green's functions for gradient point (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DF_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW | MAGDATA_FLG_DF_EW))
            {
              /* compute prior model for gradient point (B_main + B_crust + B_ext) */
              magdata_prior_grad(j, B_prior_grad, mptr);

              /* compute Y_{nm} = P_{nm}(cos(theta)) exp(i m phi) and d/dtheta Y_{nm} functions,
               * for gradient point, stored in w->green_grad_p[thread_id]->{Ynm,dYnm} */
              green_complex_Ynm_deriv(mptr->theta_ns[j], mptr->phi_ns[j], w->green_grad_p[thread_id]);

              if (w->p_pint > 0)
                {
                  vx_grad_int = gsl_matrix_complex_subrow(w->omp_dX_grad, thread_id, w->pint_offset, w->p_pint);
                  vy_grad_int = gsl_matrix_complex_subrow(w->omp_dY_grad, thread_id, w->pint_offset, w->p_pint);
                  vz_grad_int = gsl_matrix_complex_subrow(w->omp_dZ_grad, thread_id, w->pint_offset, w->p_pint);

                  /* calculate internal Green's functions */
                  green_complex_calc_int(w->nmax_int, w->mmax_int, w->R, mptr->r_ns[j], mptr->theta_ns[j],
                                         w->green_grad_p[thread_id]->Ynm,
                                         w->green_grad_p[thread_id]->dYnm,
                                         (complex double *) vx_grad_int.vector.data,
                                         (complex double *) vy_grad_int.vector.data,
                                         (complex double *) vz_grad_int.vector.data);
                }

              if (w->p_pext > 0)
                {
                  vx_grad_ext = gsl_matrix_complex_subrow(w->omp_dX_grad, thread_id, w->pext_offset, w->p_pext);
                  vy_grad_ext = gsl_matrix_complex_subrow(w->omp_dY_grad, thread_id, w->pext_offset, w->p_pext);
                  vz_grad_ext = gsl_matrix_complex_subrow(w->omp_dZ_grad, thread_id, w->pext_offset, w->p_pext);

                  /* calculate external Green's functions */
                  green_complex_calc_ext(w->nmax_ext, w->mmax_ext, w->R, mptr->r_ns[j], mptr->theta_ns[j],
                                         w->green_grad_p[thread_id]->Ynm,
                                         w->green_grad_p[thread_id]->dYnm,
                                         (complex double *) vx_grad_ext.vector.data,
                                         (complex double *) vy_grad_ext.vector.data,
                                         (complex double *) vz_grad_ext.vector.data);
                }

              if (w->p_psh > 0)
                {
                  vx_grad_sh = gsl_matrix_complex_subrow(w->omp_dX_grad, thread_id, w->psh_offset, w->p_psh);
                  vy_grad_sh = gsl_matrix_complex_subrow(w->omp_dY_grad, thread_id, w->psh_offset, w->p_psh);
                  vz_grad_sh = gsl_matrix_complex_subrow(w->omp_dZ_grad, thread_id, w->psh_offset, w->p_psh);

                  /* calculate shell Green's functions */
                  poltor_calc_psh(mptr->r_ns[j], mptr->theta_ns[j],
                                  w->green_grad_p[thread_id]->Ynm,
                                  w->green_grad_p[thread_id]->dYnm,
                                  (complex double *) vx_grad_sh.vector.data,
                                  (complex double *) vy_grad_sh.vector.data,
                                  (complex double *) vz_grad_sh.vector.data,
                                  w);
                }

              if (w->p_tor > 0)
                {
                  vx_grad_tor = gsl_matrix_complex_subrow(w->omp_dX_grad, thread_id, w->tor_offset, w->p_tor);
                  vy_grad_tor = gsl_matrix_complex_subrow(w->omp_dY_grad, thread_id, w->tor_offset, w->p_tor);
                  vz_grad_tor = gsl_matrix_complex_subrow(w->omp_dZ_grad, thread_id, w->tor_offset, w->p_tor);

                  /* calculate shell Green's functions */
                  poltor_calc_tor(mptr->r_ns[j], mptr->theta_ns[j],
                                  w->green_grad_p[thread_id]->Ynm,
                                  w->green_grad_p[thread_id]->dYnm,
                                  (complex double *) vx_grad_tor.vector.data,
                                  (complex double *) vy_grad_tor.vector.data,
                                  (complex double *) vz_grad_tor.vector.data,
                                  w);
                }

              if (x != NULL)
                {
                  /* compute ionosphere part of B_model_grad = x . dB_grad */
                  v = gsl_matrix_complex_row(w->omp_dX_grad, thread_id);
                  gsl_blas_zdotu(&v.vector, x, &B_model_grad[0]);

                  v = gsl_matrix_complex_row(w->omp_dY_grad, thread_id);
                  gsl_blas_zdotu(&v.vector, x, &B_model_grad[1]);

                  v = gsl_matrix_complex_row(w->omp_dZ_grad, thread_id);
                  gsl_blas_zdotu(&v.vector, x, &B_model_grad[2]);

                  /* compute total model: B_model = alpha*B_iono + B_prior */
                  for (k = 0; k < 3; ++k)
                    {
                      GSL_REAL(B_model_grad[k]) *= alpha_grad;
                      GSL_REAL(B_model_grad[k]) += B_prior_grad[k];
                    }

                  F_model_grad = gsl_hypot3(GSL_REAL(B_model_grad[0]), GSL_REAL(B_model_grad[1]), GSL_REAL(B_model_grad[2]));

                  /* compute unit vector */
                  for (k = 0; k < 3; ++k)
                    {
                      assert(fabs(GSL_IMAG(B_model_grad[k])) < 1.0e-8); /* sanity check */
                      b_model_grad[k] = gsl_complex_div_real(B_model_grad[k], F_model_grad);
                    }
                }
              else
                {
                  /* set B_iono = 0, so B_model_grad = B_prior_grad */
                  GSL_SET_COMPLEX(&B_model_grad[0], B_prior_grad[0], 0.0);
                  GSL_SET_COMPLEX(&B_model_grad[1], B_prior_grad[1], 0.0);
                  GSL_SET_COMPLEX(&B_model_grad[2], B_prior_grad[2], 0.0);
                  F_model_grad = gsl_hypot3(GSL_REAL(B_model_grad[0]), GSL_REAL(B_model_grad[1]), GSL_REAL(B_model_grad[2]));
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (w->p_pint > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pint_offset, w->p_pint);
                      gsl_vector_complex_memcpy(&v.vector, &vx_int.vector);
                    }

                  if (w->p_pext > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pext_offset, w->p_pext);
                      gsl_vector_complex_memcpy(&v.vector, &vx_ext.vector);
                    }

                  if (w->p_psh > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->psh_offset, w->p_psh);
                      gsl_vector_complex_memcpy(&v.vector, &vx_sh.vector);
                    }

                  if (w->p_tor > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->tor_offset, w->p_tor);
                      gsl_vector_complex_memcpy(&v.vector, &vx_tor.vector);
                    }

                  /* compute: alpha * sqrt(w_j) * row */
                  v = gsl_matrix_complex_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]);
                  gsl_blas_zdscal(alpha * sqrt(wj), &v.vector);

                  GSL_SET_REAL(&fval, sqrt(wj) * (GSL_REAL(B_model[0]) - mptr->Bx_nec[j]));
                  gsl_vector_complex_set(w->omp_f[thread_id], w->omp_rowidx[thread_id]++, fval);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (w->p_pint > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pint_offset, w->p_pint);
                      gsl_vector_complex_memcpy(&v.vector, &vy_int.vector);
                    }

                  if (w->p_pext > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pext_offset, w->p_pext);
                      gsl_vector_complex_memcpy(&v.vector, &vy_ext.vector);
                    }

                  if (w->p_psh > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->psh_offset, w->p_psh);
                      gsl_vector_complex_memcpy(&v.vector, &vy_sh.vector);
                    }

                  if (w->p_tor > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->tor_offset, w->p_tor);
                      gsl_vector_complex_memcpy(&v.vector, &vy_tor.vector);
                    }

                  /* compute: alpha * sqrt(w_j) * row */
                  v = gsl_matrix_complex_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]);
                  gsl_blas_zdscal(alpha * sqrt(wj), &v.vector);

                  GSL_SET_REAL(&fval, sqrt(wj) * (GSL_REAL(B_model[1]) - mptr->By_nec[j]));
                  gsl_vector_complex_set(w->omp_f[thread_id], w->omp_rowidx[thread_id]++, fval);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (w->p_pint > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pint_offset, w->p_pint);
                      gsl_vector_complex_memcpy(&v.vector, &vz_int.vector);
                    }

                  if (w->p_pext > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pext_offset, w->p_pext);
                      gsl_vector_complex_memcpy(&v.vector, &vz_ext.vector);
                    }

                  if (w->p_psh > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->psh_offset, w->p_psh);
                      gsl_vector_complex_memcpy(&v.vector, &vz_sh.vector);
                    }

                  if (w->p_tor > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->tor_offset, w->p_tor);
                      gsl_vector_complex_memcpy(&v.vector, &vz_tor.vector);
                    }

                  /* compute: alpha * sqrt(w_j) * row */
                  v = gsl_matrix_complex_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]);
                  gsl_blas_zdscal(alpha * sqrt(wj), &v.vector);

                  GSL_SET_REAL(&fval, sqrt(wj) * (GSL_REAL(B_model[2]) - mptr->Bz_nec[j]));
                  gsl_vector_complex_set(w->omp_f[thread_id], w->omp_rowidx[thread_id]++, fval);
                }
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double sfac = gsl_vector_get(w->solar_flux, ridx);
              double wj = gsl_vector_get(weights, ridx++);
              if (x != NULL)
                {
                  gsl_complex dBnm[3];
                  gsl_vector_complex_view vdBnm = gsl_vector_complex_view_array((double *) dBnm, 3);
                  gsl_complex val;

                  /* get a view of the current row of the Jacobian */
                  v = gsl_matrix_complex_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]);

                  for (k = 0; k < w->p_pint; ++k)
                    {
                      dBnm[0] = gsl_vector_complex_get(&vx_int.vector, k);
                      dBnm[1] = gsl_vector_complex_get(&vy_int.vector, k);
                      dBnm[2] = gsl_vector_complex_get(&vz_int.vector, k);

                      gsl_blas_zdotu(&vdBnm.vector, &vb_model.vector, &val);
                      gsl_vector_complex_set(&v.vector, k + w->pint_offset, val);
                    }

                  for (k = 0; k < w->p_pext; ++k)
                    {
                      dBnm[0] = gsl_vector_complex_get(&vx_ext.vector, k);
                      dBnm[1] = gsl_vector_complex_get(&vy_ext.vector, k);
                      dBnm[2] = gsl_vector_complex_get(&vz_ext.vector, k);

                      gsl_blas_zdotu(&vdBnm.vector, &vb_model.vector, &val);
                      gsl_vector_complex_set(&v.vector, k + w->pext_offset, val);
                    }

                  for (k = 0; k < w->p_psh; ++k)
                    {
                      dBnm[0] = gsl_vector_complex_get(&vx_sh.vector, k);
                      dBnm[1] = gsl_vector_complex_get(&vy_sh.vector, k);
                      dBnm[2] = gsl_vector_complex_get(&vz_sh.vector, k);

                      gsl_blas_zdotu(&vdBnm.vector, &vb_model.vector, &val);
                      gsl_vector_complex_set(&v.vector, k + w->psh_offset, val);
                    }

                  for (k = 0; k < w->p_tor; ++k)
                    {
                      dBnm[0] = gsl_vector_complex_get(&vx_tor.vector, k);
                      dBnm[1] = gsl_vector_complex_get(&vy_tor.vector, k);
                      dBnm[2] = gsl_vector_complex_get(&vz_tor.vector, k);

                      gsl_blas_zdotu(&vdBnm.vector, &vb_model.vector, &val);
                      gsl_vector_complex_set(&v.vector, k + w->tor_offset, val);
                    }

                  /* scale row by: alpha * sqrt(wj) */
                  gsl_blas_zdscal(alpha * sqrt(wj), &v.vector);
                  assert(alpha == sfac);

                  GSL_SET_REAL(&fval, sqrt(wj) * (F_model - mptr->F[j]));
                  gsl_vector_complex_set(w->omp_f[thread_id], w->omp_rowidx[thread_id]++, fval);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              double wj = gsl_vector_get(weights, ridx++);
              double sqrt_wj = sqrt(wj);
              gsl_complex fac = gsl_complex_rect(-sqrt_wj * alpha, 0.0);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (w->p_pint > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pint_offset, w->p_pint);
                      gsl_vector_complex_memcpy(&v.vector, &vx_grad_int.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vx_int.vector, &v.vector);
                    }

                  if (w->p_pext > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pext_offset, w->p_pext);
                      gsl_vector_complex_memcpy(&v.vector, &vx_grad_ext.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vx_ext.vector, &v.vector);
                    }

                  if (w->p_psh > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->psh_offset, w->p_psh);
                      gsl_vector_complex_memcpy(&v.vector, &vx_grad_sh.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vx_sh.vector, &v.vector);
                    }

                  if (w->p_tor > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->tor_offset, w->p_tor);
                      gsl_vector_complex_memcpy(&v.vector, &vx_grad_tor.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vx_tor.vector, &v.vector);
                    }

                  GSL_SET_REAL(&fval, sqrt_wj * (GSL_REAL(B_model_grad[0]) - GSL_REAL(B_model[0]) - (mptr->Bx_nec_ns[j] - mptr->Bx_nec[j])));
                  gsl_vector_complex_set(w->omp_f[thread_id], w->omp_rowidx[thread_id]++, fval);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              double wj = gsl_vector_get(weights, ridx++);
              double sqrt_wj = sqrt(wj);
              gsl_complex fac = gsl_complex_rect(-sqrt_wj * alpha, 0.0);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (w->p_pint > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pint_offset, w->p_pint);
                      gsl_vector_complex_memcpy(&v.vector, &vy_grad_int.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vy_int.vector, &v.vector);
                    }

                  if (w->p_pext > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pext_offset, w->p_pext);
                      gsl_vector_complex_memcpy(&v.vector, &vy_grad_ext.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vy_ext.vector, &v.vector);
                    }

                  if (w->p_psh > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->psh_offset, w->p_psh);
                      gsl_vector_complex_memcpy(&v.vector, &vy_grad_sh.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vy_sh.vector, &v.vector);
                    }

                  if (w->p_tor > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->tor_offset, w->p_tor);
                      gsl_vector_complex_memcpy(&v.vector, &vy_grad_tor.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vy_tor.vector, &v.vector);
                    }

                  GSL_SET_REAL(&fval, sqrt_wj * (GSL_REAL(B_model_grad[1]) - GSL_REAL(B_model[1]) - (mptr->By_nec_ns[j] - mptr->By_nec[j])));
                  gsl_vector_complex_set(w->omp_f[thread_id], w->omp_rowidx[thread_id]++, fval);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              double wj = gsl_vector_get(weights, ridx++);
              double sqrt_wj = sqrt(wj);
              gsl_complex fac = gsl_complex_rect(-sqrt_wj * alpha, 0.0);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (w->p_pint > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pint_offset, w->p_pint);
                      gsl_vector_complex_memcpy(&v.vector, &vz_grad_int.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vz_int.vector, &v.vector);
                    }

                  if (w->p_pext > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->pext_offset, w->p_pext);
                      gsl_vector_complex_memcpy(&v.vector, &vz_grad_ext.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vz_ext.vector, &v.vector);
                    }

                  if (w->p_psh > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->psh_offset, w->p_psh);
                      gsl_vector_complex_memcpy(&v.vector, &vz_grad_sh.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vz_sh.vector, &v.vector);
                    }

                  if (w->p_tor > 0)
                    {
                      v = gsl_matrix_complex_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id], w->tor_offset, w->p_tor);
                      gsl_vector_complex_memcpy(&v.vector, &vz_grad_tor.vector);
                      gsl_blas_zdscal(sqrt_wj * alpha_grad, &v.vector);
                      gsl_blas_zaxpy(fac, &vz_tor.vector, &v.vector);
                    }

                  GSL_SET_REAL(&fval, sqrt_wj * (GSL_REAL(B_model_grad[2]) - GSL_REAL(B_model[2]) - (mptr->Bz_nec_ns[j] - mptr->Bz_nec[j])));
                  gsl_vector_complex_set(w->omp_f[thread_id], w->omp_rowidx[thread_id]++, fval);
                }
            }

          if (MAGDATA_FitMF(mptr->flags[j]) && mptr->flags[j] & (MAGDATA_FLG_DF_NS | MAGDATA_FLG_DF_EW))
            {
              double wj = gsl_vector_get(weights, ridx++);
              double sqrt_wj = sqrt(wj);
              if (x != NULL)
                {
                  gsl_complex dBnm[3], dBnm_grad[3];
                  gsl_vector_complex_view vdBnm = gsl_vector_complex_view_array((double *) dBnm, 3);
                  gsl_vector_complex_view vdBnm_grad = gsl_vector_complex_view_array((double *) dBnm_grad, 3);
                  gsl_complex val, val_grad;

                  /* get a view of the current row of the Jacobian */
                  v = gsl_matrix_complex_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]);

                  for (k = 0; k < w->p_pint; ++k)
                    {
                      dBnm[0] = gsl_vector_complex_get(&vx_int.vector, k);
                      dBnm[1] = gsl_vector_complex_get(&vy_int.vector, k);
                      dBnm[2] = gsl_vector_complex_get(&vz_int.vector, k);

                      dBnm_grad[0] = gsl_vector_complex_get(&vx_grad_int.vector, k);
                      dBnm_grad[1] = gsl_vector_complex_get(&vy_grad_int.vector, k);
                      dBnm_grad[2] = gsl_vector_complex_get(&vz_grad_int.vector, k);

                      gsl_blas_zdotu(&vdBnm.vector, &vb_model.vector, &val);
                      gsl_blas_zdotu(&vdBnm_grad.vector, &vb_model_grad.vector, &val_grad);

                      gsl_complex_mul_real(val, sqrt_wj * alpha);
                      gsl_complex_mul_real(val_grad, sqrt_wj * alpha_grad);

                      gsl_vector_complex_set(&v.vector, k + w->pint_offset, gsl_complex_sub(val_grad, val));
                    }

                  for (k = 0; k < w->p_pext; ++k)
                    {
                      dBnm[0] = gsl_vector_complex_get(&vx_ext.vector, k);
                      dBnm[1] = gsl_vector_complex_get(&vy_ext.vector, k);
                      dBnm[2] = gsl_vector_complex_get(&vz_ext.vector, k);

                      dBnm_grad[0] = gsl_vector_complex_get(&vx_grad_ext.vector, k);
                      dBnm_grad[1] = gsl_vector_complex_get(&vy_grad_ext.vector, k);
                      dBnm_grad[2] = gsl_vector_complex_get(&vz_grad_ext.vector, k);

                      gsl_blas_zdotu(&vdBnm.vector, &vb_model.vector, &val);
                      gsl_blas_zdotu(&vdBnm_grad.vector, &vb_model_grad.vector, &val_grad);

                      gsl_complex_mul_real(val, sqrt_wj * alpha);
                      gsl_complex_mul_real(val_grad, sqrt_wj * alpha_grad);

                      gsl_vector_complex_set(&v.vector, k + w->pext_offset, gsl_complex_sub(val_grad, val));
                    }

                  for (k = 0; k < w->p_psh; ++k)
                    {
                      dBnm[0] = gsl_vector_complex_get(&vx_sh.vector, k);
                      dBnm[1] = gsl_vector_complex_get(&vy_sh.vector, k);
                      dBnm[2] = gsl_vector_complex_get(&vz_sh.vector, k);

                      dBnm_grad[0] = gsl_vector_complex_get(&vx_grad_sh.vector, k);
                      dBnm_grad[1] = gsl_vector_complex_get(&vy_grad_sh.vector, k);
                      dBnm_grad[2] = gsl_vector_complex_get(&vz_grad_sh.vector, k);

                      gsl_blas_zdotu(&vdBnm.vector, &vb_model.vector, &val);
                      gsl_blas_zdotu(&vdBnm_grad.vector, &vb_model_grad.vector, &val_grad);

                      gsl_complex_mul_real(val, sqrt_wj * alpha);
                      gsl_complex_mul_real(val_grad, sqrt_wj * alpha_grad);

                      gsl_vector_complex_set(&v.vector, k + w->psh_offset, gsl_complex_sub(val_grad, val));
                    }

                  for (k = 0; k < w->p_tor; ++k)
                    {
                      dBnm[0] = gsl_vector_complex_get(&vx_tor.vector, k);
                      dBnm[1] = gsl_vector_complex_get(&vy_tor.vector, k);
                      dBnm[2] = gsl_vector_complex_get(&vz_tor.vector, k);

                      dBnm_grad[0] = gsl_vector_complex_get(&vx_grad_tor.vector, k);
                      dBnm_grad[1] = gsl_vector_complex_get(&vy_grad_tor.vector, k);
                      dBnm_grad[2] = gsl_vector_complex_get(&vz_grad_tor.vector, k);

                      gsl_blas_zdotu(&vdBnm.vector, &vb_model.vector, &val);
                      gsl_blas_zdotu(&vdBnm_grad.vector, &vb_model_grad.vector, &val_grad);

                      gsl_complex_mul_real(val, sqrt_wj * alpha);
                      gsl_complex_mul_real(val_grad, sqrt_wj * alpha_grad);

                      gsl_vector_complex_set(&v.vector, k + w->tor_offset, gsl_complex_sub(val_grad, val));
                    }

                  GSL_SET_REAL(&fval, sqrt_wj * (F_model_grad - mptr->F_ns[j] - (F_model - mptr->F[j])));
                  gsl_vector_complex_set(w->omp_f[thread_id], w->omp_rowidx[thread_id]++, fval);
                }
            }

          /*
           * check if omp_J[thread_id] is full and should be folded into JHJ; the
           * 8 is just some slop to prevent trying to fill rows past the matrix buffer
           * in the loop above
           */
          if (w->omp_rowidx[thread_id] >= w->nblock - 8)
            {
              /* fold current matrix block into JHJ, one thread at a time */
              gsl_matrix_complex_view m = gsl_matrix_complex_submatrix(w->omp_J[thread_id], 0, 0, w->omp_rowidx[thread_id], w->p);
              gsl_vector_complex_view v = gsl_vector_complex_subvector(w->omp_f[thread_id], 0, w->omp_rowidx[thread_id]);

              /* keep cumulative total of rows processed by this thread for progress bar */
              w->omp_nrows[thread_id] += w->omp_rowidx[thread_id];
              w->omp_rowidx[thread_id] = 0;

#pragma omp critical
              {
                /* JHJ += m^H m */
                gsl_blas_zherk(CblasUpper, CblasConjTrans, 1.0, &m.matrix, 1.0, w->JHJ);

                /* JHf += m^H v */
                gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, &m.matrix, &v.vector, GSL_COMPLEX_ONE, w->JHf);
              }

              if (thread_id == 0)
                {
                  double progress = 0.0;
                  size_t k;

                  for (k = 0; k < w->max_threads; ++k)
                    progress += (double) w->omp_nrows[k];
    
                  progress /= (double) w->n;

                  fprintf(stderr, "\t");
                  progress_bar(stderr, progress, 70);
                }
            }
        }
    } /* for (i = 0; i < list->n; ++i) */

  /* now loop through to see if any rows were not accumulated into JTJ_vec */
  for (i = 0; i < w->max_threads; ++i)
    {
      if (w->omp_rowidx[i] > 0)
        {
          /* accumulate final Green's functions into JHJ */
          gsl_matrix_complex_view m = gsl_matrix_complex_submatrix(w->omp_J[i], 0, 0, w->omp_rowidx[i], w->p);
          gsl_vector_complex_view v = gsl_vector_complex_subvector(w->omp_f[i], 0, w->omp_rowidx[i]);

          gsl_blas_zherk(CblasUpper, CblasConjTrans, 1.0, &m.matrix, 1.0, w->JHJ);
          gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, &m.matrix, &v.vector, GSL_COMPLEX_ONE, w->JHf);
        }
    }

  fprintf(stderr, "\t");
  progress_bar(stderr, 1.0, 70);

#if 0
  printherm_octave(w->JHJ, "JHJ");
  printcv_octave(w->JHf, "JHf");
  exit(1);
#endif

  return s;
}

/*
poltor_calc_f()
  Calculate residual vector for a given set of coefficients

Inputs: x      - model coefficients
        params - workspace
        f      - (output) residual vector, size n
*/

static int
poltor_calc_f(gsl_vector_complex * x, void * params, gsl_vector * f)
{
  poltor_workspace * w = (poltor_workspace *) params;
  magdata_list *list = w->data;
  const gsl_vector *weights = w->wts_final;
  size_t i;

  gsl_vector_set_zero(f);

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);
      size_t j;

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index of this data point in [0,w->n-1] */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          gsl_complex fval;
          double B_prior[4], B_prior_grad[4]; /* a priori field models (main, crust, magnetosphere) */
          gsl_complex B_model[3];       /* B_model (ionosphere + prior) */
          gsl_complex B_model_grad[3];  /* B_model (ionosphere + prior) for gradient point */
          double F_model, F_model_grad; /* || B_model ||, || B_model_grad || */
          double alpha = gsl_vector_get(w->solar_flux, ridx); /* alpha factor = 1 + N*EUVAC (same for all residuals for this data point) */
          size_t k;

          GSL_SET_IMAG(&fval, 0.0);

          gsl_vector_complex_view vx_int, vy_int, vz_int;
          gsl_vector_complex_view vx_ext, vy_ext, vz_ext;
          gsl_vector_complex_view vx_sh, vy_sh, vz_sh;
          gsl_vector_complex_view vx_tor, vy_tor, vz_tor;
          gsl_vector_complex_view vx_grad_int, vy_grad_int, vz_grad_int;
          gsl_vector_complex_view vx_grad_ext, vy_grad_ext, vz_grad_ext;
          gsl_vector_complex_view vx_grad_sh, vy_grad_sh, vz_grad_sh;
          gsl_vector_complex_view vx_grad_tor, vy_grad_tor, vz_grad_tor;
          gsl_vector_complex_view v;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute prior model (B_main + B_crust + B_ext) */
          magdata_prior(j, B_prior, mptr);

          /* compute Y_{nm} = P_{nm}(cos(theta)) exp(i m phi) and d/dtheta Y_{nm} functions,
           * stored in w->green_p[thread_id]->{Ynm,dYnm} */
          green_complex_Ynm_deriv(theta, phi, w->green_p[thread_id]);

          if (w->p_pint > 0)
            {
              vx_int = gsl_matrix_complex_subrow(w->omp_dX, thread_id, w->pint_offset, w->p_pint);
              vy_int = gsl_matrix_complex_subrow(w->omp_dY, thread_id, w->pint_offset, w->p_pint);
              vz_int = gsl_matrix_complex_subrow(w->omp_dZ, thread_id, w->pint_offset, w->p_pint);

              /* calculate internal Green's functions */
              green_complex_calc_int(w->nmax_int, w->mmax_int, w->R, r, theta,
                                     w->green_p[thread_id]->Ynm,
                                     w->green_p[thread_id]->dYnm,
                                     (complex double *) vx_int.vector.data,
                                     (complex double *) vy_int.vector.data,
                                     (complex double *) vz_int.vector.data);
            }

          if (w->p_pext > 0)
            {
              vx_ext = gsl_matrix_complex_subrow(w->omp_dX, thread_id, w->pext_offset, w->p_pext);
              vy_ext = gsl_matrix_complex_subrow(w->omp_dY, thread_id, w->pext_offset, w->p_pext);
              vz_ext = gsl_matrix_complex_subrow(w->omp_dZ, thread_id, w->pext_offset, w->p_pext);

              /* calculate external Green's functions */
              green_complex_calc_ext(w->nmax_ext, w->mmax_ext, w->R, r, theta,
                                     w->green_p[thread_id]->Ynm,
                                     w->green_p[thread_id]->dYnm,
                                     (complex double *) vx_ext.vector.data,
                                     (complex double *) vy_ext.vector.data,
                                     (complex double *) vz_ext.vector.data);
            }

          if (w->p_psh > 0)
            {
              vx_sh = gsl_matrix_complex_subrow(w->omp_dX, thread_id, w->psh_offset, w->p_psh);
              vy_sh = gsl_matrix_complex_subrow(w->omp_dY, thread_id, w->psh_offset, w->p_psh);
              vz_sh = gsl_matrix_complex_subrow(w->omp_dZ, thread_id, w->psh_offset, w->p_psh);

              /* calculate shell Green's functions */
              poltor_calc_psh(r, theta,
                              w->green_p[thread_id]->Ynm,
                              w->green_p[thread_id]->dYnm,
                              (complex double *) vx_sh.vector.data,
                              (complex double *) vy_sh.vector.data,
                              (complex double *) vz_sh.vector.data,
                              w);
            }

          if (w->p_tor > 0)
            {
              vx_tor = gsl_matrix_complex_subrow(w->omp_dX, thread_id, w->tor_offset, w->p_tor);
              vy_tor = gsl_matrix_complex_subrow(w->omp_dY, thread_id, w->tor_offset, w->p_tor);
              vz_tor = gsl_matrix_complex_subrow(w->omp_dZ, thread_id, w->tor_offset, w->p_tor);

              /* calculate shell Green's functions */
              poltor_calc_tor(r, theta,
                              w->green_p[thread_id]->Ynm,
                              w->green_p[thread_id]->dYnm,
                              (complex double *) vx_tor.vector.data,
                              (complex double *) vy_tor.vector.data,
                              (complex double *) vz_tor.vector.data,
                              w);
            }

          /* compute ionsphere part of B_model = x . dB */
          v = gsl_matrix_complex_row(w->omp_dX, thread_id);
          gsl_blas_zdotu(&v.vector, x, &B_model[0]);

          v = gsl_matrix_complex_row(w->omp_dY, thread_id);
          gsl_blas_zdotu(&v.vector, x, &B_model[1]);

          v = gsl_matrix_complex_row(w->omp_dZ, thread_id);
          gsl_blas_zdotu(&v.vector, x, &B_model[2]);

          /* compute total model: B_model = alpha*B_iono + B_prior */
          for (k = 0; k < 3; ++k)
            {
              GSL_REAL(B_model[k]) *= alpha;
              GSL_REAL(B_model[k]) += B_prior[k];
            }

          F_model = gsl_hypot3(GSL_REAL(B_model[0]), GSL_REAL(B_model[1]), GSL_REAL(B_model[2]));

          /* calculate internal Green's functions for gradient point (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DF_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW | MAGDATA_FLG_DF_EW))
            {
              double alpha_grad = gsl_vector_get(w->solar_flux_grad, ridx);

              /* compute prior model for gradient point (B_main + B_crust + B_ext) */
              magdata_prior_grad(j, B_prior_grad, mptr);

              /* compute Y_{nm} = P_{nm}(cos(theta)) exp(i m phi) and d/dtheta Y_{nm} functions,
               * for gradient point, stored in w->green_grad_p[thread_id]->{Ynm,dYnm} */
              green_complex_Ynm_deriv(mptr->theta_ns[j], mptr->phi_ns[j], w->green_grad_p[thread_id]);

              if (w->p_pint > 0)
                {
                  vx_grad_int = gsl_matrix_complex_subrow(w->omp_dX_grad, thread_id, w->pint_offset, w->p_pint);
                  vy_grad_int = gsl_matrix_complex_subrow(w->omp_dY_grad, thread_id, w->pint_offset, w->p_pint);
                  vz_grad_int = gsl_matrix_complex_subrow(w->omp_dZ_grad, thread_id, w->pint_offset, w->p_pint);

                  /* calculate internal Green's functions */
                  green_complex_calc_int(w->nmax_int, w->mmax_int, w->R, mptr->r_ns[j], mptr->theta_ns[j],
                                         w->green_grad_p[thread_id]->Ynm,
                                         w->green_grad_p[thread_id]->dYnm,
                                         (complex double *) vx_grad_int.vector.data,
                                         (complex double *) vy_grad_int.vector.data,
                                         (complex double *) vz_grad_int.vector.data);
                }

              if (w->p_pext > 0)
                {
                  vx_grad_ext = gsl_matrix_complex_subrow(w->omp_dX_grad, thread_id, w->pext_offset, w->p_pext);
                  vy_grad_ext = gsl_matrix_complex_subrow(w->omp_dY_grad, thread_id, w->pext_offset, w->p_pext);
                  vz_grad_ext = gsl_matrix_complex_subrow(w->omp_dZ_grad, thread_id, w->pext_offset, w->p_pext);

                  /* calculate external Green's functions */
                  green_complex_calc_ext(w->nmax_ext, w->mmax_ext, w->R, mptr->r_ns[j], mptr->theta_ns[j],
                                         w->green_grad_p[thread_id]->Ynm,
                                         w->green_grad_p[thread_id]->dYnm,
                                         (complex double *) vx_grad_ext.vector.data,
                                         (complex double *) vy_grad_ext.vector.data,
                                         (complex double *) vz_grad_ext.vector.data);
                }

              if (w->p_psh > 0)
                {
                  vx_grad_sh = gsl_matrix_complex_subrow(w->omp_dX_grad, thread_id, w->psh_offset, w->p_psh);
                  vy_grad_sh = gsl_matrix_complex_subrow(w->omp_dY_grad, thread_id, w->psh_offset, w->p_psh);
                  vz_grad_sh = gsl_matrix_complex_subrow(w->omp_dZ_grad, thread_id, w->psh_offset, w->p_psh);

                  /* calculate shell Green's functions */
                  poltor_calc_psh(mptr->r_ns[j], mptr->theta_ns[j],
                                  w->green_grad_p[thread_id]->Ynm,
                                  w->green_grad_p[thread_id]->dYnm,
                                  (complex double *) vx_grad_sh.vector.data,
                                  (complex double *) vy_grad_sh.vector.data,
                                  (complex double *) vz_grad_sh.vector.data,
                                  w);
                }

              if (w->p_tor > 0)
                {
                  vx_grad_tor = gsl_matrix_complex_subrow(w->omp_dX_grad, thread_id, w->tor_offset, w->p_tor);
                  vy_grad_tor = gsl_matrix_complex_subrow(w->omp_dY_grad, thread_id, w->tor_offset, w->p_tor);
                  vz_grad_tor = gsl_matrix_complex_subrow(w->omp_dZ_grad, thread_id, w->tor_offset, w->p_tor);

                  /* calculate shell Green's functions */
                  poltor_calc_tor(mptr->r_ns[j], mptr->theta_ns[j],
                                  w->green_grad_p[thread_id]->Ynm,
                                  w->green_grad_p[thread_id]->dYnm,
                                  (complex double *) vx_grad_tor.vector.data,
                                  (complex double *) vy_grad_tor.vector.data,
                                  (complex double *) vz_grad_tor.vector.data,
                                  w);
                }

              /* compute ionosphere part of B_model_grad = x . dB_grad */
              v = gsl_matrix_complex_row(w->omp_dX_grad, thread_id);
              gsl_blas_zdotu(&v.vector, x, &B_model_grad[0]);

              v = gsl_matrix_complex_row(w->omp_dY_grad, thread_id);
              gsl_blas_zdotu(&v.vector, x, &B_model_grad[1]);

              v = gsl_matrix_complex_row(w->omp_dZ_grad, thread_id);
              gsl_blas_zdotu(&v.vector, x, &B_model_grad[2]);

              /* compute total model: B_model = alpha*B_iono + B_prior */
              for (k = 0; k < 3; ++k)
                {
                  GSL_REAL(B_model_grad[k]) *= alpha_grad;
                  GSL_REAL(B_model_grad[k]) += B_prior_grad[k];
                }

              F_model_grad = gsl_hypot3(GSL_REAL(B_model_grad[0]), GSL_REAL(B_model_grad[1]), GSL_REAL(B_model_grad[2]));
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double wj = gsl_vector_get(weights, ridx);
              gsl_vector_set(f, ridx++, sqrt(wj) * (GSL_REAL(B_model[0]) - mptr->Bx_nec[j]));
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double wj = gsl_vector_get(weights, ridx);
              gsl_vector_set(f, ridx++, sqrt(wj) * (GSL_REAL(B_model[1]) - mptr->By_nec[j]));
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double wj = gsl_vector_get(weights, ridx);
              gsl_vector_set(f, ridx++, sqrt(wj) * (GSL_REAL(B_model[2]) - mptr->Bz_nec[j]));
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double wj = gsl_vector_get(weights, ridx);
              gsl_vector_set(f, ridx++, sqrt(wj) * (F_model - mptr->F[j]));
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              double wj = gsl_vector_get(weights, ridx);
              gsl_vector_set(f, ridx++, sqrt(wj) * (GSL_REAL(B_model_grad[0]) - GSL_REAL(B_model[0]) - (mptr->Bx_nec_ns[j] - mptr->Bx_nec[j])));
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              double wj = gsl_vector_get(weights, ridx);
              gsl_vector_set(f, ridx++, sqrt(wj) * (GSL_REAL(B_model_grad[1]) - GSL_REAL(B_model[1]) - (mptr->By_nec_ns[j] - mptr->By_nec[j])));
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              double wj = gsl_vector_get(weights, ridx);
              gsl_vector_set(f, ridx++, sqrt(wj) * (GSL_REAL(B_model_grad[2]) - GSL_REAL(B_model[2]) - (mptr->Bz_nec_ns[j] - mptr->Bz_nec[j])));
            }

          if (MAGDATA_FitMF(mptr->flags[j]) && mptr->flags[j] & (MAGDATA_FLG_DF_NS | MAGDATA_FLG_DF_EW))
            {
              double wj = gsl_vector_get(weights, ridx);
              gsl_vector_set(f, ridx++, sqrt(wj) * (F_model_grad - mptr->F_ns[j] - (F_model - mptr->F[j])));
            }
        }
    }

  return GSL_SUCCESS;
}

static int
poltor_nonlinear_regularize(poltor_workspace *w)
{
  int s = 0;
  const poltor_parameters *params = &(w->params);
  gsl_vector_complex_view diag = gsl_matrix_complex_diagonal(w->JHJ);
  size_t n;

  if (w->p_pint > 0)
    {
      /* regularize internal (Sq) field by minimizing <Br^2> at Earth surface */

      const double alpha_sq = params->alpha_int * params->alpha_int;

      for (n = 1; n <= w->nmax_int; ++n)
        {
          double fac = (n + 1.0) * (n + 1.0) / (2.0 * n + 1.0);
          int M = (int) GSL_MIN(n, w->mmax_int);
          int m;

          for (m = -M; m <= M; ++m)
            {
              size_t cidx = green_idx(n, m, w->mmax_int) + w->pint_offset;
              gsl_complex *dnm = gsl_vector_complex_ptr(&diag.vector, cidx);

              GSL_REAL(*dnm) += alpha_sq * fac;
            }
        }
    }

  if (w->p_psh > 0)
    {
      const double alpha_sq = params->alpha_sh * params->alpha_sh;

#if 0
      /* simple Tikhonov regularization with L = I */
      gsl_complex val = gsl_complex_rect(alpha_sq, 0.0);
      gsl_vector_complex_view v = gsl_vector_complex_subvector(&diag.vector, w->psh_offset, w->p_psh);
      gsl_vector_complex_add_constant(&v.vector, val);

#else

      size_t j;

      for (j = 0; j <= w->shell_J; ++j)
        {
          for (n = 1; n <= w->nmax_sh; ++n)
            {
              double fac = (n + 1.0) * (n + 1.0) / (2.0 * n + 1.0);
              int M = (int) GSL_MIN(n, w->mmax_sh);
              int m;

              for (m = -M; m <= M; ++m)
                {
                  size_t cidx = poltor_jnmidx(j, n, m, w);
                  gsl_complex *dnm = gsl_vector_complex_ptr(&diag.vector, cidx);

                  GSL_REAL(*dnm) += alpha_sq * fac;
                }
            }
        }
#endif
    }

  return s;
}

static int
poltor_robust_print_stat(const char *str, const double sigma, const gsl_rstat_workspace *rstat_p)
{
  const size_t n = gsl_rstat_n(rstat_p);

  if (n > 0)
    {
      const double mean = gsl_rstat_mean(rstat_p);
      fprintf(stderr, "\t %18s = %.2f [nT], Robust weight mean = %.4f\n", str, sigma, mean);
    }

  return 0;
}

/*
poltor_robust_weights()
  Compute robust weights

Inputs: f       - residual vector, size n
        weights - (output) weight vector, size n
        w       - workspace

*/

static int
poltor_robust_weights(const gsl_vector * f, gsl_vector * weights, poltor_workspace * w)
{
  int s = 0;
  const double tune = 1.0;
  const double qdlat_cutoff = 55.0;
  size_t i, j;
  size_t idx = 0;
  magdata_list *list = w->data;
  gsl_rstat_workspace **rstat_x = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_y = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_z = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_f = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_dx_ns = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_dy_ns = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_low_dz_ns = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_high_dz_ns = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_df_ns = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_dx_ew = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_dy_ew = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_low_dz_ew = malloc(list->n * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_high_dz_ew = malloc(list->n * sizeof(gsl_rstat_workspace *));

  for (i = 0; i < list->n; ++i)
    {
      rstat_x[i] = gsl_rstat_alloc();
      rstat_y[i] = gsl_rstat_alloc();
      rstat_z[i] = gsl_rstat_alloc();
      rstat_f[i] = gsl_rstat_alloc();
      rstat_dx_ns[i] = gsl_rstat_alloc();
      rstat_dy_ns[i] = gsl_rstat_alloc();
      rstat_low_dz_ns[i] = gsl_rstat_alloc();
      rstat_high_dz_ns[i] = gsl_rstat_alloc();
      rstat_df_ns[i] = gsl_rstat_alloc();
      rstat_dx_ew[i] = gsl_rstat_alloc();
      rstat_dy_ew[i] = gsl_rstat_alloc();
      rstat_low_dz_ew[i] = gsl_rstat_alloc();
      rstat_high_dz_ew[i] = gsl_rstat_alloc();
    }

  fprintf(stderr, "\n");

  /*
   * first loop through the residuals and compute statistics for each residual type
   * (X,Y,Z,F,DX,DY,DZ)
   */
  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      for (j = 0; j < mptr->n; ++j)
        {
          /* check if data point is discarded due to time interval */
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_ExistX(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                gsl_rstat_add(fi, rstat_x[i]);
            }

          if (MAGDATA_ExistY(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                gsl_rstat_add(fi, rstat_y[i]);
            }

          if (MAGDATA_ExistZ(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                gsl_rstat_add(fi, rstat_z[i]);
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);
              gsl_rstat_add(fi, rstat_f[i]);
            }

          if (MAGDATA_ExistDX_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                gsl_rstat_add(fi, rstat_dx_ns[i]);
            }

          if (MAGDATA_ExistDY_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                gsl_rstat_add(fi, rstat_dy_ns[i]);
            }

          if (MAGDATA_ExistDZ_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                    gsl_rstat_add(fi, rstat_low_dz_ns[i]);
                  else
                    gsl_rstat_add(fi, rstat_high_dz_ns[i]);
                }
            }

          if (MAGDATA_ExistDF_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                gsl_rstat_add(fi, rstat_df_ns[i]);
            }

          if (MAGDATA_ExistDX_EW(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                gsl_rstat_add(fi, rstat_dx_ew[i]);
            }

          if (MAGDATA_ExistDY_EW(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                gsl_rstat_add(fi, rstat_dy_ew[i]);
            }

          if (MAGDATA_ExistDZ_EW(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                    gsl_rstat_add(fi, rstat_low_dz_ew[i]);
                  else
                    gsl_rstat_add(fi, rstat_high_dz_ew[i]);
                }
            }
        }
    }

  assert(idx == w->n);

  /* loop through again and compute robust weights and the mean of the weights */

  idx = 0;
  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);
      const double alpha = 1.0; /* constant to multiply sigma so that mean(weights) = 0.95 */
      double sigma_X = alpha * gsl_rstat_sd(rstat_x[i]);
      double sigma_Y = alpha * gsl_rstat_sd(rstat_y[i]);
      double sigma_Z = alpha * gsl_rstat_sd(rstat_z[i]);
      double sigma_F = alpha * gsl_rstat_sd(rstat_f[i]);
      double sigma_DX_NS = alpha * gsl_rstat_sd(rstat_dx_ns[i]);
      double sigma_DY_NS = alpha * gsl_rstat_sd(rstat_dy_ns[i]);
      double sigma_low_DZ_NS = alpha * gsl_rstat_sd(rstat_low_dz_ns[i]);
      double sigma_high_DZ_NS = alpha * gsl_rstat_sd(rstat_high_dz_ns[i]);
      double sigma_DF_NS = alpha * gsl_rstat_sd(rstat_df_ns[i]);
      double sigma_DX_EW = alpha * gsl_rstat_sd(rstat_dx_ew[i]);
      double sigma_DY_EW = alpha * gsl_rstat_sd(rstat_dy_ew[i]);
      double sigma_low_DZ_EW = alpha * gsl_rstat_sd(rstat_low_dz_ew[i]);
      double sigma_high_DZ_EW = alpha * gsl_rstat_sd(rstat_high_dz_ew[i]);

      gsl_rstat_reset(rstat_x[i]);
      gsl_rstat_reset(rstat_y[i]);
      gsl_rstat_reset(rstat_z[i]);
      gsl_rstat_reset(rstat_f[i]);
      gsl_rstat_reset(rstat_dx_ns[i]);
      gsl_rstat_reset(rstat_dy_ns[i]);
      gsl_rstat_reset(rstat_low_dz_ns[i]);
      gsl_rstat_reset(rstat_high_dz_ns[i]);
      gsl_rstat_reset(rstat_df_ns[i]);
      gsl_rstat_reset(rstat_dx_ew[i]);
      gsl_rstat_reset(rstat_dy_ew[i]);
      gsl_rstat_reset(rstat_low_dz_ew[i]);
      gsl_rstat_reset(rstat_high_dz_ew[i]);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_ExistX(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = huber(fi / (tune * sigma_X));
              gsl_vector_set(weights, idx++, wi);
              gsl_rstat_add(wi, rstat_x[i]);
            }

          if (MAGDATA_ExistY(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = huber(fi / (tune * sigma_Y));
              gsl_vector_set(weights, idx++, wi);
              gsl_rstat_add(wi, rstat_y[i]);
            }

          if (MAGDATA_ExistZ(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = huber(fi / (tune * sigma_Z));
              gsl_vector_set(weights, idx++, wi);
              gsl_rstat_add(wi, rstat_z[i]);
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = huber(fi / (tune * sigma_F));
              gsl_vector_set(weights, idx++, wi);
              gsl_rstat_add(wi, rstat_f[i]);
            }

          if (MAGDATA_ExistDX_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = huber(fi / (tune * sigma_DX_NS));
              gsl_vector_set(weights, idx++, wi);
              gsl_rstat_add(wi, rstat_dx_ns[i]);
            }

          if (MAGDATA_ExistDY_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = huber(fi / (tune * sigma_DY_NS));
              gsl_vector_set(weights, idx++, wi);
              gsl_rstat_add(wi, rstat_dy_ns[i]);
            }

          if (MAGDATA_ExistDZ_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi;
              
              if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                {
                  wi = huber(fi / (tune * sigma_low_DZ_NS));
                  gsl_rstat_add(wi, rstat_low_dz_ns[i]);
                }
              else
                {
                  wi = huber(fi / (tune * sigma_high_DZ_NS));
                  gsl_rstat_add(wi, rstat_high_dz_ns[i]);
                }

              gsl_vector_set(weights, idx++, wi);
            }

          if (MAGDATA_ExistDF_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = huber(fi / (tune * sigma_DF_NS));
              gsl_vector_set(weights, idx++, wi);
              gsl_rstat_add(wi, rstat_df_ns[i]);
            }

          if (MAGDATA_ExistDX_EW(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = huber(fi / (tune * sigma_DX_EW));
              gsl_vector_set(weights, idx++, wi);
              gsl_rstat_add(wi, rstat_dx_ew[i]);
            }

          if (MAGDATA_ExistDY_EW(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = huber(fi / (tune * sigma_DY_EW));
              gsl_vector_set(weights, idx++, wi);
              gsl_rstat_add(wi, rstat_dy_ew[i]);
            }

          if (MAGDATA_ExistDZ_EW(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi;

              if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                {
                  wi = huber(fi / (tune * sigma_low_DZ_EW));
                  gsl_rstat_add(wi, rstat_low_dz_ew[i]);
                }
              else
                {
                  wi = huber(fi / (tune * sigma_high_DZ_EW));
                  gsl_rstat_add(wi, rstat_high_dz_ew[i]);
                }

              gsl_vector_set(weights, idx++, wi);
            }
        }

      fprintf(stderr, "\t === SATELLITE %zu (robust sigma) ===\n", i);

      poltor_robust_print_stat("sigma X", sigma_X, rstat_x[i]);
      poltor_robust_print_stat("sigma Y", sigma_Y, rstat_y[i]);
      poltor_robust_print_stat("sigma Z", sigma_Z, rstat_z[i]);
      poltor_robust_print_stat("sigma F", sigma_F, rstat_f[i]);

      poltor_robust_print_stat("sigma DX_NS", sigma_DX_NS, rstat_dx_ns[i]);
      poltor_robust_print_stat("sigma DY_NS", sigma_DY_NS, rstat_dy_ns[i]);
      poltor_robust_print_stat("sigma low DZ_NS", sigma_low_DZ_NS, rstat_low_dz_ns[i]);
      poltor_robust_print_stat("sigma high DZ_NS", sigma_high_DZ_NS, rstat_high_dz_ns[i]);
      poltor_robust_print_stat("sigma DF_NS", sigma_DF_NS, rstat_df_ns[i]);

      poltor_robust_print_stat("sigma DX_EW", sigma_DX_EW, rstat_dx_ew[i]);
      poltor_robust_print_stat("sigma DY_EW", sigma_DY_EW, rstat_dy_ew[i]);
      poltor_robust_print_stat("sigma low DZ_EW", sigma_low_DZ_EW, rstat_low_dz_ew[i]);
      poltor_robust_print_stat("sigma high DZ_EW", sigma_high_DZ_EW, rstat_high_dz_ew[i]);
    }

  assert(idx == w->n);

  for (i = 0; i < list->n; ++i)
    {
      gsl_rstat_free(rstat_x[i]);
      gsl_rstat_free(rstat_y[i]);
      gsl_rstat_free(rstat_z[i]);
      gsl_rstat_free(rstat_f[i]);
      gsl_rstat_free(rstat_dx_ns[i]);
      gsl_rstat_free(rstat_dy_ns[i]);
      gsl_rstat_free(rstat_low_dz_ns[i]);
      gsl_rstat_free(rstat_high_dz_ns[i]);
      gsl_rstat_free(rstat_df_ns[i]);
      gsl_rstat_free(rstat_dx_ew[i]);
      gsl_rstat_free(rstat_dy_ew[i]);
      gsl_rstat_free(rstat_low_dz_ew[i]);
      gsl_rstat_free(rstat_high_dz_ew[i]);
    }

  free(rstat_x);
  free(rstat_y);
  free(rstat_z);
  free(rstat_f);
  free(rstat_dx_ns);
  free(rstat_dy_ns);
  free(rstat_low_dz_ns);
  free(rstat_high_dz_ns);
  free(rstat_dx_ew);
  free(rstat_dy_ew);
  free(rstat_low_dz_ew);
  free(rstat_high_dz_ew);

  return s;
}

static double
huber(const double x)
{
  const double ax = fabs(x);

  if (ax <= 1.0)
    return 1.0;
  else
    return (1.0 / ax);
}

static double
bisquare(const double x)
{
  if (fabs(x) <= 1.0)
    {
      double f = 1.0 - x*x;
      return (f * f);
    }
  else
    return 0.0;
}
