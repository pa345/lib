/*
 * mfield_multilarge.c
 *
 * Contains routines for fitting magnetic field module using gsl_multilarge_nlinear framework
 */

static int mfield_nonlinear_vector_precompute(const gsl_vector *sqrt_weights, mfield_workspace *w);
static int mfield_vector_green(const gsl_vector * N, const size_t istart, const double sqrt_weight,
                               const gsl_vector *dB, gsl_vector *J, mfield_workspace *w);
static int mfield_calc_df3(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                           void *params, gsl_vector * v, gsl_matrix * JTJ);
static inline int mfield_jacobian_J1Tu(const gsl_vector * N, const size_t istart, const double sqrt_wj, const double uj,
                                       const gsl_vector * dB_int, gsl_vector *JTu, const mfield_workspace *w);
static inline int mfield_jacobian_J1u(const double t, const double sqrt_wt, const gsl_vector * u,
                                      const size_t ridx, const gsl_vector * dB_int, gsl_vector *Ju,
                                      const mfield_workspace *w);
static inline int mfield_jacobian_F(CBLAS_TRANSPOSE_t TransJ, const size_t istart, const gsl_vector * u,
                                    const size_t ridx, const gsl_spmatrix * J2, const gsl_vector *J_int,
                                    gsl_vector *v, gsl_matrix * J2TJ1, const mfield_workspace *w);
static int mfield_jacobian_J2TJ1(const gsl_vector * N, const size_t istart, const double sqrt_wj, const size_t ridx,
                                 const gsl_vector * dB_int, const gsl_spmatrix * J2, gsl_matrix * J2TJ1,
                                 const mfield_workspace * w);

/*
mfield_nonlinear_vector_precompute()
  Precompute J_int^T W J_int for vector measurements, since
this submatrix is independent of the model parameters and
only needs to be computed once per iteration. This function
uses OpenMP to speed up the calculation

Inputs: sqrt_weights - sqrt weight vector
        w            - workspace

Notes:
1) w->JTJ_vec is updated with J_int^T W J_int for vector
residuals

2) J_int^T W J_int is stored in lower triangle of w->JTJ_vec
*/

static int
mfield_nonlinear_vector_precompute(const gsl_vector *sqrt_weights, mfield_workspace *w)
{
  int s = GSL_SUCCESS;
  size_t i, j;
  size_t nres_vec = 0;
  size_t nres_completed = 0;

  gsl_matrix_set_zero(w->JTJ_vec);

  /*
   * w->nres_vec includes vector residuals which don't have a core/crustal
   * field component (i.e. Euler angles). So recalculate nres_vec restricting
   * count to main field fitting
   */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (!MAGDATA_FitMF(mptr->flags[j]))
            continue;

          if (MAGDATA_ExistX(mptr->flags[j]))
            ++nres_vec;
          if (MAGDATA_ExistY(mptr->flags[j]))
            ++nres_vec;
          if (MAGDATA_ExistZ(mptr->flags[j]))
            ++nres_vec;
        }
    }

  /* check for quick return */
  if (nres_vec == 0)
    return GSL_SUCCESS;

  /*
   * omp_rowidx[thread_id] contains the number of currently filled rows
   * of omp_J[thread_id]. When omp_J[thread_id] is filled, it is folded
   * into the JTJ_vec matrix and then omp_rowidx[thread_id] is reset to 0
   */
  for (i = 0; i < w->max_threads; ++i)
    w->omp_rowidx[i] = 0;

  fprintf(stderr, "\n");

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          gsl_vector *N_gauss = w->gauss_spline_workspace_p[thread_id]->B;
          size_t istart_gauss;

          gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
          gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
          gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);

          gsl_vector_view vx_grad = gsl_matrix_row(w->omp_dX_grad, thread_id);
          gsl_vector_view vy_grad = gsl_matrix_row(w->omp_dY_grad, thread_id);
          gsl_vector_view vz_grad = gsl_matrix_row(w->omp_dZ_grad, thread_id);

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              /* calculate internal Green's functions */
              green_calc_int2(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, w->green_array_p[thread_id]);

              /* calculate non-zero B-splines for the Gauss coefficients for this timestamp */
              gsl_bspline2_eval_basis_nonzero(epoch2year(mptr->t[j]), N_gauss, &istart_gauss, w->gauss_spline_workspace_p[thread_id]);
            }

          /* calculate internal Green's functions for gradient point (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
            {
              green_calc_int2(mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j],
                              &vx_grad.vector, &vy_grad.vector, &vz_grad.vector,
                              w->green_array_p[thread_id]);
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green(N_gauss, istart_gauss, sqrt_wj, &vx.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green(N_gauss, istart_gauss, sqrt_wj, &vy.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green(N_gauss, istart_gauss, sqrt_wj, &vz.vector, &v.vector, w);
                }
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              ++ridx;
            }

#if 0 /*XXX*/
          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              if (mptr->flags[j] & MAGDATA_FLG_DXDT)
                {
                  double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green_SV(t, sqrt_wj, &vx.vector, &v.vector, w);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DYDT)
                {
                  double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green_SV(t, sqrt_wj, &vy.vector, &v.vector, w);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DZDT)
                {
                  double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green_SV(t, sqrt_wj, &vz.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green_grad(t, mptr->ts_ns[j], sqrt_wj, &vx.vector, &vx_grad.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green_grad(t, mptr->ts_ns[j], sqrt_wj, &vy.vector, &vy_grad.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green_grad(t, mptr->ts_ns[j], sqrt_wj, &vz.vector, &vz_grad.vector, &v.vector, w);
                }
            }
#endif /*XXX*/

          /*
           * check if omp_J[thread_id] is full and should be folded into JTJ; the
           * 15 is just some slop to prevent trying to fill rows past the matrix buffer
           * in the loop above
           */
          if (w->omp_rowidx[thread_id] >= w->omp_J[thread_id]->size1 - 15)
            {
              /* fold current matrix block into JTJ_vec, one thread at a time */
              gsl_matrix_view m = gsl_matrix_submatrix(w->omp_J[thread_id], 0, 0, w->omp_rowidx[thread_id], w->p_int);

#pragma omp critical
              {
                gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &m.matrix, 1.0, w->JTJ_vec);
              }

              /* keep cumulative total of rows processed for progress bar */
#pragma omp atomic
              nres_completed += w->omp_rowidx[thread_id];

              /*XXX*/
              gsl_matrix_set_zero(&m.matrix);
              w->omp_rowidx[thread_id] = 0;

#pragma omp critical
              {
                fprintf(stderr, "\t");
                progress_bar(stderr, (double) nres_completed / (double) nres_vec, 70);
              }
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  /* now loop through to see if any rows were not accumulated into JTJ_vec */
  for (i = 0; i < w->max_threads; ++i)
    {
      if (w->omp_rowidx[i] > 0)
        {
          /* accumulate final Green's functions into JTJ_vec */
          gsl_matrix_view m = gsl_matrix_submatrix(w->omp_J[i], 0, 0, w->omp_rowidx[i], w->p_int);
          gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &m.matrix, 1.0, w->JTJ_vec);
          nres_completed += w->omp_rowidx[i];
        }
    }

  fprintf(stderr, "\t");
  progress_bar(stderr, (double) nres_completed / (double) nres_vec, 70);
  fprintf(stderr, "\n");

#if 0
  printsym_octave(w->JTJ_vec, "JTJ_vec");
  exit(1);
#endif

  return s;
}

/*
mfield_vector_green()
  Function to compute sqrt(w) J_1(ridx,:) for a given set of
vector Green's functions, where

                          core                    crust
J_1(ridx,:) = [ -N_k(t(ridx)) B_n^m(r(ridx)) ; -B_n^m(ridx) ]

Inputs: N           - B-spline functions for t(ridx)
        istart      - index of first non-zero B-spline in [0, ncontrol-1]
        sqrt_weight - sqrt(weight) for this measurement
        dB          - vector Green's functions B_n^m, size nnm_tot
        J           - (output) combined vector G = sqrt(w) J_1(ridx,:), size w->p_int
        w           - workspace
*/

static int
mfield_vector_green(const gsl_vector * N, const size_t istart, const double sqrt_weight,
                    const gsl_vector *dB, gsl_vector *J, mfield_workspace *w)
{
  gsl_vector_const_view dB_core = gsl_vector_const_subvector(dB, 0, w->nnm_core);
  size_t k;

  for (k = 0; k < w->params.gauss_spline_order; ++k)
    {
      double Nk = gsl_vector_get(N, k);
      gsl_vector_view v = gsl_vector_subvector(J, (k + istart) * w->nnm_core, w->nnm_core);
      gsl_vector_memcpy_scale(&v.vector, &dB_core.vector, -Nk * sqrt_weight);
    }

  if (w->nnm_crust > 0)
    {
      gsl_vector_const_view dB_crust = gsl_vector_const_subvector(dB, w->nnm_core, w->nnm_crust);
      gsl_vector_view J_crust = gsl_vector_subvector(J, w->p_core, w->nnm_crust);
      gsl_vector_memcpy_scale(&J_crust.vector, &dB_crust.vector, -sqrt_weight);
    }

  return GSL_SUCCESS;
}

/*
mfield_calc_df()
  Compute J^T J matrix and J^T u or J u vector using OpenMP
for speed improvement

The Jacobian is structured as:

J = [ J_1 J_2(x) ]

where J_1 is nres-by-p_int and corresponds to the Gauss coefficients
and J_2(x) is nres-by-p_sparse and corresponds to the Euler angles,
fluxgate calibration, external field and crustal biases.

J^T u = [   J_1^T u  ]
        [ J_2(x)^T u ]

J u = [ J_1 J_2(x) ] [ u_1 ] = [ J_1 u_1 + J_2(x) u_2 ]
                     [ u_2 ]

where u_1 is length p_int and u_2 is length p_sparse.

J^T J = [   J_1^T J_1     J_1^T J_2(x)   ]
        [ J_2(x)^T J_1   J_2(x)^T J_2(x) ]
*/

static int
mfield_calc_df3(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                void *params, gsl_vector * v, gsl_matrix * JTJ)
{
  mfield_workspace *w = (mfield_workspace *) params;
  const size_t gauss_spline_order = gsl_bspline2_order(w->gauss_spline_workspace_p[0]);
  size_t i, j;
  gsl_matrix_view JTJ_int; /* internal field portion of J^T J */
  gsl_matrix_view vJ2TJ1;   /* J_2^T J_1, p_sparse-by-p_int */
  gsl_matrix *J2TJ1 = NULL;
  struct timeval tv0, tv1, tv2, tv3;

  gettimeofday(&tv0, NULL);
  mfield_debug("mfield_calc_df3: entering function...\n");
  mfield_debug("mfield_calc_df3: TransJ = %s...\n",
               TransJ == CblasTrans ? "trans" : "notrans");

  /* initialize outputs to 0 */
  if (v)
    gsl_vector_set_zero(v);

  if (JTJ)
    {
      gsl_matrix_set_zero(JTJ);

      /* copy previously computed vector internal field portion of J^T J
       * (doesn't depend on x) */
      JTJ_int = gsl_matrix_submatrix(JTJ, 0, 0, w->p_int, w->p_int);
      gsl_matrix_tricpy('L', 1, &JTJ_int.matrix, w->JTJ_vec);
    }

  if (w->p_sparse > 0)
    {
      /* calculate sparse part of Jacobian J_2 */
      mfield_calc_J2(x, w->J2, w);

      if (w->J2_csr == NULL)
        w->J2_csr = gsl_spmatrix_alloc_nzmax(w->J2->size1, w->J2->size2, w->J2->nz, GSL_SPMATRIX_CSR);

      mfield_debug("mfield_calc_df3: compressing J_2 to CSR...");
      gettimeofday(&tv2, NULL);
      gsl_spmatrix_csr(w->J2_csr, w->J2);
      gettimeofday(&tv3, NULL);
      mfield_debug("done (%g seconds)\n", time_diff(tv2, tv3));

      mfield_debug("mfield_calc_df3: computing sqrt(W) J_2...");
      gettimeofday(&tv2, NULL);
      gsl_spmatrix_scale_rows(w->J2_csr, w->sqrt_wts_final);
      gettimeofday(&tv3, NULL);
      mfield_debug("done (%g seconds)\n", time_diff(tv2, tv3));

      if (JTJ)
        {
          /* store J_2^T J_2 in lower right portion of JTJ */
          gsl_matrix_view J2TJ2 = gsl_matrix_submatrix(JTJ, w->p_int, w->p_int, w->p_sparse, w->p_sparse);

          mfield_debug("mfield_calc_df3: computing J_2^T W J_2...");
          gettimeofday(&tv2, NULL);
          gsl_spblas_dussyrk(CblasLower, CblasTrans, 1.0, w->J2_csr, 0.0, &J2TJ2.matrix);
          gettimeofday(&tv3, NULL);
          mfield_debug("done (%g seconds)\n", time_diff(tv2, tv3));

          /* set this view which will be used in the loop below to update J_2^T J_1 part of JTJ */
          vJ2TJ1 = gsl_matrix_submatrix(JTJ, w->p_int, 0, w->p_sparse, w->p_int);
          J2TJ1 = &vJ2TJ1.matrix;
        }

      if (TransJ == CblasTrans)
        {
          /* store J_2(x)^T sqrt(W) u in lower part of v */
          gsl_vector_view tmp = gsl_vector_subvector(v, w->p_int, w->p_sparse);
          gsl_spblas_dusmv(CblasTrans, 1.0, w->J2_csr, u, 0.0, &tmp.vector);
        }
      else
        {
          /* store J_2(x) u_2 in v */
          gsl_vector_const_view u2 = gsl_vector_const_subvector(u, w->p_int, w->p_sparse);
          gsl_spblas_dusmv(CblasNoTrans, 1.0, w->J2_csr, &u2.vector, 0.0, v);
        }
    }

  /*
   * omp_rowidx[thread_id] contains the number of currently filled rows
   * of omp_J[thread_id]. When omp_J[thread_id] is filled, it is folded
   * into the JTJ_vec matrix and then omp_rowidx[thread_id] is reset to 0
   */
  for (i = 0; i < w->max_threads; ++i)
    w->omp_rowidx[i] = 0;

  /* loop over satellites */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      /* loop over data for individual satellite */
#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t k;
          double t = mptr->ts[j];       /* use scaled time */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          gsl_vector *N_gauss = w->gauss_spline_workspace_p[thread_id]->B;
          size_t istart_gauss;

          gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
          gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
          gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);

          gsl_vector_view vx_grad = gsl_matrix_row(w->omp_dX_grad, thread_id);
          gsl_vector_view vy_grad = gsl_matrix_row(w->omp_dY_grad, thread_id);
          gsl_vector_view vz_grad = gsl_matrix_row(w->omp_dZ_grad, thread_id);

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              /* compute internal Green's functions for this point */
              green_calc_int2(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, w->green_array_p[thread_id]);

              /* calculate non-zero B-splines for the Gauss coefficients for this timestamp */
              gsl_bspline2_eval_basis_nonzero(epoch2year(mptr->t[j]), N_gauss, &istart_gauss, w->gauss_spline_workspace_p[thread_id]);
            }

          /* calculate internal Green's functions for gradient point (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
            {
              green_calc_int2(mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j],
                              &vx_grad.vector, &vy_grad.vector, &vz_grad.vector,
                              w->green_array_p[thread_id]);
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx);
                  double uj = gsl_vector_get(u, ridx);

                  if (TransJ == CblasTrans)
                    {
#pragma omp critical
                      {
                        mfield_jacobian_J1Tu(N_gauss, istart_gauss, sqrt_wj, uj, &vx.vector, v, w);
                      }
                    }
                  else
                    {
                      mfield_jacobian_J1u(t, sqrt_wj, u, ridx, &vx.vector, v, w);
                    }

                  if (J2TJ1)
                    {
#pragma omp critical
                      {
                        mfield_jacobian_J2TJ1(N_gauss, istart_gauss, sqrt_wj, ridx, &vx.vector, w->J2_csr, J2TJ1, w);
                      }
                    }
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx);
                  double uj = gsl_vector_get(u, ridx);

                  if (TransJ == CblasTrans)
                    {
#pragma omp critical
                      {
                        mfield_jacobian_J1Tu(N_gauss, istart_gauss, sqrt_wj, uj, &vy.vector, v, w);
                      }
                    }
                  else
                    {
                      mfield_jacobian_J1u(t, sqrt_wj, u, ridx, &vy.vector, v, w);
                    }

                  if (J2TJ1)
                    {
#pragma omp critical
                      {
                        mfield_jacobian_J2TJ1(N_gauss, istart_gauss, sqrt_wj, ridx, &vy.vector, w->J2_csr, J2TJ1, w);
                      }
                    }
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx);
                  double uj = gsl_vector_get(u, ridx);

                  if (TransJ == CblasTrans)
                    {
#pragma omp critical
                      {
                        mfield_jacobian_J1Tu(N_gauss, istart_gauss, sqrt_wj, uj, &vz.vector, v, w);
                      }
                    }
                  else
                    {
                      mfield_jacobian_J1u(t, sqrt_wj, u, ridx, &vz.vector, v, w);
                    }

                  if (J2TJ1)
                    {
#pragma omp critical
                      {
                        mfield_jacobian_J2TJ1(N_gauss, istart_gauss, sqrt_wj, ridx, &vz.vector, w->J2_csr, J2TJ1, w);
                      }
                    }
                }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              gsl_vector_view vf = gsl_matrix_row(w->omp_dF, thread_id);
              gsl_vector_view vf_core = gsl_vector_subvector(&vf.vector, 0, w->nnm_core);
              double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx);
              gsl_vector_view Jv = gsl_matrix_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id]++, 0, w->p_int);
              double X, Y, Z, F;

              /* compute internal field model */
              X = mfield_nonlinear_model_int(N_gauss, istart_gauss, &vx.vector, x, w);
              Y = mfield_nonlinear_model_int(N_gauss, istart_gauss, &vy.vector, x, w);
              Z = mfield_nonlinear_model_int(N_gauss, istart_gauss, &vz.vector, x, w);

              /* add apriori model of external (and possibly crustal) field */
              X += mptr->Bx_model[j];
              Y += mptr->By_model[j];
              Z += mptr->Bz_model[j];

              F = gsl_hypot3(X, Y, Z);

              X /= F;
              Y /= F;
              Z /= F;

              /* compute 1/F * (X dX + Y dY + Z dZ) */
              for (k = 0; k < w->nnm_tot; ++k)
                {
                  double dXk = gsl_vector_get(&vx.vector, k);
                  double dYk = gsl_vector_get(&vy.vector, k);
                  double dZk = gsl_vector_get(&vz.vector, k);
                  double val = X * dXk + Y * dYk + Z * dZk;

                  gsl_vector_set(&vf.vector, k, val);
                }

             /* construct row of Jacobian: -sqrt_wt_i * N_k(t_i) * (X dX + Y dY + Z dZ) / F */
             for (k = 0; k < gauss_spline_order; ++k)
               {
                 double Nk = gsl_vector_get(N_gauss, k);
                 gsl_vector_view tmp = gsl_vector_subvector(&Jv.vector, (k + istart_gauss) * w->nnm_core, w->nnm_core);
                 gsl_vector_memcpy_scale(&tmp.vector, &vf_core.vector, -Nk * sqrt_wj);
               }

             if (w->nnm_crust > 0)
               {
                 gsl_vector_view vf_crust = gsl_vector_subvector(&vf.vector, w->nnm_core, w->nnm_crust);
                 gsl_vector_view tmp = gsl_vector_subvector(&Jv.vector, w->p_core, w->nnm_crust);
                 gsl_vector_memcpy_scale(&tmp.vector, &vf_crust.vector, -sqrt_wj);
               }

#pragma omp critical
              {
                mfield_jacobian_F(TransJ, istart_gauss, u, ridx, w->J2_csr, &Jv.vector, v, J2TJ1, w);
              }

              ++ridx;
            }

#if 0
          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              if (mptr->flags[j] & MAGDATA_FLG_DXDT)
                {
                  double wj = gsl_vector_get(w->wts_final, ridx);

                  if (TransJ == CblasTrans)
                    {
#pragma omp critical
                      {
                        mfield_jacobian_SV_JTu(t, mptr->flags[j], wj, u, ridx, &vx.vector, v, w);
                      }
                    }
                  else
                    {
#if 0
#pragma omp critical
                      {
                        mfield_jacobian_Ju(t, mptr->flags[j], wj, u, ridx, &vx.vector,
                                           extidx, dB_ext[0], euler_idx, B_nec_alpha[0],
                                           B_nec_beta[0], B_nec_gamma[0], v, w);
                      }
#endif
                    }

                  ++ridx;
                }

              if (mptr->flags[j] & MAGDATA_FLG_DYDT)
                {
                  double wj = gsl_vector_get(w->wts_final, ridx);

                  if (TransJ == CblasTrans)
                    {
#pragma omp critical
                      {
                        mfield_jacobian_SV_JTu(t, mptr->flags[j], wj, u, ridx, &vy.vector, v, w);
                      }
                    }
                  else
                    {
                    }

                  ++ridx;
                }

              if (mptr->flags[j] & MAGDATA_FLG_DZDT)
                {
                  double wj = gsl_vector_get(w->wts_final, ridx);

                  if (TransJ == CblasTrans)
                    {
#pragma omp critical
                      {
                        mfield_jacobian_SV_JTu(t, mptr->flags[j], wj, u, ridx, &vz.vector, v, w);
                      }
                    }
                  else
                    {
                    }

                  ++ridx;
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              double wj = gsl_vector_get(w->wts_final, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_grad_JTu(t, mptr->ts_ns[j], mptr->flags[j], wj, u, ridx, &vx.vector,
                                             &vx_grad.vector, v, w);
                  }
                }

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              double wj = gsl_vector_get(w->wts_final, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_grad_JTu(t, mptr->ts_ns[j], mptr->flags[j], wj, u, ridx, &vy.vector,
                                             &vy_grad.vector, v, w);
                  }
                }

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              double wj = gsl_vector_get(w->wts_final, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_grad_JTu(t, mptr->ts_ns[j], mptr->flags[j], wj, u, ridx, &vz.vector,
                                             &vz_grad.vector, v, w);
                  }
                }

              ++ridx;
            }
#endif /* 0 */

          /* check if omp_J[thread_id] is full and should be folded into JTJ */
          if (w->omp_rowidx[thread_id] >= w->omp_J[thread_id]->size1)
            {
              if (JTJ)
                {
                  /* accumulate scalar J_int^T J_int into J^T J; it is much faster to do this
                   * with blocks and dsyrk() rather than individual rows with dsyr() */
                  gsl_matrix_view Jm = gsl_matrix_submatrix(w->omp_J[thread_id], 0, 0, w->omp_rowidx[thread_id], w->p_int);

#pragma omp critical
                  {
                    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &Jm.matrix, 1.0, &JTJ_int.matrix);
                  }

                  /*XXX*/
                  /*gsl_matrix_set_zero(&Jm.matrix);*/
                }

              /* reset for new block of rows */
              w->omp_rowidx[thread_id] = 0;
            }
        }
    }

  /* accumulate any last rows of internal field Green's functions */
  for (i = 0; i < w->max_threads; ++i)
    {
      if (JTJ && w->omp_rowidx[i] > 0)
        {
          gsl_matrix_view Jm = gsl_matrix_submatrix(w->omp_J[i], 0, 0, w->omp_rowidx[i], w->p_int);
          gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &Jm.matrix, 1.0, &JTJ_int.matrix);
        }
    }

  /* regularize by adding L^T L to J^T J */
  if (JTJ)
    {
      gsl_spmatrix_add_to_dense(JTJ, w->LTL);
    }

#if 0
  if (u)
    printv_octave(u, "u");

  if (v)
    printv_octave(v, "JTu");

  if (JTJ)
    printsym_octave(JTJ, "JTJ");

  exit(1);
#endif

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_df3: leaving function... (%g seconds)\n", time_diff(tv0, tv1));

  return GSL_SUCCESS;
}

/*
mfield_jacobian_J1Tu()
  Update the J_1^T sqrt(W) u vector with a new row of the Jacobian matrix,
corresponding to a vector residual.

J^T sqrt(W) u = [   J_1^T sqrt(W) u  ]
                [ J_2^(x)T sqrt(W) u ]

The J_2(x)^T sqrt(W) u update is handled via sparse BLAS operations in
mfield_calc_df3().

Inputs: N       - B-spline functions for t(j)
        istart  - index of first non-zero B-spline in [0, ncontrol-1]
        sqrt_wj - sqrt(weight) for this data point - sqrt(W(j,j))
        uj      - input vector element u(j)
        dB      - Green's functions B_n^m for desired vector component of
                  internal SH expansion, nnm_tot-by-1
        JTu     - (output) J_1^T sqrt(W) u vector, length p_int
        w       - workspace
*/

static inline int
mfield_jacobian_J1Tu(const gsl_vector * N, const size_t istart, const double sqrt_wj, const double uj,
                     const gsl_vector * dB, gsl_vector *JTu, const mfield_workspace *w)
{
  const double y = sqrt_wj * uj;
  gsl_vector_const_view dB_core = gsl_vector_const_subvector(dB, 0, w->nnm_core);
  size_t k;

  /* update J^T u */
  for (k = 0; k < w->params.gauss_spline_order; ++k)
    {
      double Nk = gsl_vector_get(N, k);
      gsl_vector_view v = gsl_vector_subvector(JTu, (k + istart) * w->nnm_core, w->nnm_core);
      gsl_blas_daxpy(-Nk * y, &dB_core.vector, &v.vector);
    }

  if (w->nnm_crust > 0)
    {
      gsl_vector_const_view dB_crust = gsl_vector_const_subvector(dB, w->nnm_core, w->nnm_crust);
      gsl_vector_view v = gsl_vector_subvector(JTu, w->p_core, w->nnm_crust);
      gsl_blas_daxpy(-y, &dB_crust.vector, &v.vector);
    }

  return GSL_SUCCESS;
}

/*
mfield_jacobian_Ju()
  Update the J u vector with a new row of the Jacobian matrix,
corresponding to a vector residual.

For this row of J u (specified by ridx), the value is given
by:

(J u)_i = sqrt(w_i) J_1(i,:)^T u

where J_i^T is row i of the Jacobian (1-by-p)

Inputs: t           - scaled timestamp
        sqrt_wj     - sqrt(weight) for this data point
        u           - input vector u, size p
        ridx        - residual index of this row in [0,nres-1]
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_tot-by-1
        Ju          - (output) J u vector
        w           - workspace

Notes:
1) This routine accesses only the element Ju(ridx) so it can be called
by multiple threads with different ridx parameters
*/

static inline int
mfield_jacobian_J1u(const double t, const double sqrt_wj, const gsl_vector * u,
                    const size_t ridx, const gsl_vector * dB_int, gsl_vector *Ju,
                    const mfield_workspace *w)
{
  double *Ju_ptr = gsl_vector_ptr(Ju, ridx);
  double tmp;

  /* update J u */

  if (w->nnm_mf > 0)
    {
      gsl_vector_const_view g_mf = gsl_vector_const_subvector(dB_int, 0, w->nnm_mf);
      gsl_vector_const_view u_mf = gsl_vector_const_subvector(u, 0, w->nnm_mf);

      gsl_blas_ddot(&g_mf.vector, &u_mf.vector, &tmp);
      *Ju_ptr = tmp;
    }

  if (w->nnm_sv > 0)
    {
      gsl_vector_const_view g_sv = gsl_vector_const_subvector(dB_int, 0, w->nnm_sv);
      gsl_vector_const_view u_sv = gsl_vector_const_subvector(u, w->sv_offset, w->nnm_sv);

      gsl_blas_ddot(&g_sv.vector, &u_sv.vector, &tmp);
      *Ju_ptr += t * tmp;
    }

  if (w->nnm_sa > 0)
    {
      gsl_vector_const_view g_sa = gsl_vector_const_subvector(dB_int, 0, w->nnm_sa);
      gsl_vector_const_view u_sa = gsl_vector_const_subvector(u, w->sa_offset, w->nnm_sa);

      gsl_blas_ddot(&g_sa.vector, &u_sa.vector, &tmp);
      *Ju_ptr += 0.5 * t * t * tmp;
    }

  *Ju_ptr *= sqrt_wj;

  return GSL_SUCCESS;
}

/*
mfield_jacobian_F()
  Construct a row of the Jacobian matrix corresponding to
a scalar measurement and update op(J)*u vector. Optionally
update the J_2^T J_1 block according to:

J2TJ1 += J_2(ridx,:) J_1(ridx,:)

where J_2 is passed in via 'J2' and J_1(ridx,:) is
computed in this function.

Inputs: TransJ      - op(J)
        istart      - index in [0,ncontrol-1] of this timestamp
        sqrt_wj     - sqrt(weight) for this data point, sqrt(W(ridx,ridx))
        u           - input vector
        ridx        - index of this row in [0,nres-1]
        J2          - J_2 matrix, CSR format (may be NULL)
        J_int       - row of Jacobian (weighted) for
                      scalar residual and internal Green's functions, p_int-by-1
        v           - (output) updated op(J) u vector
        J2TJ1       - J_2^T J_1 matrix, p_sparse-by-p_int (may be NULL)
        w           - workspace
*/

static inline int
mfield_jacobian_F(CBLAS_TRANSPOSE_t TransJ, const size_t istart, const gsl_vector * u,
                  const size_t ridx, const gsl_spmatrix * J2, const gsl_vector *J_int,
                  gsl_vector *v, gsl_matrix * J2TJ1, const mfield_workspace *w)
{
  const size_t order = w->params.gauss_spline_order;
  gsl_vector_const_view vJ_core = gsl_vector_const_subvector(J_int, istart * w->nnm_core, order * w->nnm_core);

  if (TransJ == CblasTrans)
    {
      /* update J^T u */
      double ui = gsl_vector_get(u, ridx);
      gsl_vector_view vJTu;
      
      vJTu = gsl_vector_subvector(v, istart * w->nnm_core, order * w->nnm_core);
      gsl_blas_daxpy(ui, &vJ_core.vector, &vJTu.vector);

      if (w->nnm_crust > 0)
        {
          gsl_vector_const_view vJ_crust = gsl_vector_const_subvector(J_int, w->p_core, w->nnm_crust);
          vJTu = gsl_vector_subvector(v, w->p_core, w->nnm_crust);
          gsl_blas_daxpy(ui, &vJ_crust.vector, &vJTu.vector);
        }
    }
  else
    {
      /* update (J u)_i = J_int . u(1:pint) */
      double *Ju_ptr = gsl_vector_ptr(v, ridx);
      gsl_vector_const_view z = gsl_vector_const_subvector(u, 0, w->p_int);
      gsl_blas_ddot(J_int, &z.vector, Ju_ptr);
    }

  if (J2TJ1)
    {
      const int *Ap = J2->p;
      const int *Aj = J2->i;
      const double *Ad = J2->data;
      int p;

      for (p = Ap[ridx]; p < Ap[ridx + 1]; ++p)
        {
          gsl_vector_view out = gsl_matrix_subrow(J2TJ1, Aj[p], istart * w->nnm_core, order * w->nnm_core);
          gsl_blas_daxpy(Ad[p], &vJ_core.vector, &out.vector);
        }

      if (w->nnm_crust > 0)
        {
          gsl_vector_const_view vJ_crust = gsl_vector_const_subvector(J_int, w->p_core, w->nnm_crust);

          for (p = Ap[ridx]; p < Ap[ridx + 1]; ++p)
            {
              gsl_vector_view out = gsl_matrix_subrow(J2TJ1, Aj[p], w->p_core, w->nnm_crust);
              gsl_blas_daxpy(Ad[p], &vJ_crust.vector, &out.vector);
            }
        }
    }

  return GSL_SUCCESS;
}

/*
mfield_jacobian_J2TJ1()
  Update J_2^T W J_1 block of J^T W J matrix with:

J2TJ1 += sqrt(W(ridx,ridx)) J_2(ridx,:) J_1(ridx,:)^T

Inputs: N       - non-zero B-spline functions
        istart  - index in [0, ncontrol-1] of first non-zero B-spline
        sqrt_wj - sqrt(weight) for this data point, sqrt(W(ridx,ridx))
        ridx    - index of this row in [0,nres-1]
        dB      - internal Green's functions B_n^m for desired component, length nnm_tot
        J2      - J_2 matrix, CSR format
        J2TJ1   - (output) J_2^T J_1 block of JTJ matrix, size p_sparse-by-p_int
        w       - workspace

Notes:
1) J2 matrix already contains sqrt(W) * J_2, so we need to multiply by another factor
of sqrt(W)
*/

static int
mfield_jacobian_J2TJ1(const gsl_vector * N, const size_t istart, const double sqrt_wj, const size_t ridx,
                      const gsl_vector * dB, const gsl_spmatrix * J2, gsl_matrix * J2TJ1,
                      const mfield_workspace * w)
{
  if (!GSL_SPMATRIX_ISCSR(J2))
    {
      GSL_ERROR("require matrix in CSR format", GSL_EINVAL);
    }
  else
    {
      const size_t order = w->params.gauss_spline_order;
      const int *Ap = J2->p;
      const int *Aj = J2->i;
      const double *Ad = J2->data;
      gsl_vector_const_view dB_core = gsl_vector_const_subvector(dB, 0, w->nnm_core);
      int p;
      size_t k;

      for (k = 0; k < order; ++k)
        {
          double y = gsl_vector_get(N, k) * sqrt_wj;

          for (p = Ap[ridx]; p < Ap[ridx + 1]; ++p)
            {
              gsl_vector_view out = gsl_matrix_subrow(J2TJ1, Aj[p], (istart + k) * w->nnm_core, w->nnm_core);
              gsl_blas_daxpy(-y * Ad[p], &dB_core.vector, &out.vector);
            }
        }

      if (w->nnm_crust > 0)
        {
          for (p = Ap[ridx]; p < Ap[ridx + 1]; ++p)
            {
              gsl_vector_const_view dB_crust = gsl_vector_const_subvector(dB, w->nnm_core, w->nnm_crust);
              gsl_vector_view out = gsl_matrix_subrow(J2TJ1, Aj[p], w->p_core, w->nnm_crust);
              gsl_blas_daxpy(-sqrt_wj * Ad[p], &dB_crust.vector, &out.vector);
            }
        }

      return GSL_SUCCESS;
    }
}
