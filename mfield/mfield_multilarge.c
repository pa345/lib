/*
 * mfield_multilarge.c
 *
 * Contains routines for fitting magnetic field module using gsl_multilarge_nlinear framework
 */

static int mfield_calc_df3(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                           void *params, gsl_vector * v, gsl_matrix * JTJ);
static inline int mfield_jacobian_J1Tu(const double t, const double sqrt_wj, const double uj,
                                       const gsl_vector * dB_int, gsl_vector *JTu, const mfield_workspace *w);
static inline int mfield_jacobian_J1u(const double t, const double sqrt_wt, const gsl_vector * u,
                                      const size_t ridx, const gsl_vector * dB_int, gsl_vector *Ju,
                                      const mfield_workspace *w);
static inline int mfield_jacobian_F(CBLAS_TRANSPOSE_t TransJ, const double t, const double sqrt_wj,
                                    const gsl_vector * u, const size_t ridx,
                                    gsl_vector * dX, gsl_vector * dY, gsl_vector * dZ,
                                    const double B_model[4], const gsl_spmatrix * J2,
                                    gsl_vector *J_int, gsl_vector *v, gsl_matrix * J2TJ1, const mfield_workspace *w);
static int mfield_jacobian_J2TJ1(const double t, const double sqrt_wj, const size_t ridx,
                                 const gsl_vector * dB_int, const gsl_spmatrix * J2, gsl_matrix * J2TJ1,
                                 const mfield_workspace * w);

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
  size_t i, j;
  gsl_matrix_view JTJ_int; /* internal field portion of J^T J */
  gsl_matrix_view vJ2TJ1;   /* J_2^T J_1, p_sparse-by-p_int */
  gsl_matrix *J2TJ1 = NULL;
  gsl_spmatrix *J2_crs = NULL;
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

      mfield_debug("mfield_calc_df3: compressing J_2 to CRS...");
      gettimeofday(&tv2, NULL);
      J2_crs = gsl_spmatrix_crs(w->J2);
      gettimeofday(&tv3, NULL);
      mfield_debug("done (%g seconds)\n", time_diff(tv2, tv3));

      if (JTJ)
        {
          /* store J_2^T J_2 in lower right portion of JTJ */
          gsl_matrix_view J2TJ2 = gsl_matrix_submatrix(JTJ, w->p_int, w->p_int, w->p_sparse, w->p_sparse);

          mfield_debug("mfield_calc_df3: computing J_2^T J_2...");
          gettimeofday(&tv2, NULL);
          gsl_spblas_dusrk(CblasLower, CblasTrans, 1.0, J2_crs, 0.0, &J2TJ2.matrix);
          gettimeofday(&tv3, NULL);
          mfield_debug("done (%g seconds)\n", time_diff(tv2, tv3));

          /* set this view which will be used in the loop below to update J_2^T J_1 part of JTJ */
          vJ2TJ1 = gsl_matrix_submatrix(JTJ, w->p_int, 0, w->p_sparse, w->p_int);
          J2TJ1 = &vJ2TJ1.matrix;
        }

      if (TransJ == CblasTrans)
        {
          /* store J_2(x)^T u in lower part of v */
          gsl_vector_view tmp = gsl_vector_subvector(v, w->p_int, w->p_sparse);
          gsl_spblas_dusmv(CblasTrans, 1.0, J2_crs, u, 0.0, &tmp.vector);
        }
      else
        {
          /* store J_2(x) u_2 in v */
          gsl_vector_const_view u2 = gsl_vector_const_subvector(u, w->p_int, w->p_sparse);
          gsl_spblas_dusmv(CblasNoTrans, 1.0, J2_crs, &u2.vector, 0.0, v);
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

          gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
          gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
          gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);

          gsl_vector_view vx_grad = gsl_matrix_row(w->omp_dX_grad, thread_id);
          gsl_vector_view vy_grad = gsl_matrix_row(w->omp_dY_grad, thread_id);
          gsl_vector_view vz_grad = gsl_matrix_row(w->omp_dZ_grad, thread_id);

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute internal Green's functions for this point */
          green_calc_int2(r, theta, phi, &vx.vector, &vy.vector, &vz.vector,
                          w->green_array_p[thread_id]);

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
              double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx);
              double uj = gsl_vector_get(u, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_J1Tu(t, sqrt_wj, uj, &vx.vector, v, w);
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
                    mfield_jacobian_J2TJ1(t, sqrt_wj, ridx, &vx.vector, J2_crs, J2TJ1, w);
                  }
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx);
              double uj = gsl_vector_get(u, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_J1Tu(t, sqrt_wj, uj, &vy.vector, v, w);
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
                    mfield_jacobian_J2TJ1(t, sqrt_wj, ridx, &vy.vector, J2_crs, J2TJ1, w);
                  }
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx);
              double uj = gsl_vector_get(u, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_J1Tu(t, sqrt_wj, uj, &vz.vector, v, w);
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
                    mfield_jacobian_J2TJ1(t, sqrt_wj, ridx, &vz.vector, J2_crs, J2TJ1, w);
                  }
                }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx);
              gsl_vector_view Jv = gsl_matrix_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id]++, 0, w->p_int);
              double B_int[3];    /* internal field model */
              double B_model[3];  /* a priori model (crustal/external) */
              double B_total[4];  /* internal + external */

              /* compute internal field model */
              B_int[0] = mfield_nonlinear_model_int(t, &vx.vector, x, w);
              B_int[1] = mfield_nonlinear_model_int(t, &vy.vector, x, w);
              B_int[2] = mfield_nonlinear_model_int(t, &vz.vector, x, w);

              /* load apriori model of external (and possibly crustal) field */
              B_model[0] = mptr->Bx_model[j];
              B_model[1] = mptr->By_model[j];
              B_model[2] = mptr->Bz_model[j];

              /* compute total modeled field (internal + external) */
              for (k = 0; k < 3; ++k)
                B_total[k] = B_int[k] + B_model[k];

              B_total[3] = gsl_hypot3(B_total[0], B_total[1], B_total[2]);

#pragma omp critical
              {
                mfield_jacobian_F(TransJ, t, sqrt_wj, u, ridx, &vx.vector, &vy.vector, &vz.vector,
                                  B_total, J2_crs, &Jv.vector, v, J2TJ1, w);
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

  /* regularize by adding L^T L to diag(J^T J) */
  if (JTJ)
    {
      gsl_vector_view v = gsl_matrix_diagonal(JTJ);
      gsl_vector_add(&v.vector, w->LTL);
    }

  if (J2_crs)
    gsl_spmatrix_free(J2_crs);

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
  Update the J_1^T u vector with a new row of the Jacobian matrix,
corresponding to a vector residual.

J^T u = [   J_1^T u  ]
        [ J_2^(x)T u ]

The J_2(x)^T u update is handled via sparse BLAS operations in
mfield_calc_df3().

Inputs: t           - scaled timestamp
        sqrt_wj     - sqrt(weight) for this data point - sqrt(W(j,j))
        uj          - input vector element u(j)
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_max-by-1
        JTu         - (output) J_1^T y vector, length p_int
        w           - workspace
*/

static inline int
mfield_jacobian_J1Tu(const double t, const double sqrt_wj, const double uj,
                     const gsl_vector * dB_int, gsl_vector *JTu, const mfield_workspace *w)
{
  const double y = -sqrt_wj * uj;
  gsl_vector_view v;

  /* update J^T y */

  if (w->nnm_mf > 0)
    {
      gsl_vector_const_view g_mf = gsl_vector_const_subvector(dB_int, 0, w->nnm_mf);

      v = gsl_vector_subvector(JTu, 0, w->nnm_mf);
      gsl_blas_daxpy(y, &g_mf.vector, &v.vector);
    }

  if (w->nnm_sv > 0)
    {
      gsl_vector_const_view g_sv = gsl_vector_const_subvector(dB_int, 0, w->nnm_sv);

      v = gsl_vector_subvector(JTu, w->sv_offset, w->nnm_sv);
      gsl_blas_daxpy(t * y, &g_sv.vector, &v.vector);
    }

  if (w->nnm_sa > 0)
    {
      gsl_vector_const_view g_sa = gsl_vector_const_subvector(dB_int, 0, w->nnm_sa);

      v = gsl_vector_subvector(JTu, w->sa_offset, w->nnm_sa);
      gsl_blas_daxpy(0.5 * t * t * y, &g_sa.vector, &v.vector);
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
                      internal SH expansion, nnm_max-by-1
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
        t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        sqrt_wj     - sqrt(weight) for this data point, sqrt(W(ridx,ridx))
        u           - input vector
        ridx        - index of this row in [0,nres-1]
        dX          - Green's functions for X component
        dY          - Green's functions for Y component
        dZ          - Green's functions for Z component
        B_model     - total model vector
                      B_model[0] = X model
                      B_model[1] = Y model
                      B_model[2] = Z model
                      B_model[3] = F model
        J2          - J_2 matrix, CRS format (may be NULL)
        J_int       - (output) row of Jacobian (weighted) for
                      scalar residual and internal Green's functions, p_int-by-1
        v           - (output) updated op(J) u vector
        J2TJ1       - J_2^T J_1 matrix, p_sparse-by-p_int (may be NULL)
        w           - workspace
*/

static inline int
mfield_jacobian_F(CBLAS_TRANSPOSE_t TransJ, const double t, const double sqrt_wj,
                  const gsl_vector * u, const size_t ridx,
                  gsl_vector * dX, gsl_vector * dY, gsl_vector * dZ,
                  const double B_model[4], const gsl_spmatrix * J2,
                  gsl_vector *J_int, gsl_vector *v, gsl_matrix * J2TJ1, const mfield_workspace *w)
{
  size_t k;
  double b[3];

  /* compute unit vector in model direction */
  for (k = 0; k < 3; ++k)
    b[k] = B_model[k] / B_model[3];

  /* compute (X dX + Y dY + Z dZ) */
  for (k = 0; k < w->nnm_max; ++k)
    {
      double dXk = gsl_vector_get(dX, k);
      double dYk = gsl_vector_get(dY, k);
      double dZk = gsl_vector_get(dZ, k);
      double val = -sqrt_wj * (b[0] * dXk +
                               b[1] * dYk +
                               b[2] * dZk);

      mfield_set_mf(J_int, k, val, w);
      mfield_set_sv(J_int, k, t * val, w);
      mfield_set_sa(J_int, k, 0.5 * t * t * val, w);
    }

  if (TransJ == CblasTrans)
    {
      /* update J^T u */
      double ui = gsl_vector_get(u, ridx);
      gsl_vector_view vJTu = gsl_vector_subvector(v, 0, w->p_int);
      gsl_blas_daxpy(ui, J_int, &vJTu.vector);
    }
  else
    {
      /* update (J u)_i = J_int . u(1:pint) */
      double *Ju_ptr = gsl_vector_ptr(v, ridx);
      gsl_vector_const_view z = gsl_vector_const_subvector(u, 0, w->p_int);
      gsl_blas_ddot(J_int, &z.vector, Ju_ptr);
    }

  if (J2 && J2TJ1)
    {
      const size_t *Ap = J2->p;
      const size_t *Aj = J2->i;
      const double *Ad = J2->data;
      size_t k, p;

      for (k = 0; k < J_int->size; ++k)
        {
          double Jk = gsl_vector_get(J_int, k);

          for (p = Ap[ridx]; p < Ap[ridx + 1]; ++p)
            {
              double *ptr = gsl_matrix_ptr(J2TJ1, Aj[p], k);
              *ptr += Ad[p] * Jk;
            }
        }
    }

  return GSL_SUCCESS;
}

/*
mfield_jacobian_J2TJ1()
  Update J_2^T J_1 block of J^T J matrix with:

J2TJ1 += J_2(ridx,:) J_1(ridx,:)^T

Inputs: t       - scaled timestamp
        sqrt_wj - sqrt(weight) for this data point, sqrt(W(ridx,ridx))
        ridx    - index of this row in [0,nres-1]
        dB_int  - internal Green's functions for desired component, length nnm_max
        J2      - J_2 matrix, CRS format
        J2TJ1   - (output) J_2^T J_1 block of JTJ matrix, size p_sparse-by-p_int
        w       - workspace
*/

static int
mfield_jacobian_J2TJ1(const double t, const double sqrt_wj, const size_t ridx,
                      const gsl_vector * dB_int, const gsl_spmatrix * J2, gsl_matrix * J2TJ1,
                      const mfield_workspace * w)
{
  if (!GSL_SPMATRIX_ISCRS(J2))
    {
      GSL_ERROR("require matrix in CRS format", GSL_EINVAL);
    }
  else
    {
      const size_t *Ap = J2->p;
      const size_t *Aj = J2->i;
      const double *Ad = J2->data;
      size_t p;

      if (w->nnm_mf > 0)
        {
          gsl_vector_const_view in = gsl_vector_const_subvector(dB_int, 0, w->nnm_mf);

          for (p = Ap[ridx]; p < Ap[ridx + 1]; ++p)
            {
              gsl_vector_view out = gsl_matrix_subrow(J2TJ1, Aj[p], 0, w->nnm_mf);
              gsl_blas_daxpy(-Ad[p], &in.vector, &out.vector);
            }
        }

      if (w->nnm_sv > 0)
        {
          gsl_vector_const_view in = gsl_vector_const_subvector(dB_int, 0, w->nnm_sv);

          for (p = Ap[ridx]; p < Ap[ridx + 1]; ++p)
            {
              gsl_vector_view out = gsl_matrix_subrow(J2TJ1, Aj[p], w->sv_offset, w->nnm_sv);
              gsl_blas_daxpy(-t * Ad[p], &in.vector, &out.vector);
            }
        }

      if (w->nnm_sa > 0)
        {
          gsl_vector_const_view in = gsl_vector_const_subvector(dB_int, 0, w->nnm_sa);

          for (p = Ap[ridx]; p < Ap[ridx + 1]; ++p)
            {
              gsl_vector_view out = gsl_matrix_subrow(J2TJ1, Aj[p], w->sa_offset, w->nnm_sa);
              gsl_blas_daxpy(-0.5 * t * t * Ad[p], &in.vector, &out.vector);
            }
        }

      return GSL_SUCCESS;
    }
}
