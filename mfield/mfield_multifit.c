/*
 * mfield_multifit.c
 *
 * Contains routines for fitting magnetic field module using gsl_multifit_nlinear framework
 */

static int mfield_calc_f(const gsl_vector *x, void *params, gsl_vector *f);
static int mfield_calc_Wf(const gsl_vector *x, void *params, gsl_vector *f);
static int mfield_calc_J2(const gsl_vector *x, gsl_spmatrix *J, mfield_workspace *w);
static int mfield_calc_df(const gsl_vector *x, void *params, gsl_matrix *J);
static int mfield_nonlinear_model(const int res_flag, const gsl_vector * x, const magdata * mptr, const size_t src_idx,
                                  const size_t data_idx, const size_t thread_id, double B_model[3], mfield_workspace *w);
static double mfield_nonlinear_model_int(const gsl_vector * N, const size_t istart, const gsl_vector *dB,
                                         const gsl_vector *g, const mfield_workspace *w);
static int mfield_nonlinear_model_SV(const gsl_vector * x, const magdata * mptr, const size_t idx,
                                     const size_t thread_id, double dBdt_model[3], mfield_workspace *w);
static double mfield_nonlinear_model_int_SV(const gsl_vector * dN, const size_t istart, const gsl_vector *dB,
                                            const gsl_vector *g, const mfield_workspace *w);
static inline int jacobian_row_int(const gsl_vector * N, const size_t istart, const gsl_vector * dB,
                                   gsl_vector * J, mfield_workspace * w);
static inline int jacobian_row_int_SV(const gsl_vector * dN, const size_t istart, const gsl_vector * dB,
                                      gsl_vector * J, mfield_workspace * w);
static inline int jacobian_row_int_grad(const double t, const double t2, const gsl_vector * dB, const gsl_vector * dB2,
                                        gsl_vector * J, mfield_workspace * w);

/*
mfield_calc_f()
  Construct residual vector f(x) using OpenMP parallel
processing to compute all the Green's functions quickly.

Inputs: x      - model coefficients
        params - parameters
        f      - (output) residual vector

Notes:
1) For the histograms, w->wts_final must be initialized prior
to calling this function
*/

static int
mfield_calc_f(const gsl_vector *x, void *params, gsl_vector *f)
{
  int s = GSL_SUCCESS;
  mfield_workspace *w = (mfield_workspace *) params;
  const mfield_parameters *mparams = &(w->params);
  size_t i, j;
  struct timeval tv0, tv1;

  mfield_debug("mfield_calc_f: entering function (||x|| = %g)...\n", gsl_blas_dnrm2(x));
  gettimeofday(&tv0, NULL);

  gsl_vector_set_zero(f);

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      const double *t_year = w->data_workspace_p->t_year[i];

      /* alignment parameters */
      int fit_align = mparams->fit_align && (mptr->global_flags & MAGDATA_GLOBFLG_ALIGN);
      size_t align_ncontrol = fit_align ? gsl_bspline2_ncontrol(w->align_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)]) : 0;
      size_t align_idx = w->align_offset + w->offset_align[i];

      /* fluxgate calibration */
      int fit_fluxcal = mparams->fit_fluxcal && (mptr->global_flags & MAGDATA_GLOBFLG_FLUXCAL);
      size_t fluxcal_ncontrol = fit_fluxcal ? gsl_bspline2_ncontrol(w->fluxcal_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)]) : 0;
      size_t fluxcal_idx = w->fluxcal_offset + w->offset_fluxcal[i];

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          double B_model[3];            /* internal + external */
          double B_obs[3];              /* observation vector NEC frame */
          double B_model_ns[3];         /* N/S internal + external */
          double B_obs_ns[3];           /* N/S observation vector NEC frame */
          double dBdt_model[3];         /* dB/dt (SV of internal model) */
          double dBdt_obs[3];           /* SV observation vector (NEC frame) */
          double F_obs;                 /* scalar field measurement */
          gsl_bspline2_workspace *align_spline_p = fit_align ? w->align_spline_workspace_p[CIDX2(i, w->nsat, thread_id, w->max_threads)] : NULL;
          gsl_bspline2_workspace *fluxcal_spline_p = fit_fluxcal ? w->fluxcal_spline_workspace_p[CIDX2(i, w->nsat, thread_id, w->max_threads)] : NULL;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute vector model for this residual */
          mfield_nonlinear_model(0, x, mptr, i, j, thread_id, B_model, w);

          if (!gsl_finite(B_model[0]) || !gsl_finite(B_model[1]) || !gsl_finite(B_model[2]))
            {
              fprintf(stderr, "HERE\n");
              exit(1);
            }

          /* compute vector SV model for this residual */
          if (MAGDATA_FitMF(mptr->flags[j]) && (mptr->flags[j] & (MAGDATA_FLG_DXDT | MAGDATA_FLG_DYDT | MAGDATA_FLG_DZDT)))
            {
              mfield_nonlinear_model_SV(x, mptr, j, thread_id, dBdt_model, w);
              dBdt_obs[0] = mptr->dXdt_nec[j];
              dBdt_obs[1] = mptr->dYdt_nec[j];
              dBdt_obs[2] = mptr->dZdt_nec[j];
            }

          /* compute vector model for gradient residual (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
            mfield_nonlinear_model(1, x, mptr, i, j, thread_id, B_model_ns, w);

          if (fit_align)
            {
              /*
               * get the alignment parameters for this satellite and time period,
               * and apply rotation
               */
              double align_data[ALIGN_P];
              gsl_vector_view align_params = gsl_vector_view_array(align_data, ALIGN_P);

              /* alignment control points for this dataset */
              gsl_vector_const_view v1 = gsl_vector_const_subvector(x, align_idx, ALIGN_P * align_ncontrol);
              gsl_matrix_const_view align_control_pts = gsl_matrix_const_view_vector(&v1.vector, ALIGN_P, align_ncontrol);

              double *q = &(mptr->q[4*j]);
              double B_vfm[4];

              /* compute alignment parameters for this timestamp */
              gsl_bspline2_vector_eval(t_year[j], &align_control_pts.matrix, &align_params.vector, align_spline_p);

              B_vfm[0] = mptr->Bx_vfm[j];
              B_vfm[1] = mptr->By_vfm[j];
              B_vfm[2] = mptr->Bz_vfm[j];

              if (fit_fluxcal)
                {
                  int status;
                  double cal_data[FLUXCAL_P];
                  gsl_vector_view cal_params = gsl_vector_view_array(cal_data, FLUXCAL_P);

                  /* fluxgate calibration control points for this dataset */
                  gsl_vector_const_view v2 = gsl_vector_const_subvector(x, fluxcal_idx, FLUXCAL_P * fluxcal_ncontrol);
                  gsl_matrix_const_view fluxcal_control_pts = gsl_matrix_const_view_vector(&v2.vector, FLUXCAL_P, fluxcal_ncontrol);

                  gsl_bspline2_vector_eval(t_year[j], &fluxcal_control_pts.matrix, &cal_params.vector, fluxcal_spline_p);

                  status = mfield_fluxcal_apply_datum(&cal_params.vector, B_vfm, B_vfm);
                  if (status)
                    {
                      /* XXX non-orthogonality angles caused a singularity in the P^{-1} matrix */
                      B_vfm[0] = mptr->Bx_vfm[j];
                      B_vfm[1] = mptr->By_vfm[j];
                      B_vfm[2] = mptr->Bz_vfm[j];
                    }
                }

              /* rotate VFM vector to NEC, B_obs = R_q R_3(alpha) B_VFM(c) */
              mfield_align_vfm2nec(&align_params.vector, q, B_vfm, B_obs, w->att_workspace_p[i], w);

              if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                    MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
                {
                  double *q_ns = &(mptr->q_ns[4*j]);
                  double B_vfm_ns[3];

                  B_vfm_ns[0] = mptr->Bx_vfm_ns[j];
                  B_vfm_ns[1] = mptr->By_vfm_ns[j];
                  B_vfm_ns[2] = mptr->Bz_vfm_ns[j];

                  mfield_align_vfm2nec(&align_params.vector, q_ns, B_vfm_ns, B_obs_ns, w->att_workspace_p[i], w);
                }
            }
          else
            {
              assert(fit_fluxcal == 0);

              /* use supplied NEC vector */
              B_obs[0] = mptr->Bx_nec[j];
              B_obs[1] = mptr->By_nec[j];
              B_obs[2] = mptr->Bz_nec[j];

              if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                    MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
                {
                  B_obs_ns[0] = mptr->Bx_nec_ns[j];
                  B_obs_ns[1] = mptr->By_nec_ns[j];
                  B_obs_ns[2] = mptr->Bz_nec_ns[j];
                }
            }

          if (fit_align && fit_fluxcal && (mptr->flags[j] & (MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z)))
            {
              /* F_obs = || R_q R_3(alpha) B_vfm(c) || */
              F_obs = gsl_hypot3(B_obs[0], B_obs[1], B_obs[2]);
            }
          else
            {
              F_obs = mptr->F[j];
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              gsl_vector_set(f, ridx++, B_obs[0] - B_model[0]);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              gsl_vector_set(f, ridx++, B_obs[1] - B_model[1]);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              gsl_vector_set(f, ridx++, B_obs[2] - B_model[2]);
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double F_mod = gsl_hypot3(B_model[0], B_model[1], B_model[2]);
              gsl_vector_set(f, ridx++, F_obs - F_mod);
            }

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              if (mptr->flags[j] & MAGDATA_FLG_DXDT)
                {
                  gsl_vector_set(f, ridx++, dBdt_obs[0] - dBdt_model[0]);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DYDT)
                {
                  gsl_vector_set(f, ridx++, dBdt_obs[1] - dBdt_model[1]);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DZDT)
                {
                  gsl_vector_set(f, ridx++, dBdt_obs[2] - dBdt_model[2]);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              gsl_vector_set(f, ridx++, (B_model_ns[0] - B_model[0]) - (B_obs_ns[0] - B_obs[0]));
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              gsl_vector_set(f, ridx++, (B_model_ns[1] - B_model[1]) - (B_obs_ns[1] - B_obs[1]));
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              gsl_vector_set(f, ridx++, (B_model_ns[2] - B_model[2]) - (B_obs_ns[2] - B_obs[2]));
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    }

  if (mparams->regularize && !mparams->synth_data)
    {
      /* store L^T*x in bottom of f for regularization */
      gsl_vector_view v = gsl_vector_subvector(f, w->nres, w->p);
      gsl_spblas_dusmv(CblasTrans, 1.0, w->L, x, 0.0, &v.vector);
    }

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_f: leaving function (%g seconds, ||f|| = %g)\n",
               time_diff(tv0, tv1), gsl_blas_dnrm2(f));

#if 0
  printv_octave(x, "x");
  printv_octave(f, "f");
  exit(1);
#endif

  return s;
}

/*
mfield_calc_Wf()
  Compute weighted residuals:

f~(x) = sqrt(W) f(x)

Inputs: x      - model parameters
        params - parameters
        f      - (output) f~(x)
*/

static int
mfield_calc_Wf(const gsl_vector *x, void *params, gsl_vector *f)
{
  int s;
  mfield_workspace *w = (mfield_workspace *) params;
  size_t i;
  struct timeval tv0, tv1;

  mfield_debug("mfield_calc_Wf: entering function...\n");
  gettimeofday(&tv0, NULL);

  s = mfield_calc_f(x, params, f);
  if (s)
    return s;

  for (i = 0; i < w->nres; ++i)
    {
      double swi = gsl_vector_get(w->sqrt_wts_final, i);
      double *fi = gsl_vector_ptr(f, i);

      *fi *= swi;
    }

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_Wf: leaving function (%g seconds, ||sqrt(W) f|| = %g)\n",
               time_diff(tv0, tv1), gsl_blas_dnrm2(f));

#if 0
  printv_octave(f, "f");
  printv_octave(x, "x");
  printv_octave(w->wts_final, "wts");
  exit(1);
#endif

  return GSL_SUCCESS;
}

/*
mfield_calc_fvv()
  Construct D_v^2 f(x) using OpenMP parallel
processing to compute all the Green's functions quickly.

Inputs: x      - model coefficients, length p
        v      - velocity vector, length p
        params - parameters
        f      - (output) residual vector

Notes:
1) For the histograms, w->wts_final must be initialized prior
to calling this function
*/

static int
mfield_calc_fvv(const gsl_vector *x, const gsl_vector * v, void *params, gsl_vector *fvv)
{
  int s = GSL_SUCCESS;
#if 1
  (void) x;
  (void) v;
  (void) params;
  (void) fvv;
#else
  mfield_workspace *w = (mfield_workspace *) params;
  const mfield_parameters *mparams = &(w->params);
  size_t i, j;
  struct timeval tv0, tv1;

  mfield_debug("mfield_calc_fvv: entering function...\n");
  gettimeofday(&tv0, NULL);

  gsl_vector_set_zero(fvv);

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      int fit_euler = mparams->fit_euler && (mptr->global_flags & MAGDATA_GLOBFLG_ALIGN);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          double Dvv[3] = { 0.0, 0.0, 0.0 }; /* second directional derivative of vector residuals */

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (fit_euler)
            {
              /*
               * get the Euler angles for this satellite and time period,
               * and apply rotation
               */
              size_t euler_idx = mfield_euler_idx(i, mptr->t[j], w);
              double alpha = gsl_vector_get(x, euler_idx);
              double beta = gsl_vector_get(x, euler_idx + 1);
              double gamma = gsl_vector_get(x, euler_idx + 2);
              double v_alpha = gsl_vector_get(v, euler_idx);
              double v_beta = gsl_vector_get(v, euler_idx + 1);
              double v_gamma = gsl_vector_get(v, euler_idx + 2);
              double *q = &(mptr->q[4*j]);
              double B_vfm[3], B_nec[3];
              size_t k;

              B_vfm[0] = mptr->Bx_vfm[j];
              B_vfm[1] = mptr->By_vfm[j];
              B_vfm[2] = mptr->Bz_vfm[j];

              /* compute diagonal terms first */

              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV2_ALPHA, alpha, beta, gamma, q, B_vfm, B_nec);
              for (k = 0; k < 3; ++k)
                Dvv[k] += v_alpha * v_alpha * B_nec[k];

              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV2_BETA, alpha, beta, gamma, q, B_vfm, B_nec);
              for (k = 0; k < 3; ++k)
                Dvv[k] += v_beta * v_beta * B_nec[k];

              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV2_GAMMA, alpha, beta, gamma, q, B_vfm, B_nec);
              for (k = 0; k < 3; ++k)
                Dvv[k] += v_gamma * v_gamma * B_nec[k];

              /* now the cross terms */

              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_ALPHA | EULER_FLG_DERIV_BETA, alpha, beta, gamma, q, B_vfm, B_nec);
              for (k = 0; k < 3; ++k)
                Dvv[k] += 2.0 * v_alpha * v_beta * B_nec[k];

              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_ALPHA | EULER_FLG_DERIV_GAMMA, alpha, beta, gamma, q, B_vfm, B_nec);
              for (k = 0; k < 3; ++k)
                Dvv[k] += 2.0 * v_alpha * v_gamma * B_nec[k];

              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_BETA | EULER_FLG_DERIV_GAMMA, alpha, beta, gamma, q, B_vfm, B_nec);
              for (k = 0; k < 3; ++k)
                Dvv[k] += 2.0 * v_beta * v_gamma * B_nec[k];
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              gsl_vector_set(fvv, ridx, Dvv[0]);
              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              gsl_vector_set(fvv, ridx, Dvv[1]);
              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              gsl_vector_set(fvv, ridx, Dvv[2]);
              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              const double t = mptr->ts[j];       /* use scaled time */
              const double t2 = t * t;
              const double t3 = t * t2;
              const double t4 = t2 * t2;
              double B_model[4], b_model[3];
              double Dvvf = 0.0;
              gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
              gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
              gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);
              size_t k, kp;

              /* compute vector model for this residual (this call will fill in vx, vy, vz) */
              mfield_nonlinear_model(0, x, mptr, i, j, thread_id, B_model, w);
              B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

              /* b_model = B_model / || B_model || */
              for (k = 0; k < 3; ++k)
                b_model[k] = B_model[k] / B_model[3];

#if 1
              for (k = 0; k < w->nnm_max; ++k)
                {
                  double v_mf = mfield_get_mf(v, k, w);
                  double v_sv = mfield_get_sv(v, k, w);
                  double v_sa = mfield_get_sa(v, k, w);
                  double dXk = gsl_vector_get(&vx.vector, k);
                  double dYk = gsl_vector_get(&vy.vector, k);
                  double dZk = gsl_vector_get(&vz.vector, k);
                  double term0 = b_model[0] * dXk + b_model[1] * dYk + b_model[2] * dZk; /* b_model . dB^{int}_{nm} */
                  double tmp[3];

                  /* tmp = (b_model . dB^{int}_{nm}) * b_model + dB^{int}_{nm} */
                  tmp[0] = term0 * b_model[0] + dXk;
                  tmp[1] = term0 * b_model[1] + dYk;
                  tmp[2] = term0 * b_model[2] + dZk;

                  for (kp = 0; kp < w->nnm_max; ++kp)
                    {
                      double vp_mf = mfield_get_mf(v, kp, w);
                      double vp_sv = mfield_get_sv(v, kp, w);
                      double vp_sa = mfield_get_sa(v, kp, w);
                      double dXkp = gsl_vector_get(&vx.vector, kp);
                      double dYkp = gsl_vector_get(&vy.vector, kp);
                      double dZkp = gsl_vector_get(&vz.vector, kp);
                      double xi = tmp[0] * dXkp + tmp[1] * dYkp + tmp[2] * dZkp; /* ||B_model|| * (d/dg_{nm}) (d/dg_{n'm'}) f_i */

                      Dvvf += xi * (v_mf * vp_mf +
                                    t2 * v_sv * vp_sv +
                                    0.25 * t4 * v_sa * vp_sa +
                                    t * (v_mf * vp_sv + v_sv * vp_mf) +
                                    0.5 * t2 * (v_mf * vp_sa + v_sa * vp_mf) +
                                    0.5 * t3 * (v_sv * vp_sa + v_sa * vp_sv));
                    }
                }
#endif

              Dvvf /= B_model[3];

              gsl_vector_set(fvv, ridx, Dvvf);
              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              ++ridx;
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    }

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_fvv: leaving function (%g seconds)\n", time_diff(tv0, tv1));
#endif

  return s;
}

/*
mfield_calc_Wfvv()
  Compute weighted residuals:

fvv~(x) = sqrt(W) fvv(x)

Inputs: x      - model parameters
        v      - velocity vector
        params - parameters
        fvv    - (output) fvv~(x)
*/

static int
mfield_calc_Wfvv(const gsl_vector *x, const gsl_vector *v, void *params, gsl_vector *fvv)
{
  int s;
  mfield_workspace *w = (mfield_workspace *) params;
  size_t i;
  struct timeval tv0, tv1;

  mfield_debug("mfield_calc_Wfvv: entering function...\n");
  gettimeofday(&tv0, NULL);

  s = mfield_calc_fvv(x, v, params, fvv);
  if (s)
    return s;

  for (i = 0; i < w->nres; ++i)
    {
      double wi = gsl_vector_get(w->wts_final, i);
      double *fi = gsl_vector_ptr(fvv, i);

      *fi *= sqrt(wi);
    }

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_Wfvv: leaving function (%g seconds)\n", time_diff(tv0, tv1));

  return GSL_SUCCESS;
}

/*
mfield_calc_J2()
  Compute Jacobian matrix J_2(x) containing Euler angles,
fluxgate calibration, and external field parameters. This matrix
has significant sparse structure.

Inputs: x - parameter vector, length p
        J - (output) J_2(x) matrix, nres-by-p_sparse
        w - workspace
*/

static int
mfield_calc_J2(const gsl_vector *x, gsl_spmatrix *J, mfield_workspace *w)
{
  int s = GSL_SUCCESS;
  gsl_vector_const_view y = gsl_vector_const_subvector(x, w->p_int, w->p - w->p_int);
  size_t i, j;
  struct timeval tv0, tv1;
  int parallel;

  mfield_debug("mfield_calc_J2: entering function (J_2: %zu-by-%zu)...\n", J->size1, J->size2);
  gettimeofday(&tv0, NULL);

  if (J->nz == 0)
    {
      /*
       * on first call, use single thread to determine all non-zero matrix elements;
       * on subsequent calls, we can use multiple threads since the sparsity pattern
       * will be fixed and each thread will access a matrix element in a fixed tree
       * structure, without worrying about race conditions
       */
      parallel = 0;
    }
  else
    {
      parallel = 1;
      J->spflags |= GSL_SPMATRIX_FLG_FIXED;
    }

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      const double *t_year = w->data_workspace_p->t_year[i];

      /* Euler angles */
      int fit_align = w->params.fit_align && (mptr->global_flags & MAGDATA_GLOBFLG_ALIGN);
      size_t align_ncontrol = fit_align ? gsl_bspline2_ncontrol(w->align_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)]) : 0;
      size_t align_idx = w->align_offset - w->p_int + w->offset_align[i];

      /* fluxgate calibration */
      int fit_fluxcal = w->params.fit_fluxcal && (mptr->global_flags & MAGDATA_GLOBFLG_FLUXCAL);
      size_t fluxcal_ncontrol = fit_fluxcal ? gsl_bspline2_ncontrol(w->fluxcal_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)]) : 0;
      size_t fluxcal_idx = w->fluxcal_offset - w->p_int + w->offset_fluxcal[i];

#pragma omp parallel for private(j) if (parallel)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          size_t k;

          /* Euler angle parameters */
          gsl_bspline2_workspace *align_spline_p = fit_align ? w->align_spline_workspace_p[CIDX2(i, w->nsat, thread_id, w->max_threads)] : NULL;
          double align_data[ALIGN_P];
          gsl_vector_view align_params = gsl_vector_view_array(align_data, ALIGN_P);
          gsl_vector *N_align = fit_align ? align_spline_p->B : NULL;
          size_t istart_align;
          double B_vfm[4], B_nec_alpha[3], B_nec_beta[3], B_nec_gamma[3];
          double B_nec_alpha_ns[3], B_nec_beta_ns[3], B_nec_gamma_ns[3];

          /* fluxgate calibration parameters */
          gsl_bspline2_workspace *fluxcal_spline_p = fit_fluxcal ? w->fluxcal_spline_workspace_p[CIDX2(i, w->nsat, thread_id, w->max_threads)] : NULL;
          double cal_data[FLUXCAL_P];
          gsl_vector_view cal_params = gsl_vector_view_array(cal_data, FLUXCAL_P);
          double jac_fluxcal_data[3 * FLUXCAL_P];
          gsl_matrix_view jac_fluxcal = gsl_matrix_view_array(jac_fluxcal_data, 3, FLUXCAL_P);
          double jac_fluxcal_align_data[3 * FLUXCAL_P];
          gsl_matrix_view jac_fluxcal_align = gsl_matrix_view_array(jac_fluxcal_align_data, 3, FLUXCAL_P);
          gsl_vector *N_fluxcal = fit_fluxcal ? fluxcal_spline_p->B : NULL;
          size_t istart_fluxcal;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute alignment parameter derivatives of residual vector */
          if (fit_align)
            {
              double *q = &(mptr->q[4*j]);
              gsl_vector_const_view v1 = gsl_vector_const_subvector(&y.vector, align_idx, ALIGN_P * align_ncontrol);
              gsl_matrix_const_view align_control_pts = gsl_matrix_const_view_vector(&v1.vector, ALIGN_P, align_ncontrol);

              /* compute Euler angles for this timestamp */
              gsl_bspline2_vector_eval(t_year[j], &align_control_pts.matrix, &align_params.vector, align_spline_p);

              /* evaluate non-zero basis splines for time t */
              gsl_bspline2_eval_basis_nonzero(t_year[j], N_align, &istart_align, align_spline_p);

              /* get vector in VFM frame */
              B_vfm[0] = mptr->Bx_vfm[j];
              B_vfm[1] = mptr->By_vfm[j];
              B_vfm[2] = mptr->Bz_vfm[j];

              if (fit_fluxcal)
                {
                  int status;

                  /* fluxgate calibration control points for this dataset */
                  gsl_vector_const_view v2 = gsl_vector_const_subvector(&y.vector, fluxcal_idx, FLUXCAL_P * fluxcal_ncontrol);
                  gsl_matrix_const_view fluxcal_control_pts = gsl_matrix_const_view_vector(&v2.vector, FLUXCAL_P, fluxcal_ncontrol);

                  gsl_bspline2_vector_eval(t_year[j], &fluxcal_control_pts.matrix, &cal_params.vector, fluxcal_spline_p);

                  /* compute jac_fluxcal := d/dm B_vfm(m) */
                  status = mfield_fluxcal_jac(&cal_params.vector, B_vfm, &jac_fluxcal.matrix);
                  if (status)
                    {
                      /* XXX non-orthogonality angles caused a singularity in the P^{-1} matrix */
                      gsl_matrix_set_zero(&jac_fluxcal.matrix);
                    }

                  /* apply R_q * R_3(alpha) to d/dm B_vfm(m) vectors */
                  mfield_align_matrix_vfm2nec(&align_params.vector, q, &jac_fluxcal.matrix, &jac_fluxcal_align.matrix,
                                              w->att_workspace_p[i], w);

                  /* compute B_vfm := P^{-1} S (B_vfm - O) */
                  status = mfield_fluxcal_apply_datum(&cal_params.vector, B_vfm, B_vfm);
                  if (status)
                    {
                      /* XXX non-orthogonality angles caused a singularity in the P^{-1} matrix */
                      B_vfm[0] = mptr->Bx_vfm[j];
                      B_vfm[1] = mptr->By_vfm[j];
                      B_vfm[2] = mptr->Bz_vfm[j];
                    }

                  /* evaluate non-zero basis splines for time t */
                  gsl_bspline2_eval_basis_nonzero(t_year[j], N_fluxcal, &istart_fluxcal, fluxcal_spline_p);
                }

              /* compute alpha derivative of: R_q R_3 B_vfm */
              mfield_align_deriv_vfm2nec(0, &align_params.vector, q, B_vfm, B_nec_alpha, w->att_workspace_p[i], w);

              /* compute beta derivative of: R_q R_3 B_vfm */
              mfield_align_deriv_vfm2nec(1, &align_params.vector, q, B_vfm, B_nec_beta, w->att_workspace_p[i], w);

              /* compute gamma derivative of: R_q R_3 B_vfm */
              mfield_align_deriv_vfm2nec(2, &align_params.vector, q, B_vfm, B_nec_gamma, w->att_workspace_p[i], w);

              if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                    MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
                {
                  q = &(mptr->q_ns[4*j]);

                  /* get vector in VFM frame */
                  B_vfm[0] = mptr->Bx_vfm_ns[j];
                  B_vfm[1] = mptr->By_vfm_ns[j];
                  B_vfm[2] = mptr->Bz_vfm_ns[j];

                  /* compute alpha derivative of: R_q R_3 B_vfm_ns */
                  mfield_align_deriv_vfm2nec(0, &align_params.vector, q, B_vfm, B_nec_alpha_ns, w->att_workspace_p[i], w);

                  /* compute beta derivative of: R_q R_3 B_vfm_ns */
                  mfield_align_deriv_vfm2nec(1, &align_params.vector, q, B_vfm, B_nec_beta_ns, w->att_workspace_p[i], w);

                  /* compute gamma derivative of: R_q R_3 B_vfm_ns */
                  mfield_align_deriv_vfm2nec(2, &align_params.vector, q, B_vfm, B_nec_gamma_ns, w->att_workspace_p[i], w);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              /* check if fitting alignment parameters to this data point */
              if (fit_align && MAGDATA_FitAlign(mptr->flags[j]))
                {
                  for (k = 0; k < align_spline_p->spline_order; ++k)
                    {
                      double Nk = gsl_vector_get(N_align, k);

                      gsl_spmatrix_set(J, ridx, align_idx + CIDX2(0, ALIGN_P, istart_align + k, align_ncontrol), Nk * B_nec_alpha[0]);
                      gsl_spmatrix_set(J, ridx, align_idx + CIDX2(1, ALIGN_P, istart_align + k, align_ncontrol), Nk * B_nec_beta[0]);
                      gsl_spmatrix_set(J, ridx, align_idx + CIDX2(2, ALIGN_P, istart_align + k, align_ncontrol), Nk * B_nec_gamma[0]);
                    }

                  if (fit_fluxcal)
                    {
                      size_t idx;

                      for (k = 0; k < fluxcal_spline_p->spline_order; ++k)
                        {
                          double Nk = gsl_vector_get(N_fluxcal, k);

                          for (idx = 0; idx < FLUXCAL_P; ++idx)
                            {
                              gsl_spmatrix_set(J, ridx, fluxcal_idx + CIDX2(idx, FLUXCAL_P, istart_fluxcal + k, fluxcal_ncontrol),
                                               Nk * gsl_matrix_get(&jac_fluxcal_align.matrix, 0, idx));
                            }
                        }
                    }
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              /* check if fitting Euler angles to this data point */
              if (fit_align && MAGDATA_FitAlign(mptr->flags[j]))
                {
                  for (k = 0; k < align_spline_p->spline_order; ++k)
                    {
                      double Nk = gsl_vector_get(N_align, k);

                      gsl_spmatrix_set(J, ridx, align_idx + CIDX2(0, ALIGN_P, istart_align + k, align_ncontrol), Nk * B_nec_alpha[1]);
                      gsl_spmatrix_set(J, ridx, align_idx + CIDX2(1, ALIGN_P, istart_align + k, align_ncontrol), Nk * B_nec_beta[1]);
                      gsl_spmatrix_set(J, ridx, align_idx + CIDX2(2, ALIGN_P, istart_align + k, align_ncontrol), Nk * B_nec_gamma[1]);
                    }

                  if (fit_fluxcal)
                    {
                      size_t idx;

                      for (k = 0; k < fluxcal_spline_p->spline_order; ++k)
                        {
                          double Nk = gsl_vector_get(N_fluxcal, k);

                          for (idx = 0; idx < FLUXCAL_P; ++idx)
                            {
                              gsl_spmatrix_set(J, ridx, fluxcal_idx + CIDX2(idx, FLUXCAL_P, istart_fluxcal + k, fluxcal_ncontrol),
                                               Nk * gsl_matrix_get(&jac_fluxcal_align.matrix, 1, idx));
                            }
                        }
                    }
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              /* check if fitting Euler angles to this data point */
              if (fit_align && MAGDATA_FitAlign(mptr->flags[j]))
                {
                  for (k = 0; k < align_spline_p->spline_order; ++k)
                    {
                      double Nk = gsl_vector_get(N_align, k);

                      gsl_spmatrix_set(J, ridx, align_idx + CIDX2(0, ALIGN_P, istart_align + k, align_ncontrol), Nk * B_nec_alpha[2]);
                      gsl_spmatrix_set(J, ridx, align_idx + CIDX2(1, ALIGN_P, istart_align + k, align_ncontrol), Nk * B_nec_beta[2]);
                      gsl_spmatrix_set(J, ridx, align_idx + CIDX2(2, ALIGN_P, istart_align + k, align_ncontrol), Nk * B_nec_gamma[2]);
                    }

                  if (fit_fluxcal)
                    {
                      size_t idx;

                      for (k = 0; k < fluxcal_spline_p->spline_order; ++k)
                        {
                          double Nk = gsl_vector_get(N_fluxcal, k);

                          for (idx = 0; idx < FLUXCAL_P; ++idx)
                            {
                              gsl_spmatrix_set(J, ridx, fluxcal_idx + CIDX2(idx, FLUXCAL_P, istart_fluxcal + k, fluxcal_ncontrol),
                                               Nk * gsl_matrix_get(&jac_fluxcal_align.matrix, 2, idx));
                            }
                        }
                    }
                }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double F;

              if (fit_align && fit_fluxcal)
                {
                  gsl_vector_view vB_vfm = gsl_vector_view_array(B_vfm, 3);
                  size_t idx;

                  /* F = || B_vfm(c) || */
                  F = gsl_hypot3(B_vfm[0], B_vfm[1], B_vfm[2]);

                  for (idx = 0; idx < FLUXCAL_P; ++idx)
                    {
                      gsl_vector_view dBdc_vfm = gsl_matrix_column(&jac_fluxcal.matrix, idx); /* d/dc_i B_vfm(c) */
                      double result;

                      /* result = B_vfm(c) . d/dc_i B_vfm(c) */
                      gsl_blas_ddot(&vB_vfm.vector, &dBdc_vfm.vector, &result);

                      for (k = 0; k < fluxcal_spline_p->spline_order; ++k)
                        {
                          double Nk = gsl_vector_get(N_fluxcal, k);
                          gsl_spmatrix_set(J, ridx, fluxcal_idx + CIDX2(idx, FLUXCAL_P, istart_fluxcal + k, fluxcal_ncontrol),
                                           Nk * (result / F));
                        }
                    }
                }

              ++ridx;
            }
        }
    }

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_J2: leaving function (%g seconds, nnz = %zu [%.1f%%])\n",
               time_diff(tv0, tv1), gsl_spmatrix_nnz(J),
               (double) gsl_spmatrix_nnz(J) / (double) (J->size1 * J->size2) * 100.0);

#if 0
  printsp_octave(J, "J2");
  exit(1);
#endif

  return s;
}

/*
mfield_calc_df()
  Compute Jacobian matrix J(x) using OpenMP to
calculate Green's functions quickly.

Inputs: x      - parameter vector, length p
        params - mfield workspace
        J      - (output) J(x), n-by-p
*/

static int
mfield_calc_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
  int s = GSL_SUCCESS;
  mfield_workspace *w = (mfield_workspace *) params;
  const size_t gauss_spline_order = gsl_bspline2_order(w->gauss_spline_workspace_p[0]);
  const mfield_parameters *mparams = &(w->params);
  size_t i, j;
  struct timeval tv0, tv1;

  mfield_debug("mfield_calc_df: entering function...\n");
  gettimeofday(&tv0, NULL);

  /* initialize matrix to 0 (only needed on first call) */
  if (w->J2->nz == 0)
    gsl_matrix_set_zero(J);

  if (w->p_sparse > 0)
    {
      /* calculate sparse part of Jacobian J_2 */
      mfield_calc_J2(x, w->J2, w);

      /* copy J_2(x) into the full Jacobian */
      mfield_debug("mfield_calc_df: copying J_2 into J...");

#pragma omp parallel for private(i)
      for (i = 0; i < w->J2->nz; ++i)
        gsl_matrix_set(J, w->J2->i[i], w->J2->p[i] + w->p_int, w->J2->data[i]);

      mfield_debug("done\n");

#if 0
      printsp_octave(w->J2, "J2");
      exit(1);
#endif
    }

  mfield_debug("mfield_calc_df: building dense part of Jacobian...");

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      double * t_year = w->data_workspace_p->t_year[i];
      int fit_bias = w->params.fit_cbias && (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
#if 1 /*XXX*/
          double t = mptr->ts[j];       /* use scaled time */
#endif
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          gsl_vector_view N_gauss = gsl_matrix_column(w->gauss_spline_workspace_p[thread_id]->dB, 0);
          gsl_vector_view dN_gauss = gsl_matrix_column(w->gauss_spline_workspace_p[thread_id]->dB, 1);
          size_t istart_gauss;
          size_t k;

          /* internal Green's functions for current point */
          gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
          gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
          gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);

          /* internal Green's functions for gradient point */
          gsl_vector_view vx_ns = gsl_matrix_row(w->omp_dX_grad, thread_id);
          gsl_vector_view vy_ns = gsl_matrix_row(w->omp_dY_grad, thread_id);
          gsl_vector_view vz_ns = gsl_matrix_row(w->omp_dZ_grad, thread_id);

#if MFIELD_FIT_EXTFIELD
          size_t extidx = 0;
          double extcoeff = 0.0;
          double dB_ext[3];
#endif

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              /* calculate internal Green's functions */
              green_calc_int2(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, w->green_array_p[thread_id]);

              /* calculate non-zero B-splines (and possibly first derivatives) for the Gauss coefficients for this timestamp */
              if (mptr->flags[j] & (MAGDATA_FLG_DXDT | MAGDATA_FLG_DYDT | MAGDATA_FLG_DZDT))
                {
                  gsl_bspline2_eval_deriv_basis_nonzero(t_year[j], 1, w->gauss_spline_workspace_p[thread_id]->dB,
                                                        &istart_gauss, w->gauss_spline_workspace_p[thread_id]);
                }
              else
                {
                  gsl_bspline2_eval_basis_nonzero(t_year[j], &N_gauss.vector, &istart_gauss, w->gauss_spline_workspace_p[thread_id]);
                }
            }

          /* calculate internal Green's functions for gradient point (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
            {
              green_calc_int2(mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j],
                              &vx_ns.vector, &vy_ns.vector, &vz_ns.vector,
                              w->green_array_p[thread_id]);
            }

#if MFIELD_FIT_EXTFIELD
          /* external field is fitted only to data points used for main field modeling */
          if ((MAGDATA_ExistX(mptr->flags[j]) || MAGDATA_ExistY(mptr->flags[j]) ||
               MAGDATA_ExistZ(mptr->flags[j]) || MAGDATA_ExistScalar(mptr->flags[j])) &&
              (MAGDATA_FitMF(mptr->flags[j])))
            {
              extidx = mfield_extidx(mptr->t[j], w);
              extcoeff = gsl_vector_get(x, extidx);

              /* compute external field model */
              mfield_nonlinear_model_ext(mptr->r[j], mptr->theta[j], mptr->phi[j],
                                         x, dB_ext, w);
            }
#endif

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);

                  jacobian_row_int(&N_gauss.vector, istart_gauss, &vx.vector, &Jv.vector, w);

#if MFIELD_FIT_EXTFIELD
                  gsl_matrix_set(J, ridx, extidx, dB_ext[0]);
#endif

                  if (fit_bias)
                    gsl_matrix_set(J, ridx, w->bias_offset + w->bias_idx[i], -1.0);
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);

                  jacobian_row_int(&N_gauss.vector, istart_gauss, &vy.vector, &Jv.vector, w);

#if MFIELD_FIT_EXTFIELD
                  gsl_matrix_set(J, ridx, extidx, dB_ext[1]);
#endif

                  if (fit_bias)
                    gsl_matrix_set(J, ridx, w->bias_offset + w->bias_idx[i] + 1, -1.0);
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);

                  jacobian_row_int(&N_gauss.vector, istart_gauss, &vz.vector, &Jv.vector, w);

#if MFIELD_FIT_EXTFIELD
                  gsl_matrix_set(J, ridx, extidx, dB_ext[2]);
#endif

                  if (fit_bias)
                    gsl_matrix_set(J, ridx, w->bias_offset + w->bias_idx[i] + 2, -1.0);
                }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              gsl_vector_view vf = gsl_matrix_row(w->omp_dF, thread_id);
              gsl_vector_view vf_core = gsl_vector_subvector(&vf.vector, 0, w->nnm_core);
              double X, Y, Z, F;

              /* compute internal X, Y, Z */
              X = mfield_nonlinear_model_int(&N_gauss.vector, istart_gauss, &vx.vector, x, w);
              Y = mfield_nonlinear_model_int(&N_gauss.vector, istart_gauss, &vy.vector, x, w);
              Z = mfield_nonlinear_model_int(&N_gauss.vector, istart_gauss, &vz.vector, x, w);

              /* add apriori (external and crustal) field */
              X += mptr->Bx_model[j];
              Y += mptr->By_model[j];
              Z += mptr->Bz_model[j];

#if MFIELD_FIT_EXTFIELD
              /* add external field correction */
              X += extcoeff * dB_ext[0];
              Y += extcoeff * dB_ext[1];
              Z += extcoeff * dB_ext[2];
#endif

              F = gsl_hypot3(X, Y, Z);

              X /= F;
              Y /= F;
              Z /= F;

              /* compute 1 / F * (X dX + Y dY + Z dZ) */
              for (k = 0; k < w->nnm_tot; ++k)
                {
                  double dXk = gsl_vector_get(&vx.vector, k);
                  double dYk = gsl_vector_get(&vy.vector, k);
                  double dZk = gsl_vector_get(&vz.vector, k);
                  double val = X * dXk + Y * dYk + Z * dZk;

                  gsl_vector_set(&vf.vector, k, val);
                }

              /* compute -N_k(t_i) / F * (X dX + Y dY + Z dZ) */
              for (k = 0; k < gauss_spline_order; ++k)
                {
                  double Nk = gsl_vector_get(&N_gauss.vector, k);
                  gsl_vector_view tmp = gsl_matrix_subrow(J, ridx, (k + istart_gauss) * w->nnm_core, w->nnm_core);
                  gsl_vector_axpby(-Nk, &vf_core.vector, 0.0, &tmp.vector);
                }

              if (w->nnm_crust > 0)
                {
                  gsl_vector_view vf_crust = gsl_vector_subvector(&vf.vector, w->nnm_core, w->nnm_crust);
                  gsl_vector_view tmp = gsl_matrix_subrow(J, ridx, w->p_core, w->nnm_crust);
                  gsl_vector_axpby(-1.0, &vf_crust.vector, 0.0, &tmp.vector);
                }

#if MFIELD_FIT_EXTFIELD
              Jv = gsl_matrix_row(J, ridx);
              gsl_vector_set(&Jv.vector, extidx, (X * dB_ext[0] + Y * dB_ext[1] + Z * dB_ext[2]) / F);
#endif

              ++ridx;
            }

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              if (mptr->flags[j] & MAGDATA_FLG_DXDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p_int);
                  jacobian_row_int_SV(&dN_gauss.vector, istart_gauss, &vx.vector, &Jv.vector, w);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DYDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p_int);
                  jacobian_row_int_SV(&dN_gauss.vector, istart_gauss, &vy.vector, &Jv.vector, w);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DZDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p_int);
                  jacobian_row_int_SV(&dN_gauss.vector, istart_gauss, &vz.vector, &Jv.vector, w);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                  jacobian_row_int_grad(t, mptr->ts_ns[j], &vx.vector, &vx_ns.vector, &Jv.vector, w);
                }

#if OLD_EULER
              /* check if fitting Euler angles to this data point */
              if (fit_align && MAGDATA_FitAlign(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, align_idx, B_nec_alpha[0] - B_nec_alpha_ns[0]);
                  gsl_matrix_set(J, ridx, align_idx + 1, B_nec_beta[0] - B_nec_beta_ns[0]);
                  gsl_matrix_set(J, ridx, align_idx + 2, B_nec_gamma[0] - B_nec_gamma_ns[0]);
                }
#endif

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                  jacobian_row_int_grad(t, mptr->ts_ns[j], &vy.vector, &vy_ns.vector, &Jv.vector, w);
                }

#if OLD_EULER
              /* check if fitting Euler angles to this data point */
              if (fit_align && MAGDATA_FitAlign(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, align_idx, B_nec_alpha[1] - B_nec_alpha_ns[1]);
                  gsl_matrix_set(J, ridx, align_idx + 1, B_nec_beta[1] - B_nec_beta_ns[1]);
                  gsl_matrix_set(J, ridx, align_idx + 2, B_nec_gamma[1] - B_nec_gamma_ns[1]);
                }
#endif

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                  jacobian_row_int_grad(t, mptr->ts_ns[j], &vz.vector, &vz_ns.vector, &Jv.vector, w);
                }

#if OLD_EULER
              /* check if fitting Euler angles to this data point */
              if (fit_align && MAGDATA_FitAlign(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, align_idx, B_nec_alpha[2] - B_nec_alpha_ns[2]);
                  gsl_matrix_set(J, ridx, align_idx + 1, B_nec_beta[2] - B_nec_beta_ns[2]);
                  gsl_matrix_set(J, ridx, align_idx + 2, B_nec_gamma[2] - B_nec_gamma_ns[2]);
                }
#endif

              ++ridx;
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  mfield_debug("done\n");

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_df: leaving function (%g seconds)\n", time_diff(tv0, tv1));

  if (mparams->regularize && !mparams->synth_data)
    {
      /* copy L^T into lower portion of J */
      gsl_matrix_view m = gsl_matrix_submatrix(J, w->nres, 0, w->p, w->p);

      for (i = 0; i < w->L->nz; ++i)
        gsl_matrix_set(&m.matrix, w->L->p[i], w->L->i[i], w->L->data[i]);
    }

#if 0
  print_octave(J, "J");
  exit(1);
#elif 0
  {
    static int niter = 0;

    if (niter++ == 0)
      {
        const size_t p = J->size2;
        gsl_matrix * JTJ = gsl_matrix_alloc(p, p);
        gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, J, 0.0, JTJ);

        printv_octave(x, "xa");
        printsym_octave(JTJ, "JTJa");

        if (w->p_sparse > 0)
          printsp_octave(w->J2, "J2a");

        gsl_matrix_free(JTJ);
        exit(1);
      }
  }
#endif

  return s;
}

/*
mfield_nonlinear_model()
  Compute total model vector for a given residual

Inputs: res_flag  - 0 = normal residual
                    1 = gradient residual
        x         - parameter vector
        mptr      - magdata structure
        src_idx   - index of data source (mptr) in [0,w->nsat-1]
        data_idx  - index of datum in magdata
        thread_id - OpenMP thread id
        B_model   - (output) B_model (X,Y,Z) in NEC
        w         - workspace

Notes:
1) On output, w->omp_d{X,Y,Z}(thread_id,:) is filled with internal Green's functions
for dX/dg, dY/dg, dZ/dg
*/

static int
mfield_nonlinear_model(const int res_flag, const gsl_vector * x, const magdata * mptr, const size_t src_idx,
                       const size_t data_idx, const size_t thread_id, double B_model[3], mfield_workspace *w)
{
  int s = 0;
  const double * t_year = w->data_workspace_p->t_year[src_idx];
  double t, r, theta, phi;
  gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
  gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
  gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);
  gsl_vector *N_gauss = w->gauss_spline_workspace_p[thread_id]->B;
  double B_int[3], B_prior[3], B_extcorr[3], B_bias[3];
  size_t k, istart;

  if (res_flag == 0)
    {
      /* residual is for this data point specified by 'data_idx' */
      t = t_year[data_idx];
      r = mptr->r[data_idx];
      theta = mptr->theta[data_idx];
      phi = mptr->phi[data_idx];

      /* load apriori model of external (and possibly crustal) field */
      B_prior[0] = mptr->Bx_model[data_idx];
      B_prior[1] = mptr->By_model[data_idx];
      B_prior[2] = mptr->Bz_model[data_idx];
    }
  else if (res_flag == 1)
    {
      /* residual is for gradient (N/S or E/W) */
      t = epoch2year(mptr->t_ns[data_idx]);
      r = mptr->r_ns[data_idx];
      theta = mptr->theta_ns[data_idx];
      phi = mptr->phi_ns[data_idx];

      /* load apriori model of external (and possibly crustal) field */
      B_prior[0] = mptr->Bx_model_ns[data_idx];
      B_prior[1] = mptr->By_model_ns[data_idx];
      B_prior[2] = mptr->Bz_model_ns[data_idx];
    }

  green_calc_int2(r, theta, phi, &vx.vector, &vy.vector, &vz.vector,
                  w->green_array_p[thread_id]);

  /* compute internal field model */
  gsl_bspline2_eval_basis_nonzero(t, N_gauss, &istart, w->gauss_spline_workspace_p[thread_id]);
  B_int[0] = mfield_nonlinear_model_int(N_gauss, istart, &vx.vector, x, w);
  B_int[1] = mfield_nonlinear_model_int(N_gauss, istart, &vy.vector, x, w);
  B_int[2] = mfield_nonlinear_model_int(N_gauss, istart, &vz.vector, x, w);

#if MFIELD_FIT_EXTFIELD

  /* external field is fitted only to data points used for main field modeling */
  if ((MAGDATA_ExistX(mptr->flags[data_idx]) || MAGDATA_ExistY(mptr->flags[data_idx]) ||
       MAGDATA_ExistZ(mptr->flags[data_idx]) || MAGDATA_ExistScalar(mptr->flags[data_idx])) &&
      (MAGDATA_FitMF(mptr->flags[data_idx])))
    {
      size_t extidx = mfield_extidx(t, w);
      double extcoeff = gsl_vector_get(x, extidx);
      double dB_ext[3];

      /* compute external field model correction */
      mfield_nonlinear_model_ext(r, theta, phi, x, dB_ext, w);
    }

  /* add correction to POMME field */
  B_extcorr[0] = extcoeff * dB_ext[0];
  B_extcorr[1] = extcoeff * dB_ext[1];
  B_extcorr[2] = extcoeff * dB_ext[2];

#else

  B_extcorr[0] = B_extcorr[1] = B_extcorr[2] = 0.0;

#endif

  if (w->params.fit_cbias && mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
    {
      B_bias[0] = gsl_vector_get(x, w->bias_offset + w->bias_idx[src_idx]);
      B_bias[1] = gsl_vector_get(x, w->bias_offset + w->bias_idx[src_idx] + 1);
      B_bias[2] = gsl_vector_get(x, w->bias_offset + w->bias_idx[src_idx] + 2);
    }
  else
    {
      B_bias[0] = B_bias[1] = B_bias[2] = 0.0;
    }

  /* compute total modeled field (internal + external) */
  for (k = 0; k < 3; ++k)
    B_model[k] = B_int[k] + B_prior[k] + B_extcorr[k] + B_bias[k];

  return s;
}

/*
mfield_nonlinear_model_int()
  Evaluate internal field model for a given coefficient vector

Inputs: N      - vector of non-zero B-splines N_j(t), length gauss_spline_order
        istart - index of first non-zero B-spline in [0, ncontrol-1]
        dB     - vector of basis functions (dX/dg,dY/dg,dZ/dg), length nnm_tot
        g      - model coefficients, length at least p_int
        w      - workspace

Return: model = sum_{nmj} g_{nm,j} N_j(t) dB_{nm}(r)
*/

static double
mfield_nonlinear_model_int(const gsl_vector * N, const size_t istart, const gsl_vector *dB,
                           const gsl_vector *g, const mfield_workspace *w)
{
  gsl_vector_const_view dB_core = gsl_vector_const_subvector(dB, 0, w->nnm_core);
  double val = 0.0;
  double dot;
  size_t j;

  for (j = 0; j < w->params.gauss_spline_order; ++j)
    {
      double Nj = gsl_vector_get(N, j);
      gsl_vector_const_view v = gsl_vector_const_subvector(g, (j + istart) * w->nnm_core, w->nnm_core);

      gsl_blas_ddot(&v.vector, &dB_core.vector, &dot);
      val += Nj * dot;
    }

  if (w->nnm_crust > 0 && dB->size > w->nnm_core)
    {
      gsl_vector_const_view dB_crust = gsl_vector_const_subvector(dB, w->nnm_core, w->nnm_crust);
      gsl_vector_const_view g_crust = gsl_vector_const_subvector(g, w->p_core, w->nnm_crust);

      gsl_blas_ddot(&g_crust.vector, &dB_crust.vector, &dot);
      val += dot;
    }

  return val;
}

/*
mfield_nonlinear_model_SV()
  Compute total time derivative model vector d/dt B for a given residual

Inputs: x          - parameter vector
        mptr       - magdata structure
        idx        - index of datum in magdata
        thread_id  - OpenMP thread id
        dBdt_model - (output) d/dt B_model (X,Y,Z) in NEC
        w          - workspace

Notes:
1) On output, w->omp_d{X,Y,Z}(thread_id,:) is filled with internal Green's functions
for dX/dg, dY/dg, dZ/dg
*/

static int
mfield_nonlinear_model_SV(const gsl_vector * x, const magdata * mptr, const size_t idx,
                          const size_t thread_id, double dBdt_model[3], mfield_workspace *w)
{
  int s = 0;
  double r = mptr->r[idx];
  double theta = mptr->theta[idx];
  double phi = mptr->phi[idx];
  gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
  gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
  gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);
  gsl_vector_view dN_gauss = gsl_matrix_column(w->gauss_spline_workspace_p[thread_id]->dB, 1);
  size_t istart;

  /* residual is for this data point specified by 'idx' */
  green_calc_int2(r, theta, phi, &vx.vector, &vy.vector, &vz.vector,
                  w->green_array_p[thread_id]);

  gsl_bspline2_eval_deriv_basis_nonzero(epoch2year(mptr->t[idx]), 1, w->gauss_spline_workspace_p[thread_id]->dB,
                                        &istart, w->gauss_spline_workspace_p[thread_id]);

  /* compute internal field model */
  dBdt_model[0] = mfield_nonlinear_model_int_SV(&dN_gauss.vector, istart, &vx.vector, x, w);
  dBdt_model[1] = mfield_nonlinear_model_int_SV(&dN_gauss.vector, istart, &vy.vector, x, w);
  dBdt_model[2] = mfield_nonlinear_model_int_SV(&dN_gauss.vector, istart, &vz.vector, x, w);

  return s;
}

/*
mfield_nonlinear_model_int_SV()
  Evaluate time derivative of internal field model d/dt B_{int}(r,t; g) for a given
coefficient vector and a given component

Inputs: dN     - vector of non-zero B-spline derivatives d/dt N_j(t), length gauss_spline_order
        istart - index of first non-zero B-spline in [0, ncontrol-1]
        dB     - vector of basis functions (dX/dg,dY/dg,dZ/dg), length nnm_tot
        g      - model coefficients, length at least p_int
        w      - workspace

Return: model = sum_{nmk} g_{n,k}^m d/dt N_k(t) B_n^m(r)
*/

static double
mfield_nonlinear_model_int_SV(const gsl_vector * dN, const size_t istart, const gsl_vector *dB,
                              const gsl_vector *g, const mfield_workspace *w)
{
  gsl_vector_const_view dB_core = gsl_vector_const_subvector(dB, 0, w->nnm_core);
  double val = 0.0;
  size_t j;

  for (j = 0; j < w->params.gauss_spline_order; ++j)
    {
      double dNj = gsl_vector_get(dN, j);
      gsl_vector_const_view v = gsl_vector_const_subvector(g, (j + istart) * w->nnm_core, w->nnm_core);
      double dot;

      gsl_blas_ddot(&v.vector, &dB_core.vector, &dot);
      val += dNj * dot;
    }

  return val;
}

/*
jacobian_row_int()
  Fill in row of Jacobian corresponding to internal Green's functions

J(row,:) = - [ dB_core, dB_crust ]

Inputs: N      - vector of non-zero B-splines N_j(t), length gauss_spline_order
        istart - index of first non-zero B-spline in [0, ncontrol-1]
        dB     - Green's functions for desired component (X,Y,Z), size nnm_tot
        J      - (output) Jacobian row; columns 1 - w->p_int will be filled in, size p_int
        w      - workspace
*/

static inline int
jacobian_row_int(const gsl_vector * N, const size_t istart, const gsl_vector * dB,
                 gsl_vector * J, mfield_workspace * w)
{
  int s = 0;
  const size_t order = gsl_bspline2_order(w->gauss_spline_workspace_p[0]);
  gsl_vector_const_view dB_core = gsl_vector_const_subvector(dB, 0, w->nnm_core);
  size_t j;

  for (j = 0; j < order; ++j)
    {
      double Nj = gsl_vector_get(N, j);
      gsl_vector_view v = gsl_vector_subvector(J, (j + istart) * w->nnm_core, w->nnm_core);
      gsl_vector_axpby(-Nj, &dB_core.vector, 0.0, &v.vector);
    }

  if (w->nnm_crust > 0)
    {
      gsl_vector_const_view dB_crust = gsl_vector_const_subvector(dB, w->nnm_core, w->nnm_crust);
      gsl_vector_view J_crust = gsl_vector_subvector(J, w->p_core, w->nnm_crust);
      gsl_vector_axpby(-1.0, &dB_crust.vector, 0.0, &J_crust.vector);
    }

  return s;
}

/*
jacobian_row_int_SV()
  Fill in row of Jacobian corresponding to SV internal Green's functions

J(row,:) = [ 0, dB, t dB ]

Inputs: dN     - first derivative of B-spline basis functions, length gauss_spline_order
        istart - index of first non-zero B-spline in [0, ncontrol-1]
        dB     - Green's functions for desired component (X,Y,Z), size nnm_max
        J      - (output) Jacobian row; columns 1 - w->p_int will be filled in, size p_int
        w      - workspace
*/

static inline int
jacobian_row_int_SV(const gsl_vector * dN, const size_t istart, const gsl_vector * dB,
                    gsl_vector * J, mfield_workspace * w)
{
  int s = 0;
  const size_t order = gsl_bspline2_order(w->gauss_spline_workspace_p[0]);
  gsl_vector_const_view dB_core = gsl_vector_const_subvector(dB, 0, w->nnm_core);
  size_t j;

  for (j = 0; j < order; ++j)
    {
      double dNj = gsl_vector_get(dN, j);
      gsl_vector_view v = gsl_vector_subvector(J, (j + istart) * w->nnm_core, w->nnm_core);
      gsl_vector_axpby(-dNj, &dB_core.vector, 0.0, &v.vector);
    }

  return s;
}

/*
jacobian_row_int_grad()
  Fill in row of Jacobian corresponding to internal Green's functions for
a gradient residual

J(row,:) = [ dB2 - dB, t2 * dB2 - t * dB, 1/2 t2^2 dB2 - 1/2 t^2 dB ]

accounting for differences in spherical harmonic degrees for MF, SV, and SA

Inputs: t   - scaled time of measurement
        t2  - scaled time of gradient measurement
        dB  - Green's functions for desired component (X,Y,Z), size nnm_max
        dB2 - Green's functions for desired component (X,Y,Z) of gradient point (N/S or E/W), size nnm_max
        J   - (output) Jacobian row; columns 1 - w->p_int will be filled in, size p_int
        w   - workspace
*/

static inline int
jacobian_row_int_grad(const double t, const double t2, const gsl_vector * dB, const gsl_vector * dB2,
                      gsl_vector * J, mfield_workspace * w)
{
  int s = 0;
  const double fac1 = 0.5 * t * t;
  const double fac2 = 0.5 * t2 * t2;
  size_t i;

  for (i = 0; i < w->nnm_max; ++i)
    {
      double dBi = gsl_vector_get(dB, i);
      double dB2i = gsl_vector_get(dB2, i);

      /* main field portion */
      if (i < w->nnm_mf)
        gsl_vector_set(J, i, dB2i - dBi);

      /* secular variation portion */
      if (i < w->nnm_sv)
        gsl_vector_set(J, w->sv_offset + i, t2 * dB2i - t * dBi);

      if (i < w->nnm_sa)
        gsl_vector_set(J, w->sa_offset + i, fac2 * dB2i - fac1 * dBi);
    }

  return s;
}
