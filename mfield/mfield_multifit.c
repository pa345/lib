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
static double mfield_nonlinear_model_int(const double t, const gsl_vector *v,
                                         const gsl_vector *g, const mfield_workspace *w);
static int mfield_nonlinear_model_SV(const gsl_vector * x, const magdata * mptr, const size_t idx,
                                     const size_t thread_id, double dBdt_model[3], mfield_workspace *w);
static double mfield_nonlinear_model_int_SV(const double t, const gsl_vector *v,
                                            const gsl_vector *g, const mfield_workspace *w);
static inline int jacobian_row_int(const double t, const gsl_vector * dB, gsl_vector * J, mfield_workspace * w);
static inline int jacobian_row_int_SV(const double t, const gsl_vector * dB, gsl_vector * J, mfield_workspace * w);
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

  mfield_debug("mfield_calc_f: entering function...\n");
  gettimeofday(&tv0, NULL);

  gsl_vector_set_zero(f);

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      /* Euler angles */
      int fit_euler = mparams->fit_euler && (mptr->global_flags & MAGDATA_GLOBFLG_EULER);
      size_t euler_ncontrol = fit_euler ? gsl_bspline2_ncontrol(w->euler_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)]) : 0;
      size_t euler_idx = w->euler_offset + w->offset_euler[i];

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
          gsl_bspline2_workspace *euler_spline_p = fit_euler ? w->euler_spline_workspace_p[CIDX2(i, w->nsat, thread_id, w->max_threads)] : NULL;
          gsl_bspline2_workspace *fluxcal_spline_p = fit_fluxcal ? w->fluxcal_spline_workspace_p[CIDX2(i, w->nsat, thread_id, w->max_threads)] : NULL;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute vector model for this residual */
          mfield_nonlinear_model(0, x, mptr, i, j, thread_id, B_model, w);

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

          if (fit_euler)
            {
              /*
               * get the Euler angles for this satellite and time period,
               * and apply rotation
               */
              double euler_data[EULER_P];
              gsl_vector_view euler_params = gsl_vector_view_array(euler_data, EULER_P);

              /* Euler angle control points for this dataset */
              gsl_vector_const_view v1 = gsl_vector_const_subvector(x, euler_idx, EULER_P * euler_ncontrol);
              gsl_matrix_const_view euler_control_pts = gsl_matrix_const_view_vector(&v1.vector, EULER_P, euler_ncontrol);

              double *q = &(mptr->q[4*j]);
              double B_vfm[4];

              /* compute Euler angles for this timestamp */
              gsl_bspline2_vector_eval(mptr->t[j], &euler_control_pts.matrix, &euler_params.vector, euler_spline_p);

              B_vfm[0] = mptr->Bx_vfm[j];
              B_vfm[1] = mptr->By_vfm[j];
              B_vfm[2] = mptr->Bz_vfm[j];

              if (fit_fluxcal)
                {
                  double cal_data[FLUXCAL_P];
                  gsl_vector_view cal_params = gsl_vector_view_array(cal_data, FLUXCAL_P);

                  /* fluxgate calibration control points for this dataset */
                  gsl_vector_const_view v2 = gsl_vector_const_subvector(x, fluxcal_idx, FLUXCAL_P * fluxcal_ncontrol);
                  gsl_matrix_const_view fluxcal_control_pts = gsl_matrix_const_view_vector(&v2.vector, FLUXCAL_P, fluxcal_ncontrol);

                  gsl_bspline2_vector_eval(mptr->t[j], &fluxcal_control_pts.matrix, &cal_params.vector, fluxcal_spline_p);
                  mfield_fluxcal_apply_datum(&cal_params.vector, B_vfm, B_vfm);
                }

              /* rotate VFM vector to NEC, B_obs = R_q R_3(alpha) B_VFM(c) */
              euler_vfm2nec(mptr->euler_flags, euler_data[0], euler_data[1], euler_data[2], q, B_vfm, B_obs);

              if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                    MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
                {
                  double *q_ns = &(mptr->q_ns[4*j]);
                  double B_vfm_ns[3];

                  B_vfm_ns[0] = mptr->Bx_vfm_ns[j];
                  B_vfm_ns[1] = mptr->By_vfm_ns[j];
                  B_vfm_ns[2] = mptr->Bz_vfm_ns[j];

                  euler_vfm2nec(mptr->euler_flags, euler_data[0], euler_data[1], euler_data[2], q_ns, B_vfm_ns, B_obs_ns);
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

          if (fit_euler && fit_fluxcal && (mptr->flags[j] & (MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z)))
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
      /* store L*x in bottom of f for regularization */
      gsl_vector_view v = gsl_vector_subvector(f, w->nres, w->p);

      gsl_vector_memcpy(&v.vector, x);
      gsl_vector_mul(&v.vector, w->lambda_diag);
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
  mfield_debug("mfield_calc_Wf: leaving function (%g seconds)\n", time_diff(tv0, tv1));

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
#if 0
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
      int fit_euler = mparams->fit_euler && (mptr->global_flags & MAGDATA_GLOBFLG_EULER);

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

      /* Euler angles */
      int fit_euler = w->params.fit_euler && (mptr->global_flags & MAGDATA_GLOBFLG_EULER);
      size_t euler_ncontrol = fit_euler ? gsl_bspline2_ncontrol(w->euler_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)]) : 0;
      size_t euler_idx = w->euler_offset - w->p_int + w->offset_euler[i];

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
          gsl_bspline2_workspace *euler_spline_p = fit_euler ? w->euler_spline_workspace_p[CIDX2(i, w->nsat, thread_id, w->max_threads)] : NULL;
          double euler_data[EULER_P];
          gsl_vector_view euler_params = gsl_vector_view_array(euler_data, EULER_P);
          gsl_vector *N_euler = fit_euler ? euler_spline_p->B : NULL;
          size_t istart_euler;
          double B_vfm[3], B_nec_alpha[3], B_nec_beta[3], B_nec_gamma[3];
          double B_nec_alpha_ns[3], B_nec_beta_ns[3], B_nec_gamma_ns[3];

          /* fluxgate calibration parameters */
          gsl_bspline2_workspace *fluxcal_spline_p = fit_fluxcal ? w->fluxcal_spline_workspace_p[CIDX2(i, w->nsat, thread_id, w->max_threads)] : NULL;
          double cal_data[FLUXCAL_P];
          gsl_vector_view cal_params = gsl_vector_view_array(cal_data, FLUXCAL_P);
          double jac_fluxcal_data[3 * FLUXCAL_P];
          gsl_matrix_view jac_fluxcal = gsl_matrix_view_array(jac_fluxcal_data, 3, FLUXCAL_P);
          double jac_fluxcal_euler_data[3 * FLUXCAL_P];
          gsl_matrix_view jac_fluxcal_euler = gsl_matrix_view_array(jac_fluxcal_euler_data, 3, FLUXCAL_P);
          gsl_vector *N_fluxcal = fit_fluxcal ? fluxcal_spline_p->B : NULL;
          size_t istart_fluxcal;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute Euler angle derivatives of residual vector */
          if (fit_euler)
            {
              double *q = &(mptr->q[4*j]);
              gsl_vector_const_view v1 = gsl_vector_const_subvector(&y.vector, euler_idx, EULER_P * euler_ncontrol);
              gsl_matrix_const_view euler_control_pts = gsl_matrix_const_view_vector(&v1.vector, EULER_P, euler_ncontrol);

              /* compute Euler angles for this timestamp */
              gsl_bspline2_vector_eval(mptr->t[j], &euler_control_pts.matrix, &euler_params.vector, euler_spline_p);

              /* evaluate non-zero basis splines for time t */
              gsl_bspline2_eval_basis_nonzero(mptr->t[j], N_euler, &istart_euler, euler_spline_p);

              /* get vector in VFM frame */
              B_vfm[0] = mptr->Bx_vfm[j];
              B_vfm[1] = mptr->By_vfm[j];
              B_vfm[2] = mptr->Bz_vfm[j];

              if (fit_fluxcal)
                {
                  /* fluxgate calibration control points for this dataset */
                  gsl_vector_const_view v2 = gsl_vector_const_subvector(&y.vector, fluxcal_idx, FLUXCAL_P * fluxcal_ncontrol);
                  gsl_matrix_const_view fluxcal_control_pts = gsl_matrix_const_view_vector(&v2.vector, FLUXCAL_P, fluxcal_ncontrol);

                  gsl_bspline2_vector_eval(mptr->t[j], &fluxcal_control_pts.matrix, &cal_params.vector, fluxcal_spline_p);

                  /* compute jac_fluxcal := d/dm B_vfm(m) */
                  mfield_fluxcal_jac(&cal_params.vector, B_vfm, &jac_fluxcal.matrix);

                 /* apply R_q * R_3(alpha) to d/dm B_vfm(m) vectors */
                  euler_matrix_vfm2nec(mptr->euler_flags, euler_data[0], euler_data[1], euler_data[2], q, &jac_fluxcal.matrix, &jac_fluxcal_euler.matrix);

                  /* compute B_vfm := P^{-1} S (B_vfm - O) */
                  mfield_fluxcal_apply_datum(&cal_params.vector, B_vfm, B_vfm);

                  /* evaluate non-zero basis splines for time t */
                  gsl_bspline2_eval_basis_nonzero(mptr->t[j], N_fluxcal, &istart_fluxcal, fluxcal_spline_p);
                }

              /* compute alpha derivative of: R_q R_3 B_vfm */
              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_ALPHA, euler_data[0], euler_data[1], euler_data[2], q, B_vfm, B_nec_alpha);

              /* compute beta derivative of: R_q R_3 B_vfm */
              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_BETA, euler_data[0], euler_data[1], euler_data[2], q, B_vfm, B_nec_beta);

              /* compute gamma derivative of: R_q R_3 B_vfm */
              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_GAMMA, euler_data[0], euler_data[1], euler_data[2], q, B_vfm, B_nec_gamma);

              if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                    MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
                {
                  q = &(mptr->q_ns[4*j]);

                  /* get vector in VFM frame */
                  B_vfm[0] = mptr->Bx_vfm_ns[j];
                  B_vfm[1] = mptr->By_vfm_ns[j];
                  B_vfm[2] = mptr->Bz_vfm_ns[j];

                  /* compute alpha derivative of: R_q R_3 B_vfm_ns */
                  euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_ALPHA, euler_data[0], euler_data[1], euler_data[2], q, B_vfm, B_nec_alpha_ns);

                  /* compute beta derivative of: R_q R_3 B_vfm_ns */
                  euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_BETA, euler_data[0], euler_data[1], euler_data[2], q, B_vfm, B_nec_beta_ns);

                  /* compute gamma derivative of: R_q R_3 B_vfm_ns */
                  euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_GAMMA, euler_data[0], euler_data[1], euler_data[2], q, B_vfm, B_nec_gamma_ns);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              /* check if fitting Euler angles to this data point */
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  for (k = 0; k < euler_spline_p->spline_order; ++k)
                    {
                      double Nk = gsl_vector_get(N_euler, k);

                      gsl_spmatrix_set(J, ridx, euler_idx + CIDX2(0, EULER_P, istart_euler + k, euler_ncontrol), Nk * B_nec_alpha[0]);
                      gsl_spmatrix_set(J, ridx, euler_idx + CIDX2(1, EULER_P, istart_euler + k, euler_ncontrol), Nk * B_nec_beta[0]);
                      gsl_spmatrix_set(J, ridx, euler_idx + CIDX2(2, EULER_P, istart_euler + k, euler_ncontrol), Nk * B_nec_gamma[0]);
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
                                               Nk * gsl_matrix_get(&jac_fluxcal_euler.matrix, 0, idx));
                            }
                        }
                    }
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              /* check if fitting Euler angles to this data point */
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  for (k = 0; k < euler_spline_p->spline_order; ++k)
                    {
                      double Nk = gsl_vector_get(N_euler, k);

                      gsl_spmatrix_set(J, ridx, euler_idx + CIDX2(0, EULER_P, istart_euler + k, euler_ncontrol), Nk * B_nec_alpha[1]);
                      gsl_spmatrix_set(J, ridx, euler_idx + CIDX2(1, EULER_P, istart_euler + k, euler_ncontrol), Nk * B_nec_beta[1]);
                      gsl_spmatrix_set(J, ridx, euler_idx + CIDX2(2, EULER_P, istart_euler + k, euler_ncontrol), Nk * B_nec_gamma[1]);
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
                                               Nk * gsl_matrix_get(&jac_fluxcal_euler.matrix, 1, idx));
                            }
                        }
                    }
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              /* check if fitting Euler angles to this data point */
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  for (k = 0; k < euler_spline_p->spline_order; ++k)
                    {
                      double Nk = gsl_vector_get(N_euler, k);

                      gsl_spmatrix_set(J, ridx, euler_idx + CIDX2(0, EULER_P, istart_euler + k, euler_ncontrol), Nk * B_nec_alpha[2]);
                      gsl_spmatrix_set(J, ridx, euler_idx + CIDX2(1, EULER_P, istart_euler + k, euler_ncontrol), Nk * B_nec_beta[2]);
                      gsl_spmatrix_set(J, ridx, euler_idx + CIDX2(2, EULER_P, istart_euler + k, euler_ncontrol), Nk * B_nec_gamma[2]);
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
                                               Nk * gsl_matrix_get(&jac_fluxcal_euler.matrix, 2, idx));
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

              if (fit_euler && fit_fluxcal)
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
    }

  mfield_debug("mfield_calc_df: building dense part of Jacobian...");

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      int fit_bias = w->params.fit_cbias && (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          double t = mptr->ts[j];       /* use scaled time */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */
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

          /* calculate internal Green's functions */
          green_calc_int2(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, w->green_array_p[thread_id]);

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

                  jacobian_row_int(t, &vx.vector, &Jv.vector, w);

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

                  jacobian_row_int(t, &vy.vector, &Jv.vector, w);

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

                  jacobian_row_int(t, &vz.vector, &Jv.vector, w);

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
              gsl_vector_view Jv;
              double X, Y, Z, F;

              /* compute internal X, Y, Z */
              X = mfield_nonlinear_model_int(t, &vx.vector, x, w);
              Y = mfield_nonlinear_model_int(t, &vy.vector, x, w);
              Z = mfield_nonlinear_model_int(t, &vz.vector, x, w);

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

              /* compute (X dX + Y dY + Z dZ) */
              Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
              for (k = 0; k < w->nnm_max; ++k)
                {
                  double dXk = gsl_vector_get(&vx.vector, k);
                  double dYk = gsl_vector_get(&vy.vector, k);
                  double dZk = gsl_vector_get(&vz.vector, k);
                  double val = X * dXk + Y * dYk + Z * dZk;

                  mfield_set_mf(&Jv.vector, k, -val, w);
                  mfield_set_sv(&Jv.vector, k, -t * val, w);
                  mfield_set_sa(&Jv.vector, k, -0.5 * t * t * val, w);
                }

              /* scale by 1/F */
              gsl_vector_scale(&Jv.vector, 1.0 / F);

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
                  jacobian_row_int_SV(t, &vx.vector, &Jv.vector, w);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DYDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p_int);
                  jacobian_row_int_SV(t, &vy.vector, &Jv.vector, w);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DZDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p_int);
                  jacobian_row_int_SV(t, &vz.vector, &Jv.vector, w);
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
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, B_nec_alpha[0] - B_nec_alpha_ns[0]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, B_nec_beta[0] - B_nec_beta_ns[0]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, B_nec_gamma[0] - B_nec_gamma_ns[0]);
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
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, B_nec_alpha[1] - B_nec_alpha_ns[1]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, B_nec_beta[1] - B_nec_beta_ns[1]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, B_nec_gamma[1] - B_nec_gamma_ns[1]);
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
              if (fit_euler && MAGDATA_FitEuler(mptr->flags[j]))
                {
                  gsl_matrix_set(J, ridx, euler_idx, B_nec_alpha[2] - B_nec_alpha_ns[2]);
                  gsl_matrix_set(J, ridx, euler_idx + 1, B_nec_beta[2] - B_nec_beta_ns[2]);
                  gsl_matrix_set(J, ridx, euler_idx + 2, B_nec_gamma[2] - B_nec_gamma_ns[2]);
                }
#endif

              ++ridx;
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  mfield_debug("done\n");

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_df: leaving function (%g seconds)\n", time_diff(tv0, tv1));

#if 0
  print_octave(J, "J");
  exit(1);
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
  double t, ts, r, theta, phi;
  gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
  gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
  gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);
  double B_int[3], B_prior[3], B_extcorr[3], B_bias[3];
  size_t k;

  if (res_flag == 0)
    {
      /* residual is for this data point specified by 'data_idx' */
      t = mptr->t[data_idx];
      ts = mptr->ts[data_idx];
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
      t = mptr->t_ns[data_idx];
      ts = mptr->ts_ns[data_idx];
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
  B_int[0] = mfield_nonlinear_model_int(ts, &vx.vector, x, w);
  B_int[1] = mfield_nonlinear_model_int(ts, &vy.vector, x, w);
  B_int[2] = mfield_nonlinear_model_int(ts, &vz.vector, x, w);

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

Inputs: t - scaled time
        v - vector of basis functions (dX/dg,dY/dg,dZ/dg)
        g - model coefficients
        w - workspace

Return: model = v . g_mf + t*(v . g_sv) + 1/2*t^2*(v . g_sa)
*/

static double
mfield_nonlinear_model_int(const double t, const gsl_vector *v,
                           const gsl_vector *g, const mfield_workspace *w)
{
  double mf = 0.0, sv = 0.0, sa = 0.0, val;

  if (w->nnm_mf > 0)
    {
      /* compute v . x_mf */
      gsl_vector_const_view gmf = gsl_vector_const_subvector(g, 0, w->nnm_mf);
      gsl_vector_const_view vmf = gsl_vector_const_subvector(v, 0, w->nnm_mf);
      gsl_blas_ddot(&vmf.vector, &gmf.vector, &mf);
    }

  if (w->nnm_sv > 0)
    {
      /* compute v . x_sv */
      gsl_vector_const_view gsv = gsl_vector_const_subvector(g, w->sv_offset, w->nnm_sv);
      gsl_vector_const_view vsv = gsl_vector_const_subvector(v, 0, w->nnm_sv);
      gsl_blas_ddot(&vsv.vector, &gsv.vector, &sv);
    }

  if (w->nnm_sa > 0)
    {
      /* compute v . x_sa */
      gsl_vector_const_view gsa = gsl_vector_const_subvector(g, w->sa_offset, w->nnm_sa);
      gsl_vector_const_view vsa = gsl_vector_const_subvector(v, 0, w->nnm_sa);
      gsl_blas_ddot(&vsa.vector, &gsa.vector, &sa);
    }

  val = mf + t * sv + 0.5 * t * t * sa;

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
  double ts, r, theta, phi;
  gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
  gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
  gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);
  double dBdt_int[3];
  size_t k;

  /* residual is for this data point specified by 'idx' */
  ts = mptr->ts[idx];
  r = mptr->r[idx];
  theta = mptr->theta[idx];
  phi = mptr->phi[idx];

  green_calc_int2(r, theta, phi, &vx.vector, &vy.vector, &vz.vector,
                  w->green_array_p[thread_id]);

  /* compute internal field model */
  dBdt_int[0] = mfield_nonlinear_model_int_SV(ts, &vx.vector, x, w);
  dBdt_int[1] = mfield_nonlinear_model_int_SV(ts, &vy.vector, x, w);
  dBdt_int[2] = mfield_nonlinear_model_int_SV(ts, &vz.vector, x, w);

  /* compute total time derivative of internal modeled field */
  for (k = 0; k < 3; ++k)
    dBdt_model[k] = dBdt_int[k];

  return s;
}

/*
mfield_nonlinear_model_int_SV()
  Evaluate time derivative of internal field model d/dt B(r,t; g) for a given coefficient vector

Inputs: t - scaled time
        v - vector of basis functions (dX/dg,dY/dg,dZ/dg)
        g - model coefficients
        w - workspace

Return: model = v . g_sv + t*(v . g_sa)
*/

static double
mfield_nonlinear_model_int_SV(const double t, const gsl_vector *v,
                              const gsl_vector *g, const mfield_workspace *w)
{
  double sv = 0.0, sa = 0.0, val;

  if (w->nnm_sv > 0)
    {
      /* compute v . x_sv */
      gsl_vector_const_view gsv = gsl_vector_const_subvector(g, w->sv_offset, w->nnm_sv);
      gsl_vector_const_view vsv = gsl_vector_const_subvector(v, 0, w->nnm_sv);
      gsl_blas_ddot(&vsv.vector, &gsv.vector, &sv);
    }

  if (w->nnm_sa > 0)
    {
      /* compute v . x_sa */
      gsl_vector_const_view gsa = gsl_vector_const_subvector(g, w->sa_offset, w->nnm_sa);
      gsl_vector_const_view vsa = gsl_vector_const_subvector(v, 0, w->nnm_sa);
      gsl_blas_ddot(&vsa.vector, &gsa.vector, &sa);
    }

  val = sv + t * sa;

  return val;
}

/*
jacobian_row_int()
  Fill in row of Jacobian corresponding to internal Green's functions

J(row,:) = - [ dB, t * dB, 1/2 t^2 dB ]

accounting for differences in spherical harmonic degrees for MF, SV, and SA

Inputs: t  - scaled time of measurement
        dB - Green's functions for desired component (X,Y,Z), size nnm_max
        J  - (output) Jacobian row; columns 1 - w->p_int will be filled in, size p_int
        w  - workspace
*/

static inline int
jacobian_row_int(const double t, const gsl_vector * dB, gsl_vector * J, mfield_workspace * w)
{
  int s = 0;
  const double fac_sa = 0.5 * t * t;
  size_t i;

  for (i = 0; i < w->nnm_max; ++i)
    {
      double dBi = gsl_vector_get(dB, i);

      /* main field portion */
      if (i < w->nnm_mf)
        gsl_vector_set(J, i, -dBi);

      /* secular variation portion */
      if (i < w->nnm_sv)
        gsl_vector_set(J, w->sv_offset + i, -t * dBi);

      /* secular acceleration portion */
      if (i < w->nnm_sa)
        gsl_vector_set(J, w->sa_offset + i, -fac_sa * dBi);
    }

  return s;
}

/*
jacobian_row_int_SV()
  Fill in row of Jacobian corresponding to SV internal Green's functions

J(row,:) = [ 0, dB, t dB ]

accounting for differences in spherical harmonic degrees for MF, SV, and SA

Inputs: t  - scaled time of measurement
        dB - Green's functions for desired component (X,Y,Z), size nnm_max
        J  - (output) Jacobian row; columns 1 - w->p_int will be filled in, size p_int
        w  - workspace
*/

static inline int
jacobian_row_int_SV(const double t, const gsl_vector * dB, gsl_vector * J, mfield_workspace * w)
{
  int s = 0;
  size_t i;

  for (i = 0; i < w->nnm_max; ++i)
    {
      double dBi = gsl_vector_get(dB, i);

      /* main field portion */
      if (i < w->nnm_mf)
        gsl_vector_set(J, i, 0.0);

      /* secular variation portion */
      if (i < w->nnm_sv)
        gsl_vector_set(J, w->sv_offset + i, -dBi);

      /* secular acceleration portion */
      if (i < w->nnm_sa)
        gsl_vector_set(J, w->sa_offset + i, -t * dBi);
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
