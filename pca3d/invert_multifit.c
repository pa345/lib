/*
 * invert_multifit.c
 *
 * Contains routines for fitting magnetic field module using gsl_multifit_nlinear framework
 */

static int invert_calc_f(const gsl_vector *x, void *params, gsl_vector *f);
static int invert_calc_Wf(const gsl_vector *x, void *params, gsl_vector *f);
static int invert_calc_df(const gsl_vector *x, void *params, gsl_matrix *J);
static int invert_nonlinear_model(const int res_flag, const gsl_vector * x, const magdata * mptr, const size_t src_idx,
                                  const size_t data_idx, const size_t thread_id, double B_model[3], invert_workspace *w);

/*
invert_calc_f()
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
invert_calc_f(const gsl_vector *x, void *params, gsl_vector *f)
{
  int s = GSL_SUCCESS;
  invert_workspace *w = (invert_workspace *) params;
  const invert_parameters *mparams = &(w->params);
  size_t i, j;
  struct timeval tv0, tv1;

  invert_debug("invert_calc_f: entering function...\n");
  gettimeofday(&tv0, NULL);

  gsl_vector_set_zero(f);

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = invert_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          double B_model[3];            /* internal + external */
          double B_obs[3];              /* observation vector NEC frame */
          double dBdt_model[3];         /* dB/dt (SV of internal model) */
          double dBdt_obs[3];           /* SV observation vector (NEC frame) */
          double F_obs;                 /* scalar field measurement */

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* compute vector model for this residual */
          invert_nonlinear_model(0, x, mptr, i, j, thread_id, B_model, w);

          /* compute vector SV model for this residual */
          if (MAGDATA_FitMF(mptr->flags[j]) && (mptr->flags[j] & (MAGDATA_FLG_DXDT | MAGDATA_FLG_DYDT | MAGDATA_FLG_DZDT)))
            {
              dBdt_obs[0] = mptr->dXdt_nec[j];
              dBdt_obs[1] = mptr->dYdt_nec[j];
              dBdt_obs[2] = mptr->dZdt_nec[j];
            }

          /* use supplied NEC vector */
          B_obs[0] = mptr->Bx_nec[j];
          B_obs[1] = mptr->By_nec[j];
          B_obs[2] = mptr->Bz_nec[j];

          F_obs = mptr->F[j];

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

  if (mparams->regularize && !mparams->synth_data)
    {
      /* store L^T*x in bottom of f for regularization */
      gsl_vector_view v = gsl_vector_subvector(f, w->nres, w->p);
      gsl_spblas_dusmv(CblasTrans, 1.0, w->L, x, 0.0, &v.vector);
    }

  gettimeofday(&tv1, NULL);
  invert_debug("invert_calc_f: leaving function (%g seconds, ||f|| = %g)\n",
               time_diff(tv0, tv1), gsl_blas_dnrm2(f));

#if 0
  printv_octave(x, "x");
  printv_octave(f, "f");
  exit(1);
#endif

  return s;
}

/*
invert_calc_Wf()
  Compute weighted residuals:

f~(x) = sqrt(W) f(x)

Inputs: x      - model parameters
        params - parameters
        f      - (output) f~(x)
*/

static int
invert_calc_Wf(const gsl_vector *x, void *params, gsl_vector *f)
{
  int s;
  invert_workspace *w = (invert_workspace *) params;
  size_t i;
  struct timeval tv0, tv1;

  invert_debug("invert_calc_Wf: entering function...\n");
  gettimeofday(&tv0, NULL);

  s = invert_calc_f(x, params, f);
  if (s)
    return s;

  for (i = 0; i < w->nres; ++i)
    {
      double swi = gsl_vector_get(w->sqrt_wts_final, i);
      double *fi = gsl_vector_ptr(f, i);

      *fi *= swi;
    }

  gettimeofday(&tv1, NULL);
  invert_debug("invert_calc_Wf: leaving function (%g seconds, ||sqrt(W) f|| = %g)\n",
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
invert_calc_fvv()
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
invert_calc_fvv(const gsl_vector *x, const gsl_vector * v, void *params, gsl_vector *fvv)
{
  int s = GSL_SUCCESS;
  (void) x;
  (void) v;
  (void) params;
  (void) fvv;

  return s;
}

/*
invert_calc_Wfvv()
  Compute weighted residuals:

fvv~(x) = sqrt(W) fvv(x)

Inputs: x      - model parameters
        v      - velocity vector
        params - parameters
        fvv    - (output) fvv~(x)
*/

static int
invert_calc_Wfvv(const gsl_vector *x, const gsl_vector *v, void *params, gsl_vector *fvv)
{
  int s;
  invert_workspace *w = (invert_workspace *) params;
  size_t i;
  struct timeval tv0, tv1;

  invert_debug("invert_calc_Wfvv: entering function...\n");
  gettimeofday(&tv0, NULL);

  s = invert_calc_fvv(x, v, params, fvv);
  if (s)
    return s;

  for (i = 0; i < w->nres; ++i)
    {
      double wi = gsl_vector_get(w->wts_final, i);
      double *fi = gsl_vector_ptr(fvv, i);

      *fi *= sqrt(wi);
    }

  gettimeofday(&tv1, NULL);
  invert_debug("invert_calc_Wfvv: leaving function (%g seconds)\n", time_diff(tv0, tv1));

  return GSL_SUCCESS;
}

/*
invert_calc_df()
  Compute Jacobian matrix J(x) using OpenMP to
calculate Green's functions quickly.

Inputs: x      - parameter vector, length p
        params - invert workspace
        J      - (output) J(x), n-by-p
*/

static int
invert_calc_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
  int s = GSL_SUCCESS;
  invert_workspace *w = (invert_workspace *) params;
  const invert_parameters *mparams = &(w->params);
  size_t i, j;
  struct timeval tv0, tv1;

  invert_debug("invert_calc_df: entering function...\n");
  gettimeofday(&tv0, NULL);

  invert_debug("invert_calc_df: building dense part of Jacobian...");

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = invert_data_ptr(i, w->data_workspace_p);

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
          size_t k;

          /* internal Green's functions for current point */
          gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
          gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
          gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);

          /* internal Green's functions for gradient point */
          gsl_vector_view vx_ns = gsl_matrix_row(w->omp_dX_grad, thread_id);
          gsl_vector_view vy_ns = gsl_matrix_row(w->omp_dY_grad, thread_id);
          gsl_vector_view vz_ns = gsl_matrix_row(w->omp_dZ_grad, thread_id);

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              ++ridx;
            }

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              if (mptr->flags[j] & MAGDATA_FLG_DXDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p_int);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DYDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p_int);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DZDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p_int);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                }

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                }

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p_int);
                }

              ++ridx;
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  invert_debug("done\n");

  gettimeofday(&tv1, NULL);
  invert_debug("invert_calc_df: leaving function (%g seconds)\n", time_diff(tv0, tv1));

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
#elif 1
  {
    static int niter = 0;

    if (niter++ == 0)
      {
        const size_t p = J->size2;
        gsl_matrix * JTJ = gsl_matrix_alloc(p, p);
        gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, J, 0.0, JTJ);

        printv_octave(x, "xa");
        printsym_octave(JTJ, "JTJa");

        gsl_matrix_free(JTJ);
        exit(1);
      }
  }
#endif

  return s;
}

/*
invert_nonlinear_model()
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
invert_nonlinear_model(const int res_flag, const gsl_vector * x, const magdata * mptr, const size_t src_idx,
                       const size_t data_idx, const size_t thread_id, double B_model[3], invert_workspace *w)
{
  int s = 0;

  return s;
}
