/* least squares solver to use,
 *
 * FDF_SOLVER = 0   GSL multifit solver / LM
 * FDF_SOLVER = 1   GSL multilarge solver / LM
 * FDF_SOLVER = 2   GSL multilarge solver / Gauss-Newton
 */
#define FDF_SOLVER     0

#include <gsl/gsl_integration.h>
#include <spatwt/spatwt.h>

static int mfield_init_nonlinear(mfield_workspace *w);
static int mfield_calc_df2(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                           void *params, gsl_vector * v, gsl_matrix * JTJ);
static inline int mfield_jacobian_JTu(const double t, const size_t flags, const double weight,
                                      const gsl_vector * u, const size_t ridx, const gsl_vector * dB_int, const size_t extidx,
                                      const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                                      const double B_nec_beta, const double B_nec_gamma,
                                      gsl_vector *JTy, const mfield_workspace *w);
static inline int mfield_jacobian_SV_JTu(const double t, const size_t flags, const double weight,
                                         const gsl_vector * u, const size_t ridx, const gsl_vector * dB_int,
                                         gsl_vector *JTu, const mfield_workspace *w);
static inline int mfield_jacobian_grad_JTu(const double t, const double t_grad, const size_t flags, const double weight,
                                           const gsl_vector * u, const size_t ridx, const gsl_vector * dB_int,
                                           const gsl_vector * dB_int_grad, gsl_vector *JTu, const mfield_workspace *w);
static inline int mfield_jacobian_Ju(const double t, const size_t flags, const double weight,
                                     const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, const size_t extidx,
                                     const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                                     const double B_nec_beta, const double B_nec_gamma,
                                     gsl_vector *Ju, const mfield_workspace *w);
static inline int mfield_jacobian_JTJ(const double t, const size_t flags, const double weight,
                                      gsl_vector * dB_int, const size_t extidx,
                                      const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                                      const double B_nec_beta, const double B_nec_gamma,
                                      gsl_matrix *JTJ, const mfield_workspace *w);
static inline int mfield_jacobian_row_F(CBLAS_TRANSPOSE_t TransJ, const double t, const double weight,
                                        const gsl_vector * u, const size_t ridx,
                                        gsl_vector * dX, gsl_vector * dY, gsl_vector * dZ,
                                        const double B_model[4], const size_t extidx, const double dB_ext[3],
                                        gsl_vector *J_int, gsl_matrix *JTJ, gsl_vector *v,
                                        const mfield_workspace *w);
static int mfield_vector_green_SV(const double t, const double weight, const gsl_vector *g,
                                  gsl_vector *G, mfield_workspace *w);
static int mfield_vector_green_grad(const double t, const double t_grad, const double weight, const gsl_vector *g,
                                    const gsl_vector *g_grad, gsl_vector *G, mfield_workspace *w);
static int mfield_nonlinear_model_ext(const double r, const double theta,
                                      const double phi, const gsl_vector *g,
                                      double dB[3], const mfield_workspace *w);
static int mfield_nonlinear_regularize_init(mfield_workspace *w);
static int mfield_nonlinear_regularize_core(const double lambda, const size_t nderiv, mfield_workspace *w);
static void mfield_nonlinear_callback(const size_t iter, void *params,
                                      const gsl_multifit_nlinear_workspace *multifit_p);
static void mfield_nonlinear_callback2(const size_t iter, void *params,
                                       const gsl_multilarge_nlinear_workspace *multifit_p);
static int mfield_robust_weights(const gsl_vector * f, gsl_vector * wts, mfield_workspace * w);
static double huber(const double x);
static double bisquare(const double x);
static int gsl_linalg_symband_unpack(const gsl_matrix * AB, gsl_matrix * A);

#include "mfield_multifit.c"
#include "mfield_multilarge.c"
#include "mfield_gn.c"


/*
mfield_calc_nonlinear()
  Solve linear least squares system, using previously stored
satellite data

Inputs: c    - (input/output)
               on input, initial guess for coefficient vector
               on output, final coefficients
               units of nT, nT/year, nT/year^2
        w    - workspace

Notes:
1) On input, w->data_workspace_p must be filled with satellite data
2a) mfield_init() must be called first to initialize various parameters,
    including weights
2b) this includes mfield_init_nonlinear(), called from mfield_init()
3) on output, coefficients are stored in w->c with following units:
4) static coefficients have units of nT
5) SV coefficients have units of nT/dimensionless_time
6) SA coefficients have units of nT/dimensionless_time^2

Return: GSL_SUCCESS if converged, GSL_CONTINUE if not
*/

int
mfield_calc_nonlinear(gsl_vector *c, mfield_workspace *w)
{
  int s = 0;
  const mfield_parameters *params = &(w->params);
  const size_t max_iter = 50;     /* maximum iterations */
  const double xtol = 1.0e-5;
  const double gtol = 1.0e-6;
  const double ftol = 1.0e-6;
  int info;
  const size_t p = w->p;          /* number of coefficients */
  const size_t n = w->nres_tot;   /* number of residuals */
  gsl_multifit_nlinear_fdf fdf;
  gsl_vector *f;
  struct timeval tv0, tv1;
  double res0;                    /* initial residual */

  fdf.f = mfield_calc_f;
  fdf.df = mfield_calc_df;
  fdf.fvv = NULL;
  fdf.n = n;
  fdf.p = p;
  fdf.params = w;

  printv_octave(c, "c0");

  /* compute robust weights with coefficients from previous iteration */
  if (w->niter > 0)
    {
      gsl_vector_view wv = gsl_vector_subvector(w->wts_final, 0, w->nres);

      /* compute residuals f = Y_model - y_data with previous coefficients */
      mfield_calc_f(c, w, w->fvec);

      /* compute robust weights */
      fprintf(stderr, "mfield_calc_nonlinear: computing robust weights...");
      mfield_robust_weights(w->fvec, w->wts_robust, w);
      fprintf(stderr, "done\n");

      /* compute final weights = wts_robust .* wts_spatial */
      gsl_vector_memcpy(&wv.vector, w->wts_robust);
      gsl_vector_mul(&wv.vector, w->wts_spatial);
    }
  else
    {
      gsl_vector_view wv = gsl_vector_subvector(w->wts_final, 0, w->nres);
      gsl_vector_memcpy(&wv.vector, w->wts_spatial);
    }

  if (!params->use_weights || params->synth_data)
    {
      gsl_vector_set_all(w->wts_final, 1.0);
      gsl_vector_set_all(w->sqrt_wts_final, 1.0);
    }
  else
    {
      size_t i;

      for (i = 0; i < w->nres; ++i)
        {
          double wi = gsl_vector_get(w->wts_final, i);
          gsl_vector_set(w->sqrt_wts_final, i, sqrt(wi));
        }
    }

#if FDF_SOLVER == 0 /* GSL multifit */

  fprintf(stderr, "mfield_calc_nonlinear: initializing multifit...");
  gettimeofday(&tv0, NULL);
  gsl_multifit_nlinear_winit(c, w->wts_final, &fdf, w->multifit_nlinear_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* compute initial residual */
  f = gsl_multifit_nlinear_residual(w->multifit_nlinear_p);
  res0 = gsl_blas_dnrm2(f);

  fprintf(stderr, "mfield_calc_nonlinear: computing nonlinear least squares solution...");
  gettimeofday(&tv0, NULL);
  s = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                                  mfield_nonlinear_callback, (void *) w,
                                  &info, w->multifit_nlinear_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (s == GSL_SUCCESS)
    {
      fprintf(stderr, "mfield_calc_nonlinear: NITER = %zu\n",
              gsl_multifit_nlinear_niter(w->multifit_nlinear_p));
      fprintf(stderr, "mfield_calc_nonlinear: NFEV  = %zu\n", fdf.nevalf);
      fprintf(stderr, "mfield_calc_nonlinear: NJEV  = %zu\n", fdf.nevaldf);
      fprintf(stderr, "mfield_calc_nonlinear: NAEV  = %zu\n", fdf.nevalfvv);
      fprintf(stderr, "mfield_calc_nonlinear: reason for stopping: %d\n", info);
      fprintf(stderr, "mfield_calc_nonlinear: initial |f(x)|: %.12e\n", res0);
      fprintf(stderr, "mfield_calc_nonlinear: final   |f(x)|: %.12e\n",
              gsl_blas_dnrm2(f));
    }
  else
    {
      fprintf(stderr, "mfield_calc_nonlinear: multifit failed: %s\n",
              gsl_strerror(s));
    }

  /* store final coefficients in physical units */
  {
    gsl_vector *x_final = gsl_multifit_nlinear_position(w->multifit_nlinear_p);
    gsl_vector_memcpy(w->c, x_final);
  }

#elif FDF_SOLVER == 1 /* GSL multilarge */

  s = mfield_calc_nonlinear_multilarge(c, w);

#elif FDF_SOLVER == 2 /* Gauss-Newton */

  s = mfield_calc_nonlinear_gn(c, w);

#endif

  gsl_vector_memcpy(c, w->c);

  printv_octave(c, "cfinal");

  w->niter++;

  return s;
}

int
mfield_calc_JTJ_multilarge(const gsl_vector *c, gsl_matrix * JTJ, mfield_workspace *w)
{
  int s;
  gsl_vector * u = gsl_vector_alloc(w->nres_tot);
  gsl_vector * v = gsl_vector_alloc(w->p);

  s = mfield_calc_df3(CblasTrans, c, u, w, v, JTJ);

  gsl_vector_free(u);
  gsl_vector_free(v);

  return s;
}

gsl_vector *
mfield_residual(const gsl_vector *c, mfield_workspace *w)
{
#if FDF_SOLVER == 0
  gsl_vector *f = gsl_multifit_nlinear_residual(w->multifit_nlinear_p);
#elif FDF_SOLVER == 1
  gsl_vector *f = gsl_multilarge_nlinear_residual(w->nlinear_workspace_p);
#elif FDF_SOLVER == 2
  gsl_vector *f = gsl_multilarge_nlinear_residual(w->nlinear_workspace_p);
#endif

  mfield_calc_f(c, w, f);

  return f;
}

/*
mfield_init_nonlinear()
  This function is called from mfield_init() to count
total number of residuals and allocate nonlinear least
squares workspaces

Notes:
1) weight_calc() must be called prior to this function to
compute spatial weights
*/

static int
mfield_init_nonlinear(mfield_workspace *w)
{
  int s = 0;
  const mfield_parameters *params = &(w->params);
  const size_t p = w->p;                 /* number of coefficients */
  size_t ndata = 0;                      /* number of distinct data points */
  size_t nres = 0;                       /* total number of residuals (data only) */
  size_t nres_B[4] = { 0, 0, 0, 0 };     /* number of (X,Y,Z,F) residuals */
  size_t nres_dB_ns[4] = { 0, 0, 0, 0 }; /* number of north-south gradient (X,Y,Z,F) residuals */
  size_t nres_dB_ew[4] = { 0, 0, 0, 0 }; /* number of east-west gradient (X,Y,Z,F) residuals */
  size_t nres_dBdt[4] = { 0, 0, 0, 0 };  /* number of secular variation (X,Y,Z,F) residuals */
  size_t i, j, k;

  /* count total number of residuals */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          /* check if data point is discarded due to time interval */
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* store starting residual index for this data point */
          mptr->index[j] = 0;
          for (k = 0; k < 4; ++k)
            {
              mptr->index[j] += nres_B[k] + nres_dB_ns[k] + nres_dB_ew[k] + nres_dBdt[k];
            }

          if (MAGDATA_ExistX(mptr->flags[j]))
            ++nres_B[0];
          if (MAGDATA_ExistY(mptr->flags[j]))
            ++nres_B[1];
          if (MAGDATA_ExistZ(mptr->flags[j]))
            ++nres_B[2];

          /* don't increase nres_F if only fitting Euler angles */
          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            ++nres_B[3];

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              if (MAGDATA_ExistDXDT(mptr->flags[j]))
                ++nres_dBdt[0];
              if (MAGDATA_ExistDYDT(mptr->flags[j]))
                ++nres_dBdt[1];
              if (MAGDATA_ExistDZDT(mptr->flags[j]))
                ++nres_dBdt[2];
            }

          if (MAGDATA_ExistDX_NS(mptr->flags[j]))
            ++nres_dB_ns[0];

          if (MAGDATA_ExistDY_NS(mptr->flags[j]))
            ++nres_dB_ns[1];

          if (MAGDATA_ExistDZ_NS(mptr->flags[j]))
            ++nres_dB_ns[2];

          if (MAGDATA_ExistDF_NS(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            ++nres_dB_ns[3];

          if (MAGDATA_ExistDX_EW(mptr->flags[j]))
            ++nres_dB_ew[0];

          if (MAGDATA_ExistDY_EW(mptr->flags[j]))
            ++nres_dB_ew[1];

          if (MAGDATA_ExistDZ_EW(mptr->flags[j]))
            ++nres_dB_ew[2];

          if (MAGDATA_ExistDF_EW(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            ++nres_dB_ew[3];

          ++ndata;
        }
    }

  if (ndata == 0)
    {
      GSL_ERROR("no data to fit", GSL_FAILURE);
    }

  for (k = 0; k < 4; ++k)
    {
      nres += nres_B[k] + nres_dB_ns[k] + nres_dB_ew[k] + nres_dBdt[k];
    }

  w->nres_vec = nres_B[0] + nres_B[1] + nres_B[2];
  w->nres_vec_SV = nres_dBdt[0] + nres_dBdt[1] + nres_dBdt[2];
  w->nres_vec_grad = nres_dB_ns[0] + nres_dB_ns[1] + nres_dB_ns[2] +
                     nres_dB_ew[0] + nres_dB_ew[1] + nres_dB_ew[2];
  w->nres = nres;
  w->ndata = ndata;

  /* check if we can use a linear least squares approach */
  if ((nres_B[3] + nres_dB_ns[3] + nres_dB_ew[3] + nres_dBdt[3]) == 0 &&
      params->fit_euler == 0)
    {
      w->lls_solution = 1;
    }
  else
    {
      w->lls_solution = 0;
    }

  /* compute total number of residuals, including regularization terms */
  w->nres_tot = nres;
  if (params->regularize && !params->synth_data)
    w->nres_tot += w->p;

  fprintf(stderr, "mfield_init_nonlinear: %zu distinct data points\n", ndata);
  fprintf(stderr, "mfield_init_nonlinear: %zu scalar residuals\n", nres_B[3]);
  fprintf(stderr, "mfield_init_nonlinear: %zu X residuals\n", nres_B[0]);
  fprintf(stderr, "mfield_init_nonlinear: %zu Y residuals\n", nres_B[1]);
  fprintf(stderr, "mfield_init_nonlinear: %zu Z residuals\n", nres_B[2]);
  fprintf(stderr, "mfield_init_nonlinear: %zu d/dt X residuals\n", nres_dBdt[0]);
  fprintf(stderr, "mfield_init_nonlinear: %zu d/dt Y residuals\n", nres_dBdt[1]);
  fprintf(stderr, "mfield_init_nonlinear: %zu d/dt Z residuals\n", nres_dBdt[2]);
  fprintf(stderr, "mfield_init_nonlinear: %zu dX_ns residuals\n", nres_dB_ns[0]);
  fprintf(stderr, "mfield_init_nonlinear: %zu dY_ns residuals\n", nres_dB_ns[1]);
  fprintf(stderr, "mfield_init_nonlinear: %zu dZ_ns residuals\n", nres_dB_ns[2]);
  fprintf(stderr, "mfield_init_nonlinear: %zu dF_ns residuals\n", nres_dB_ns[3]);
  fprintf(stderr, "mfield_init_nonlinear: %zu dX_ew residuals\n", nres_dB_ew[0]);
  fprintf(stderr, "mfield_init_nonlinear: %zu dY_ew residuals\n", nres_dB_ew[1]);
  fprintf(stderr, "mfield_init_nonlinear: %zu dZ_ew residuals\n", nres_dB_ew[2]);
  fprintf(stderr, "mfield_init_nonlinear: %zu dF_ew residuals\n", nres_dB_ew[3]);
  fprintf(stderr, "mfield_init_nonlinear: %zu data residuals\n", w->nres);
  fprintf(stderr, "mfield_init_nonlinear: %zu total residuals (including regularization terms)\n", w->nres_tot);
  fprintf(stderr, "mfield_init_nonlinear: %zu total parameters\n", p);

#if FDF_SOLVER == 0 /* GSL multifit */

  w->old_fdf = 1;

  /* allocate fit workspace */
  {
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_parameters fdf_params =
      gsl_multifit_nlinear_default_parameters();

    fdf_params.solver = gsl_multifit_nlinear_solver_cholesky;
    fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
    fdf_params.trs = gsl_multifit_nlinear_trs_lm;
    fdf_params.fdtype = GSL_MULTIFIT_NLINEAR_CTRDIFF;
    w->multifit_nlinear_p = gsl_multifit_nlinear_alloc(T, &fdf_params, w->nres_tot, p);
  }

#elif FDF_SOLVER == 1 /* GSL multilarge */

  /* allocate fit workspace - start with Levenberg-Marquardt solver */
  mfield_nonlinear_alloc_multilarge(gsl_multilarge_nlinear_trs_lm, w);

#elif FDF_SOLVER == 2 /* Gauss-Newton */

  mfield_nonlinear_alloc_gn(w);

#endif

  /* allocate sparse Jacobian matrix for Euler angles, fluxgate calibration and external field */
  w->J2 = gsl_spmatrix_alloc(nres, GSL_MAX(w->p_sparse, 1));
  
  w->wts_spatial = gsl_vector_alloc(nres);
  w->wts_robust = gsl_vector_alloc(nres);
  w->wts_final = gsl_vector_alloc(w->nres_tot);
  w->sqrt_wts_final = gsl_vector_alloc(w->nres_tot);

  gsl_vector_set_all(w->wts_robust, 1.0);
  gsl_vector_set_all(w->wts_final, 1.0);
  gsl_vector_set_all(w->sqrt_wts_final, 1.0);

  /* to save memory, allocate p-by-p workspace since a full n-by-p isn't needed */
  w->robust_workspace_p = gsl_multifit_robust_alloc(gsl_multifit_robust_bisquare, p, p);

#if 0
  /* suggested by Nils to use 1.5 tuning constant */
  gsl_multifit_robust_tune(1.5, w->robust_workspace_p);
#endif

  w->fvec = gsl_vector_alloc(w->nres_tot);
  w->wfvec = gsl_vector_alloc(w->nres_tot);

  /* calculate spatial weights */
  {
    size_t idx = 0;
    size_t j;

    fprintf(stderr, "mfield_init_nonlinear: calculating spatial weights...");

    for (i = 0; i < w->nsat; ++i)
      {
        magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
        double global_weight = 1.0; /* global weight for this data source */

        /* down-weight DMSP data relative to CHAMP/Swarm; this should be a cfg file parameter */
        if (mptr->global_flags & MAGDATA_GLOBFLG_DMSP)
          global_weight = 0.01;
        else if (mptr->global_flags & (MAGDATA_GLOBFLG_OBSERVATORY | MAGDATA_GLOBFLG_OBSERVATORY_SV))
          global_weight = 2.0;

        for (j = 0; j < mptr->n; ++j)
          {
            double wt; /* spatial weight */

            if (MAGDATA_Discarded(mptr->flags[j]))
              continue;

            if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
              spatwt_get(mptr->theta[j], mptr->phi[j], &wt, w->spatwtMF_workspace_p);
            else if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY_SV)
              spatwt_get(mptr->theta[j], mptr->phi[j], &wt, w->spatwtSV_workspace_p);
            else if (mptr->global_flags & MAGDATA_GLOBFLG_EEJ_MAGEQ)
              wt = 1.0;
            else /* satellite data */
              track_weight_get(mptr->phi[j], mptr->theta[j], &wt, w->weight_workspace_p);

            /* include global weight factor for this satellite / observatory */
            wt *= global_weight;

            if (MAGDATA_ExistX(mptr->flags[j]))
              gsl_vector_set(w->wts_spatial, idx++, params->weight_X * wt);

            if (MAGDATA_ExistY(mptr->flags[j]))
              gsl_vector_set(w->wts_spatial, idx++, params->weight_Y * wt);

            if (MAGDATA_ExistZ(mptr->flags[j]))
              gsl_vector_set(w->wts_spatial, idx++, params->weight_Z * wt);

            if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
              gsl_vector_set(w->wts_spatial, idx++, params->weight_F * wt);

            if (MAGDATA_FitMF(mptr->flags[j]))
              {
                if (MAGDATA_ExistDXDT(mptr->flags[j]))
                  gsl_vector_set(w->wts_spatial, idx++, params->weight_DXDT * wt);

                if (MAGDATA_ExistDYDT(mptr->flags[j]))
                  gsl_vector_set(w->wts_spatial, idx++, params->weight_DYDT * wt);

                if (MAGDATA_ExistDZDT(mptr->flags[j]))
                  gsl_vector_set(w->wts_spatial, idx++, params->weight_DZDT * wt);
              }

            if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
              gsl_vector_set(w->wts_spatial, idx++, params->weight_DX * wt);

            if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
              gsl_vector_set(w->wts_spatial, idx++, params->weight_DY * wt);

            if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
              gsl_vector_set(w->wts_spatial, idx++, params->weight_DZ * wt);

            if (MAGDATA_ExistDF_NS(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
              gsl_vector_set(w->wts_spatial, idx++, params->weight_F * wt);

            if (MAGDATA_ExistDF_EW(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
              gsl_vector_set(w->wts_spatial, idx++, params->weight_F * wt);
          }
      }

    fprintf(stderr, "done\n");

    assert(idx == w->nres);
  }

  /* precompute regularization matrix */
  if (params->regularize && !params->synth_data)
    {
      fprintf(stderr, "mfield_init_nonlinear: calculating regularization matrix...");
      mfield_nonlinear_regularize_init(w);
      fprintf(stderr, "done\n");
    }
  else
    {
      gsl_spmatrix_set_zero(w->L);
      gsl_spmatrix_set_zero(w->Lambda);
    }

  return s;
} /* mfield_init_nonlinear() */

#if 0
/*
mfield_calc_df()
  Compute J^T J matrix and J^T u or J u vector using OpenMP
for speed improvement
*/

static int
mfield_calc_df2(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                void *params, gsl_vector * v, gsl_matrix * JTJ)
{
  mfield_workspace *w = (mfield_workspace *) params;
  size_t i, j;
  gsl_matrix_view JTJ_int; /* internal field portion of J^T J */
  struct timeval tv0, tv1;

  gettimeofday(&tv0, NULL);
  mfield_debug("mfield_calc_df2: entering function...\n");
  mfield_debug("mfield_calc_df2: TransJ = %s...\n",
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
      gsl_matrix_tricpy(CblasLower, CblasNonUnit, &JTJ_int.matrix, w->JTJ_vec);
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
      int fit_euler = w->params.fit_euler && (mptr->global_flags & MAGDATA_GLOBFLG_EULER);

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
          double B_int[3];              /* internal field model */
          double B_model[3];            /* a priori model (crustal/external) */
          double B_total[4];            /* internal + external */
          double dB_ext[3] = { 0.0, 0.0, 0.0 };
          double B_extcorr[3] = { 0.0, 0.0, 0.0 }; /* external field correction model */
          double B_nec_alpha[3], B_nec_beta[3], B_nec_gamma[3];
          size_t extidx = 0;
          size_t euler_idx = 0;

          gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
          gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
          gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);

          gsl_vector_view vx_grad = gsl_matrix_row(w->omp_dX_grad, thread_id);
          gsl_vector_view vy_grad = gsl_matrix_row(w->omp_dY_grad, thread_id);
          gsl_vector_view vz_grad = gsl_matrix_row(w->omp_dZ_grad, thread_id);

#if MFIELD_FIT_EXTFIELD
          double extcoeff = 0.0;
#endif
          double B_vfm[3];        /* observation vector VFM frame */

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

          /* compute internal field model */
          B_int[0] = mfield_nonlinear_model_int(t, &vx.vector, x, w);
          B_int[1] = mfield_nonlinear_model_int(t, &vy.vector, x, w);
          B_int[2] = mfield_nonlinear_model_int(t, &vz.vector, x, w);

          /* load apriori model of external (and possibly crustal) field */
          B_model[0] = mptr->Bx_model[j];
          B_model[1] = mptr->By_model[j];
          B_model[2] = mptr->Bz_model[j];

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

              /* add correction to external field model */
              for (k = 0; k < 3; ++k)
                B_extcorr[k] = extcoeff * dB_ext[k];
            }
#endif

          /* compute total modeled field (internal + external) */
          for (k = 0; k < 3; ++k)
            B_total[k] = B_int[k] + B_model[k] + B_extcorr[k];

          /* compute Euler angle derivatives of B vector */
          if (fit_euler)
            {
              const double *q = &(mptr->q[4*j]);
              double alpha, beta, gamma;

              euler_idx = mfield_euler_idx(i, mptr->t[j], w);
              alpha = gsl_vector_get(x, euler_idx);
              beta = gsl_vector_get(x, euler_idx + 1);
              gamma = gsl_vector_get(x, euler_idx + 2);

              /* get vector in VFM frame */
              B_vfm[0] = mptr->Bx_vfm[j];
              B_vfm[1] = mptr->By_vfm[j];
              B_vfm[2] = mptr->Bz_vfm[j];

              /* compute alpha derivative of: R_q R_3 B_vfm */
              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_ALPHA, alpha, beta, gamma, q, B_vfm, B_nec_alpha);

              /* compute beta derivative of: R_q R_3 B_vfm */
              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_BETA, alpha, beta, gamma, q, B_vfm, B_nec_beta);

              /* compute gamma derivative of: R_q R_3 B_vfm */
              euler_vfm2nec(mptr->euler_flags | EULER_FLG_DERIV_GAMMA, alpha, beta, gamma, q, B_vfm, B_nec_gamma);
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double wj = gsl_vector_get(w->wts_final, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_JTu(t, mptr->flags[j], wj, u, ridx, &vx.vector,
                                        extidx, dB_ext[0], euler_idx, B_nec_alpha[0],
                                        B_nec_beta[0], B_nec_gamma[0], v, w);
                  }
                }
              else
                {
#pragma omp critical
                  {
                    mfield_jacobian_Ju(t, mptr->flags[j], wj, u, ridx, &vx.vector,
                                       extidx, dB_ext[0], euler_idx, B_nec_alpha[0],
                                       B_nec_beta[0], B_nec_gamma[0], v, w);
                  }
                }

#pragma omp critical
              {
                mfield_jacobian_JTJ(t, mptr->flags[j], wj, &vx.vector,
                                    extidx, dB_ext[0], euler_idx, B_nec_alpha[0],
                                    B_nec_beta[0], B_nec_gamma[0], JTJ, w);
              }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double wj = gsl_vector_get(w->wts_final, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_JTu(t, mptr->flags[j], wj, u, ridx, &vy.vector,
                                        extidx, dB_ext[1], euler_idx, B_nec_alpha[1],
                                        B_nec_beta[1], B_nec_gamma[1], v, w);
                  }
                }
              else
                {
#pragma omp critical
                  {
                    mfield_jacobian_Ju(t, mptr->flags[j], wj, u, ridx, &vy.vector,
                                       extidx, dB_ext[1], euler_idx, B_nec_alpha[1],
                                       B_nec_beta[1], B_nec_gamma[1], v, w);
                  }
                }

#pragma omp critical
              {
                mfield_jacobian_JTJ(t, mptr->flags[j], wj, &vy.vector,
                                    extidx, dB_ext[1], euler_idx, B_nec_alpha[1],
                                    B_nec_beta[1], B_nec_gamma[1], JTJ, w);
              }

              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double wj = gsl_vector_get(w->wts_final, ridx);

              if (TransJ == CblasTrans)
                {
#pragma omp critical
                  {
                    mfield_jacobian_JTu(t, mptr->flags[j], wj, u, ridx, &vz.vector,
                                        extidx, dB_ext[2], euler_idx, B_nec_alpha[2],
                                        B_nec_beta[2], B_nec_gamma[2], v, w);
                  }
                }
              else
                {
#pragma omp critical
                  {
                    mfield_jacobian_Ju(t, mptr->flags[j], wj, u, ridx, &vz.vector,
                                       extidx, dB_ext[2], euler_idx, B_nec_alpha[2],
                                       B_nec_beta[2], B_nec_gamma[2], v, w);
                  }
                }

#pragma omp critical
              {
                mfield_jacobian_JTJ(t, mptr->flags[j], wj, &vz.vector,
                                    extidx, dB_ext[2], euler_idx, B_nec_alpha[2],
                                    B_nec_beta[2], B_nec_gamma[2], JTJ, w);
              }

              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double wj = gsl_vector_get(w->wts_final, ridx);
              gsl_vector_view Jv = gsl_matrix_subrow(w->omp_J[thread_id], w->omp_rowidx[thread_id]++, 0, w->p_int);

              B_total[3] = gsl_hypot3(B_total[0], B_total[1], B_total[2]);

#pragma omp critical
              {
                mfield_jacobian_row_F(TransJ, t, wj, u, ridx, &vx.vector, &vy.vector, &vz.vector,
                                      B_total, extidx, dB_ext, &Jv.vector, JTJ, v, w);
              }

              ++ridx;
            }

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
  mfield_debug("mfield_calc_df2: leaving function... (%g seconds)\n", time_diff(tv0, tv1));

  return GSL_SUCCESS;
}

/*
mfield_jacobian_JTu()
  Update the J^T u vector with a new row of the Jacobian matrix,
corresponding to a vector residual.

Inputs: t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
        u           - input vector u, size n
        ridx        - residual index of this row in [0,nres-1]
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_max-by-1
        extidx      - index of external field coefficient in [0,p_ext-1]
        dB_ext      - external field Green's function corresponding
                      to desired vector component
        euler_idx   - index of Euler angles
        B_nec_alpha - Green's function for alpha Euler angle
        B_nec_beta  - Green's function for beta Euler angle
        B_nec_gamma - Green's function for gamma Euler angle
        JTu         - (output) J^T y vector
        w           - workspace
*/

static inline int
mfield_jacobian_JTu(const double t, const size_t flags, const double weight,
                    const gsl_vector * u, const size_t ridx, const gsl_vector * dB_int, const size_t extidx,
                    const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                    const double B_nec_beta, const double B_nec_gamma,
                    gsl_vector *JTu, const mfield_workspace *w)
{
  const double y = gsl_vector_get(u, ridx);
  const double sWy = sqrt(weight) * y;

  (void) extidx;
  (void) dB_ext;
  (void) euler_idx;
  (void) B_nec_alpha;
  (void) B_nec_beta;
  (void) B_nec_gamma;

  /* check if fitting MF to this data point */
  if (MAGDATA_FitMF(flags))
    {
      gsl_vector_view v;

      /* update J^T y */

      if (w->nnm_mf > 0)
        {
          gsl_vector_const_view g_mf = gsl_vector_const_subvector(dB_int, 0, w->nnm_mf);

          v = gsl_vector_subvector(JTu, 0, w->nnm_mf);
          gsl_blas_daxpy(sWy, &g_mf.vector, &v.vector);
        }

      if (w->nnm_sv > 0)
        {
          gsl_vector_const_view g_sv = gsl_vector_const_subvector(dB_int, 0, w->nnm_sv);

          v = gsl_vector_subvector(JTu, w->sv_offset, w->nnm_sv);
          gsl_blas_daxpy(t * sWy, &g_sv.vector, &v.vector);
        }

      if (w->nnm_sa > 0)
        {
          gsl_vector_const_view g_sa = gsl_vector_const_subvector(dB_int, 0, w->nnm_sa);

          v = gsl_vector_subvector(JTu, w->sa_offset, w->nnm_sa);
          gsl_blas_daxpy(0.5 * t * t * sWy, &g_sa.vector, &v.vector);
        }

#if MFIELD_FIT_EXTFIELD
      {
        double *ptr = gsl_vector_ptr(JTu, extidx);

        /* update J^T y */
        *ptr += dB_ext * sWy;
      }
#endif /* MFIELD_FIT_EXTFIELD */
    }

  /* check if fitting Euler angles to this data point */
  if (w->params.fit_euler && MAGDATA_FitEuler(flags))
    {
      double x_data[3];
      gsl_vector_view vJTu = gsl_vector_subvector(JTu, euler_idx, 3);
      gsl_vector_view v = gsl_vector_view_array(x_data, 3);

      x_data[0] = -B_nec_alpha;
      x_data[1] = -B_nec_beta;
      x_data[2] = -B_nec_gamma;

      /* update J^T y */
      gsl_blas_daxpy(sWy, &v.vector, &vJTu.vector);
    }

  return GSL_SUCCESS;
}

/*
mfield_jacobian_SV_JTu()
  Update the J^T u vector with a new row of the Jacobian matrix,
corresponding to a secular variation vector residual.

Inputs: t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
        u           - input vector u, size n
        ridx        - residual index of this row in [0,nres-1]
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_max-by-1
        JTu         - (output) J^T y vector
        w           - workspace
*/

static inline int
mfield_jacobian_SV_JTu(const double t, const size_t flags, const double weight,
                       const gsl_vector * u, const size_t ridx, const gsl_vector * dB_int,
                       gsl_vector *JTu, const mfield_workspace *w)
{
  const double y = gsl_vector_get(u, ridx);
  const double sWy = sqrt(weight) * y;

  /* check if fitting MF to this data point */
  if (MAGDATA_FitMF(flags))
    {
      gsl_vector_view v;

      /* update J^T y */

      if (w->nnm_sv > 0)
        {
          gsl_vector_const_view g_sv = gsl_vector_const_subvector(dB_int, 0, w->nnm_sv);

          v = gsl_vector_subvector(JTu, w->sv_offset, w->nnm_sv);
          gsl_blas_daxpy(sWy, &g_sv.vector, &v.vector);
        }

      if (w->nnm_sa > 0)
        {
          gsl_vector_const_view g_sa = gsl_vector_const_subvector(dB_int, 0, w->nnm_sa);

          v = gsl_vector_subvector(JTu, w->sa_offset, w->nnm_sa);
          gsl_blas_daxpy(t * sWy, &g_sa.vector, &v.vector);
        }
    }

  return GSL_SUCCESS;
}

/*
mfield_jacobian_grad_JTu()
  Update the J^T u vector with a new row of the Jacobian matrix,
corresponding to a vector gradient residual.

Inputs: t           - scaled timestamp
        t_grad      - scaled timestamp of gradient point (N/S or E/W)
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
        u           - input vector u, size n
        ridx        - residual index of this row in [0,nres-1]
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_max-by-1
        dB_int_grad - Green's functions for desired vector gradient component of
                      internal SH expansion, nnm_max-by-1
        JTu         - (output) J^T y vector
        w           - workspace
*/

static inline int
mfield_jacobian_grad_JTu(const double t, const double t_grad, const size_t flags, const double weight,
                         const gsl_vector * u, const size_t ridx, const gsl_vector * dB_int, const gsl_vector * dB_int_grad,
                         gsl_vector *JTu, const mfield_workspace *w)
{
  /* check if fitting MF to this data point */
  if (MAGDATA_FitMF(flags))
    {
      const double y = gsl_vector_get(u, ridx);
      const double sWy = sqrt(weight) * y;
      gsl_vector_view v;

      /* update J^T y */

      if (w->nnm_mf > 0)
        {
          gsl_vector_const_view g_mf = gsl_vector_const_subvector(dB_int, 0, w->nnm_mf);
          gsl_vector_const_view dg_mf = gsl_vector_const_subvector(dB_int_grad, 0, w->nnm_mf);

          v = gsl_vector_subvector(JTu, 0, w->nnm_mf);
          gsl_blas_daxpy(sWy, &dg_mf.vector, &v.vector);
          gsl_blas_daxpy(-sWy, &g_mf.vector, &v.vector);
        }

      if (w->nnm_sv > 0)
        {
          gsl_vector_const_view g_sv = gsl_vector_const_subvector(dB_int, 0, w->nnm_sv);
          gsl_vector_const_view dg_sv = gsl_vector_const_subvector(dB_int_grad, 0, w->nnm_sv);

          v = gsl_vector_subvector(JTu, w->sv_offset, w->nnm_sv);
          gsl_blas_daxpy(t_grad * sWy, &dg_sv.vector, &v.vector);
          gsl_blas_daxpy(-t * sWy, &g_sv.vector, &v.vector);
        }

      if (w->nnm_sa > 0)
        {
          gsl_vector_const_view g_sa = gsl_vector_const_subvector(dB_int, 0, w->nnm_sa);
          gsl_vector_const_view dg_sa = gsl_vector_const_subvector(dB_int_grad, 0, w->nnm_sa);

          v = gsl_vector_subvector(JTu, w->sa_offset, w->nnm_sa);
          gsl_blas_daxpy(0.5 * t_grad * t_grad * sWy, &dg_sa.vector, &v.vector);
          gsl_blas_daxpy(-0.5 * t * t * sWy, &g_sa.vector, &v.vector);
        }
    }

  return GSL_SUCCESS;
}

/*
mfield_jacobian_Ju()
  Update the J u vector with a new row of the Jacobian matrix,
corresponding to a vector residual.

For this row of J u (specified by ridx), the value is given
by:

(J u)_i = sqrt(w_i) J_i^T u

where J_i^T is row i of the Jacobian (1-by-p)

Inputs: t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
        u           - input vector u, size p
        ridx        - residual index of this row in [0,nres-1]
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_max-by-1
        extidx      - index of external field coefficient in [0,p_ext-1]
        dB_ext      - external field Green's function corresponding
                      to desired vector component
        euler_idx   - index of Euler angles
        B_nec_alpha - Green's function for alpha Euler angle
        B_nec_beta  - Green's function for beta Euler angle
        B_nec_gamma - Green's function for gamma Euler angle
        Ju          - (output) J u vector
        w           - workspace
*/

static inline int
mfield_jacobian_Ju(const double t, const size_t flags, const double weight,
                   const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, const size_t extidx,
                   const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                   const double B_nec_beta, const double B_nec_gamma,
                   gsl_vector *Ju, const mfield_workspace *w)
{
  double *Ju_ptr = gsl_vector_ptr(Ju, ridx);

  (void) extidx;
  (void) dB_ext;
  (void) euler_idx;
  (void) B_nec_alpha;
  (void) B_nec_beta;
  (void) B_nec_gamma;

  /* check if fitting MF to this data point */
  if (MAGDATA_FitMF(flags))
    {
      double tmp;

      /* update J u */

      if (w->nnm_mf > 0)
        {
          gsl_vector_view g_mf = gsl_vector_subvector(dB_int, 0, w->nnm_mf);
          gsl_vector_const_view u_mf = gsl_vector_const_subvector(u, 0, w->nnm_mf);

          gsl_blas_ddot(&g_mf.vector, &u_mf.vector, &tmp);
          *Ju_ptr = tmp;
        }

      if (w->nnm_sv > 0)
        {
          gsl_vector_view g_sv = gsl_vector_subvector(dB_int, 0, w->nnm_sv);
          gsl_vector_const_view u_sv = gsl_vector_const_subvector(u, w->sv_offset, w->nnm_sv);

          gsl_blas_ddot(&g_sv.vector, &u_sv.vector, &tmp);
          *Ju_ptr += t * tmp;
        }

      if (w->nnm_sa > 0)
        {
          gsl_vector_view g_sa = gsl_vector_subvector(dB_int, 0, w->nnm_sa);
          gsl_vector_const_view u_sa = gsl_vector_const_subvector(u, w->sa_offset, w->nnm_sa);

          gsl_blas_ddot(&g_sa.vector, &u_sa.vector, &tmp);
          *Ju_ptr += 0.5 * t * t * tmp;
        }

#if MFIELD_FIT_EXTFIELD

      /* update J u */
      *Ju_ptr += dB_ext * gsl_vector_get(u, extidx);

#endif /* MFIELD_FIT_EXTFIELD */
    }

  /* check if fitting Euler angles to this data point */
  if (w->params.fit_euler && MAGDATA_FitEuler(flags))
    {
      double x_data[3];
      gsl_vector_const_view vu = gsl_vector_const_subvector(u, euler_idx, 3);
      gsl_vector_view v = gsl_vector_view_array(x_data, 3);
      double tmp;

      x_data[0] = -B_nec_alpha;
      x_data[1] = -B_nec_beta;
      x_data[2] = -B_nec_gamma;

      gsl_blas_ddot(&vu.vector, &v.vector, &tmp);
      *Ju_ptr += tmp;
    }

  *Ju_ptr *= sqrt(weight);

  return GSL_SUCCESS;
}

/*
mfield_jacobian_JTJ()
  Update the J^T J matrix with a new row of the Jacobian matrix,
corresponding to a vector residual. The internal field portion of J^T J
does not need to be computed, since it is independent of the model parameters
and is pre-computed. Only the Euler and external field
portion of J^T J must be updated.

Inputs: t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
        u           - input vector u
        ridx        - residual index of this row in [0,nres-1]
        dB_int      - Green's functions for desired vector component of
                      internal SH expansion, nnm_max-by-1
        extidx      - index of external field coefficient in [0,p_ext-1]
        dB_ext      - external field Green's function corresponding
                      to desired vector component
        euler_idx   - index of Euler angles
        B_nec_alpha - Green's function for alpha Euler angle
        B_nec_beta  - Green's function for beta Euler angle
        B_nec_gamma - Green's function for gamma Euler angle
        JTJ         - (output) J^T J matrix, possibly NULL
        w           - workspace
*/

static inline int
mfield_jacobian_JTJ(const double t, const size_t flags, const double weight,
                    gsl_vector * dB_int, const size_t extidx,
                    const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                    const double B_nec_beta, const double B_nec_gamma,
                    gsl_matrix *JTJ, const mfield_workspace *w)
{
  gsl_vector_view g_mf, g_sv, g_sa;

  (void) extidx;
  (void) dB_ext;
  (void) euler_idx;
  (void) B_nec_alpha;
  (void) B_nec_beta;
  (void) B_nec_gamma;

  /* check for quick return */
  if (JTJ == NULL)
    return GSL_SUCCESS;

  if (w->nnm_mf > 0)
    g_mf = gsl_vector_subvector(dB_int, 0, w->nnm_mf);

  if (w->nnm_sv > 0)
    g_sv = gsl_vector_subvector(dB_int, 0, w->nnm_sv);

  if (w->nnm_sa > 0)
    g_sa = gsl_vector_subvector(dB_int, 0, w->nnm_sa);

  /* check if fitting MF to this data point */
  if (MAGDATA_FitMF(flags))
    {
#if MFIELD_FIT_EXTFIELD

      /* update J^T J */
      gsl_vector_view v;
      double *ptr33 = gsl_matrix_ptr(JTJ, extidx, extidx);

      /* update (J^T J)_33 */
      *ptr33 += dB_ext * dB_ext * weight;

      /* update (J^T J)_31 = J_ext^T J_int */

      if (w->nnm_mf > 0)
        {
          v = gsl_matrix_subrow(JTJ, extidx, 0, w->nnm_mf);
          gsl_blas_daxpy(dB_ext * weight, &g_mf.vector, &v.vector);
        }

      if (w->nnm_sv > 0)
        {
          v = gsl_matrix_subrow(JTJ, extidx, w->sv_offset, w->nnm_sv);
          gsl_blas_daxpy(t * dB_ext * weight, &g_sv.vector, &v.vector);
        }

      if (w->nnm_sa > 0)
        {
          v = gsl_matrix_subrow(JTJ, extidx, w->sa_offset, w->nnm_sa);
          gsl_blas_daxpy(0.5 * t * t * dB_ext * weight, &g_sa.vector, &v.vector);
        }

#endif /* MFIELD_FIT_EXTFIELD */
    }

  /* check if fitting Euler angles to this data point */
  if (w->params.fit_euler && MAGDATA_FitEuler(flags))
    {
      double x_data[3];
      gsl_vector_view v = gsl_vector_view_array(x_data, 3);
      gsl_matrix_view m;

      x_data[0] = -B_nec_alpha;
      x_data[1] = -B_nec_beta;
      x_data[2] = -B_nec_gamma;

      /* update (J^T J)_22 */
      m = gsl_matrix_submatrix(JTJ, euler_idx, euler_idx, 3, 3);
      gsl_blas_dsyr(CblasLower, weight, &v.vector, &m.matrix);

      if (MAGDATA_FitMF(flags))
        {
          /* update (J^T J)_21 */

          if (w->nnm_mf > 0)
            {
              m = gsl_matrix_submatrix(JTJ, euler_idx, 0, 3, w->nnm_mf);
              gsl_blas_dger(weight, &v.vector, &g_mf.vector, &m.matrix);
            }

          if (w->nnm_sv > 0)
            {
              m = gsl_matrix_submatrix(JTJ, euler_idx, w->sv_offset, 3, w->nnm_sv);
              gsl_blas_dger(t * weight, &v.vector, &g_sv.vector, &m.matrix);
            }

          if (w->nnm_sa > 0)
            {
              m = gsl_matrix_submatrix(JTJ, euler_idx, w->sa_offset, 3, w->nnm_sa);
              gsl_blas_dger(0.5 * t * t * weight, &v.vector, &g_sa.vector, &m.matrix);
            }

#if MFIELD_FIT_EXTFIELD
          /* update (J^T J)_32 */
          {
            gsl_vector_view v32 = gsl_matrix_subrow(JTJ, extidx, euler_idx, 3);
            gsl_blas_daxpy(dB_ext * weight, &v.vector, &v32.vector);
          }
#endif
        }
    }

  return GSL_SUCCESS;
}

/*
mfield_jacobian_row_F()
  Construct a row of the Jacobian matrix corresponding to
a scalar measurement and update J^T J matrix and op(J) u
vector

Inputs: TransJ      - op(J)
        t           - scaled timestamp
        flags       - MAGDATA_FLG_xxx flags for this data point
        weight      - weight for this data point
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
        extidx      - index of external field coefficient
        dB_ext      - external field vector Green's functions
        J_int       - (output) row of Jacobian (weighted) for internal
                               Green's functions, p_int-by-1
        JTJ         - (output) updated J^T W J matrix, possibly NULL
        v           - (output) updated op(J) u vector
        w           - workspace
*/

static inline int
mfield_jacobian_row_F(CBLAS_TRANSPOSE_t TransJ, const double t, const double weight,
                      const gsl_vector * u, const size_t ridx,
                      gsl_vector * dX, gsl_vector * dY, gsl_vector * dZ,
                      const double B_model[4], const size_t extidx, const double dB_ext[3],
                      gsl_vector *J_int, gsl_matrix *JTJ, gsl_vector *v,
                      const mfield_workspace *w)
{
  const double sqrt_weight = sqrt(weight);
  size_t k;
  double b[3], ui, *Ju_ptr;

  (void) extidx;
  (void) dB_ext;
  (void) JTJ;

  if (TransJ == CblasTrans)
    {
      ui = gsl_vector_get(u, ridx);
    }
  else
    {
      Ju_ptr = gsl_vector_ptr(v, ridx);
    }

  /* compute unit vector in model direction */
  for (k = 0; k < 3; ++k)
    b[k] = B_model[k] / B_model[3];

  /* compute (X dX + Y dY + Z dZ) */
  for (k = 0; k < w->nnm_max; ++k)
    {
      double dXk = gsl_vector_get(dX, k);
      double dYk = gsl_vector_get(dY, k);
      double dZk = gsl_vector_get(dZ, k);
      double val = sqrt_weight * (b[0] * dXk +
                                  b[1] * dYk +
                                  b[2] * dZk);

      mfield_set_mf(J_int, k, val, w);
      mfield_set_sv(J_int, k, t * val, w);
      mfield_set_sa(J_int, k, 0.5 * t * t * val, w);
    }

  if (TransJ == CblasTrans)
    {
      /* update J^T u */
      gsl_vector_view vJTu = gsl_vector_subvector(v, 0, w->p_int);
      gsl_blas_daxpy(ui, J_int, &vJTu.vector);
    }
  else
    {
      /* update (J u)_i = J_int . u(1:pint) */
      gsl_vector_const_view z = gsl_vector_const_subvector(u, 0, w->p_int);
      gsl_blas_ddot(J_int, &z.vector, Ju_ptr);
    }

#if MFIELD_FIT_EXTFIELD
  {
    double val = sqrt_weight * (b[0] * dB_ext[0] +
                                b[1] * dB_ext[1] +
                                b[2] * dB_ext[2]);

    if (TransJ == CblasTrans)
      {
        /* update J^T u */
        double *ptr = gsl_vector_ptr(v, extidx);
        *ptr += val * ui;
      }
    else
      {
        double tmp = gsl_vector_get(u, extidx);
        *Ju_ptr += val * tmp;
      }

    /* update J^T J */
    if (JTJ)
      {
        gsl_vector_view v31 = gsl_matrix_subrow(JTJ, extidx, 0, w->p_int);
        double *ptr33 = gsl_matrix_ptr(JTJ, extidx, extidx);

        /* update (J^T J)_33 */
        *ptr33 += val * val;

        /* update (J^T J)_31 */
        gsl_blas_daxpy(val, J_int, &v31.vector);
      }
  }
#endif

  return GSL_SUCCESS;
}

/*
mfield_vector_green_SV()
  Function to compute sqrt(w) [ 0 J_sv J_sa ] for a given set of
vector Green's functions for SV calculation

Inputs: t      - scaled timestamp
        weight - weight for this measurement
        g      - vector Green's functions J_mf, size nnm_max
        G      - (output) combined vector G = sqrt(w) [ 0 ; g ; t*g ],
                 size w->p_int
        w      - workspace
*/

static int
mfield_vector_green_SV(const double t, const double weight, const gsl_vector *g,
                    gsl_vector *G, mfield_workspace *w)
{
  const double sqrt_weight = sqrt(weight);
  size_t i;

  /* form G */
  for (i = 0; i < w->nnm_max; ++i)
    {
      double gi = sqrt_weight * gsl_vector_get(g, i);

      mfield_set_mf(G, i, 0.0, w);
      mfield_set_sv(G, i, gi, w);
      mfield_set_sa(G, i, t * gi, w);
    }

  return GSL_SUCCESS;
}

/*
mfield_vector_green_grad()
  Function to compute sqrt(w) [ dJ_mf ; dJ_sv ; dJ_sa ] for a given set of
vector Green's functions for a point and a gradient point

Inputs: t      - scaled timestamp
        t_grad - scaled timestamp of gradient point (N/S or E/W)
        weight - weight for this measurement
        g      - vector Green's functions J_mf, size nnm_max
        g_grad - vector Green's functions of gradient point, size nnm_max
        G      - (output) combined vector G = sqrt(w) [ gj - gi ; tj*gj - ti*gi ; 0.5*tj*tj*gj - 0.5*ti*ti*gi ],
                 size w->p_int
        w      - workspace
*/

static int
mfield_vector_green_grad(const double t, const double t_grad, const double weight, const gsl_vector *g,
                         const gsl_vector *g_grad, gsl_vector *G, mfield_workspace *w)
{
  const double sqrt_weight = sqrt(weight);
  size_t i;

  /* form G */
  for (i = 0; i < w->nnm_max; ++i)
    {
      double gi = sqrt_weight * gsl_vector_get(g, i);
      double gj = sqrt_weight * gsl_vector_get(g_grad, i);

      mfield_set_mf(G, i, gj - gi, w);
      mfield_set_sv(G, i, t_grad * gj - t * gi, w);
      mfield_set_sa(G, i, 0.5 * t_grad * t_grad * gj - 0.5 * t * t * gi, w);
    }

  return GSL_SUCCESS;
}

#endif /* 0 */

/*
mfield_nonlinear_model_ext()
  Compute external field model:

r_k * (0.7 * external_dipole + 0.3 * internal_dipole)

where the dipole coefficients are aligned with the
direction of the main field dipole. r_k is the strength
of the ring current for a given day (k = doy) but is not
incorporated into the output of this function. This function
outputs the coefficient of r_k, ie, the "green's function" part.

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        g     - model coefficients
        dB    - (output) external field model green's functions
                dB[0] = X component of external field
                dB[1] = Y component of external field
                dB[2] = Z component of external field
        w     - workspace
*/

static int
mfield_nonlinear_model_ext(const double r, const double theta,
                           const double phi, const gsl_vector *g,
                           double dB[3], const mfield_workspace *w)
{
#if MFIELD_FIT_EXTFIELD

  int s = 0;
  double g10 = gsl_vector_get(g, mfield_coeff_nmidx(1, 0));
  double g11 = gsl_vector_get(g, mfield_coeff_nmidx(1, 1));
  double h11 = gsl_vector_get(g, mfield_coeff_nmidx(1, -1));
  double g1 = gsl_hypot3(g10, g11, h11);
  double q[3];
  mfield_green_workspace *green_p = mfield_green_alloc(1, w->R);

  /* construct unit vector along internal dipole direction */
  q[mfield_coeff_nmidx(1, 0)] = g10 / g1;
  q[mfield_coeff_nmidx(1, 1)] = g11 / g1;
  q[mfield_coeff_nmidx(1, -1)] = h11 / g1;

  /* compute internal and external dipole field components */
  mfield_green_calc(r, theta, phi, green_p);
  mfield_green_ext(r, theta, phi, green_p);

  /* add external and induced sums */
  dB[0] = 0.7*vec_dot(q, green_p->dX_ext) + 0.3*vec_dot(q, green_p->dX);
  dB[1] = 0.7*vec_dot(q, green_p->dY_ext) + 0.3*vec_dot(q, green_p->dY);
  dB[2] = 0.7*vec_dot(q, green_p->dZ_ext) + 0.3*vec_dot(q, green_p->dZ);

  mfield_green_free(green_p);

  return s;

#else

  (void) r;
  (void) theta;
  (void) phi;
  (void) g;
  (void) w;
  
  dB[0] = dB[1] = dB[2] = 0.0;
  return 0;

#endif

} /* mfield_nonlinear_model_ext() */

/*
mfield_nonlinear_regularize_init()
  Construct the p-by-p regularization matrix Lambda.

                    p_core           p_crust         p_euler          p_fluxcal
Lambda = p_core    [ G^{(n)} x C  |              |               |                ]
         p_crust   [              |      0       |               |                ]
         p_euler   [              |              | G^{(2)} x I_3 |                ]
         p_fluxcal [              |              |               | G^{(2)} x I_9  ]

Here, "x" denotes Kronecker product.
G^{(n)} denotes the Gram matrix G_{ij} = 1 / dt * \int_{min(knots)}^{max(knots)} d^n/dt^n N_i(t) d^n/dt^n N_j(t) dt
C is a diagonal matrix given by:

C_{nm,n'm'} = (a/c)^{n+2} (n+1)^2 / (2n + 1} delta_{mm'} delta_{nn'}

Inputs: w - workspace

Notes:
1) on output, w->L contains the Cholesky factor of Lambda; only the lower triangle is stored
since the matrix is symmetric

2) on output, w->Lambda contains Lambda (L L^T) stored in the lower triangle
*/

static int
mfield_nonlinear_regularize_init(mfield_workspace *w)
{
  int s = 0;

  gsl_spmatrix_set_zero(w->L);
  gsl_spmatrix_set_zero(w->Lambda);

  if (w->p_core > 0)
    {
      s += mfield_nonlinear_regularize_core(w->lambda_0, 0, w);
      s += mfield_nonlinear_regularize_core(w->lambda_1, 1, w);
      s += mfield_nonlinear_regularize_core(w->lambda_2, 2, w);
      s += mfield_nonlinear_regularize_core(w->lambda_3, 3, w);
    }

  if (w->p_fluxcal > 0)
    {
      size_t n;

      for (n = 0; n < w->nsat; ++n)
        {
          const size_t offset = w->fluxcal_offset + w->offset_fluxcal[n];
          magdata *mptr = mfield_data_ptr(n, w->data_workspace_p);
          gsl_bspline2_workspace *fluxcal_spline_p = w->fluxcal_spline_workspace_p[CIDX2(n, w->nsat, 0, w->max_threads)];
          size_t ncontrol, order;
          double t0, t1;
          gsl_matrix * G_fluxcal, *G_fluxcal_orig;
          size_t i, j;

          if (!(mptr->global_flags & MAGDATA_GLOBFLG_FLUXCAL))
            continue;

          ncontrol = gsl_bspline2_ncontrol(fluxcal_spline_p);
          order = gsl_bspline2_order(fluxcal_spline_p);
          t0 = gsl_vector_get(fluxcal_spline_p->knots, 0);
          t1 = gsl_vector_get(fluxcal_spline_p->knots, fluxcal_spline_p->knots->size - 1);

          G_fluxcal = gsl_matrix_alloc(ncontrol, order);
          G_fluxcal_orig = gsl_matrix_alloc(ncontrol, order);

          /* construct G_{ij} = 1/dt * \int N^{(2)}_i(t) N^{(2)}_j(t) dt for fluxcal splines */
          gsl_bspline2_gram(2, G_fluxcal, fluxcal_spline_p);
          printsb_octave(G_fluxcal, "G_fluxcal1");
          gsl_matrix_scale(G_fluxcal, 1.0 / (t1 - t0));

          printsb_octave(G_fluxcal, "G_fluxcal");

          gsl_matrix_memcpy(G_fluxcal_orig, G_fluxcal);

          {
            /* XXX: need to add epsilon to diagonal of G^{(2)} so LLT decomposition will work */
            gsl_vector_view diag = gsl_matrix_column(G_fluxcal, 0);
            gsl_vector_add_constant(&diag.vector, 1.0e-10);
            gsl_linalg_cholesky_band_decomp(G_fluxcal);
          }

          /*
           * construct L(1:p_fluxcal,1:p_fluxcal) = G_L^{(2)} x I_9,
           * the Kronecker product of the Cholesky factors of G^{(2)} and I_9
           */
          for (j = 0; j < ncontrol; ++j)
            {
              for (i = 0; i < order && i + j < ncontrol; ++i)
                {
                  size_t row = i + j;
                  double Lij = gsl_matrix_get(G_fluxcal, j, i); /* Cholesky factor of G^{(3)} */
                  double Gij = gsl_matrix_get(G_fluxcal_orig, j, i);
                  size_t k;

                  for (k = FLUXCAL_IDX_SX; k <= FLUXCAL_IDX_SZ; ++k)
                    {
                      gsl_spmatrix_set(w->L, offset + row * FLUXCAL_P + k, offset + j * FLUXCAL_P + k, w->lambda_s * Lij);
                      gsl_spmatrix_set(w->Lambda, offset + row * FLUXCAL_P + k, offset + j * FLUXCAL_P + k, w->lambda_s * w->lambda_s * Gij);
                    }

                  for (k = FLUXCAL_IDX_OX; k <= FLUXCAL_IDX_OZ; ++k)
                    {
                      gsl_spmatrix_set(w->L, offset + row * FLUXCAL_P + k, offset + j * FLUXCAL_P + k, w->lambda_o * Lij);
                      gsl_spmatrix_set(w->Lambda, offset + row * FLUXCAL_P + k, offset + j * FLUXCAL_P + k, w->lambda_o * w->lambda_o * Gij);
                    }

                  for (k = FLUXCAL_IDX_U1; k <= FLUXCAL_IDX_U3; ++k)
                    {
                      gsl_spmatrix_set(w->L, offset + row * FLUXCAL_P + k, offset + j * FLUXCAL_P + k, w->lambda_u * Lij);
                      gsl_spmatrix_set(w->Lambda, offset + row * FLUXCAL_P + k, offset + j * FLUXCAL_P + k, w->lambda_u * w->lambda_u * Gij);
                    }
                }
            }

          gsl_matrix_free(G_fluxcal);
          gsl_matrix_free(G_fluxcal_orig);
        }
    }

#if 0 /*XXX*/
  {
    gsl_matrix *Lambda = gsl_matrix_alloc(w->p, w->p);
    gsl_matrix *L = gsl_matrix_alloc(w->p, w->p);
    gsl_spmatrix_sp2d(Lambda, w->Lambda);
    gsl_spmatrix_sp2d(L, w->L);
    print_octave(L, "L");
    printsym_octave(Lambda, "Lambda");
    exit(1);
  }
#endif

  return s;
}

/*
mfield_nonlinear_regularize_core()
  Regularize core field by minimizing:

< | d^n / dt^n B_r |^2 >

averaged over full time interval and over the surface
of the core-mantle boundary

Inputs: lambda - regularization parameter
        nderiv - derivative to minimize
        w      - workspace

Notes:
1) On output,
*/

static int
mfield_nonlinear_regularize_core(const double lambda, const size_t nderiv, mfield_workspace *w)
{
  gsl_bspline2_workspace * gauss_spline_p = w->gauss_spline_workspace_p[0];
  const size_t order = gsl_bspline2_order(gauss_spline_p);

  if (nderiv >= order || lambda == 0.0)
    {
      /* nothing to do */
      return GSL_SUCCESS;
    }
  else
    {
      int s = 0;
      const double lambda_sq = lambda * lambda;
      const size_t nmax = (nderiv == 0) ? w->nmax : w->nmax_core;
      const size_t nnm = (nderiv == 0) ? w->nnm_tot : w->nnm_core;
      const double t0 = gsl_vector_get(gauss_spline_p->knots, 0);
      const double t1 = gsl_vector_get(gauss_spline_p->knots, gauss_spline_p->knots->size - 1);
      const size_t ncontrol = gsl_bspline2_ncontrol(gauss_spline_p);
      const double c = 3485.0;       /* Earth core radius */
      const double a = w->params.R;  /* Earth surface radius */
      const double ratio = a / c;
      double rterm = ratio * ratio * ratio; /* (a/c)^{n+2} */
      size_t n, i, j;
      gsl_vector * c_core = gsl_vector_calloc(nnm);
      gsl_matrix * G_core = gsl_matrix_alloc(ncontrol, order);
      gsl_matrix * G_core_orig = gsl_matrix_alloc(ncontrol, order);

      /*
       * construct vector c_core = diag(C) for minimizing
       *
       * < d^n/dt^n B_r > averaged over CMB. c_core contains the spatial part:
       *
       * (a/c)^{n+2} (n+1) / sqrt(2n + 1)
       */

      for (n = 1; n <= nmax; ++n)
        {
          int ni = (int) n;
          int m;
          double term = (n + 1.0) / sqrt(2.0*n + 1.0) * rterm;

          for (m = -ni; m <= ni; ++m)
            {
              size_t cidx = mfield_coeff_nmidx(n, m);
              gsl_vector_set(c_core, cidx, term);
            }

          rterm *= ratio;
        }

      /* construct G_{ij} = 1/dt * \int N^{(nderiv)}_i(t) N^{(nderiv)}_j(t) dt for Gauss coefficient splines */
      gsl_bspline2_gram(nderiv, G_core, gauss_spline_p);
      gsl_matrix_scale(G_core, 1.0 / (t1 - t0));

      gsl_matrix_memcpy(G_core_orig, G_core);

      if (nderiv > 0)
        {
          /* XXX: need to add epsilon to diagonal of G^{(nderiv)} so LLT decomposition
           * will work */
          gsl_vector_view diag = gsl_matrix_column(G_core, 0);
          gsl_vector_add_constant(&diag.vector, 1.0e-10);
          gsl_linalg_cholesky_band_decomp(G_core);
        }

      /*
       * construct L(1:p_core,1:p_core) = G_L^{(nderiv)} x C_L,
       * the Kronecker product of the Cholesky factors of G^{(nderiv)} and C
       */
      for (j = 0; j < ncontrol; ++j)
        {
          for (i = 0; i < order && i + j < ncontrol; ++i)
            {
              size_t row = i + j;
              double Lij = gsl_matrix_get(G_core, j, i); /* Cholesky factor of G^{(nderiv)} */
              double Gij = gsl_matrix_get(G_core_orig, j, i);
              size_t k;

              for (k = 0; k < w->nnm_core; ++k)
                {
                  double ck = gsl_vector_get(c_core, k);
                  size_t idx1 = row * w->nnm_core + k;
                  size_t idx2 = j * w->nnm_core + k;
                  double * Lptr = gsl_spmatrix_ptr(w->L, idx1, idx2);
                  double * Lamptr = gsl_spmatrix_ptr(w->Lambda, idx1, idx2);

                  if (Lptr)
                    *Lptr += lambda * Lij * ck;
                  else
                    gsl_spmatrix_set(w->L, idx1, idx2, lambda * Lij * ck);

                  if (Lamptr)
                    *Lamptr += lambda_sq * Gij * ck * ck;
                  else
                    gsl_spmatrix_set(w->Lambda, idx1, idx2, lambda_sq * Gij * ck * ck);
                }
            }
        }

      if (nderiv == 0)
        {
          size_t k;

          /* regularize crustal part */
          for (k = w->nnm_core; k < nnm; ++k)
            {
              double ck = gsl_vector_get(c_core, k);
              size_t idx = w->p_core + k - w->nnm_core;

              gsl_spmatrix_set(w->L, idx, idx, lambda * ck);
              gsl_spmatrix_set(w->Lambda, idx, idx, lambda_sq * ck * ck);
            }
        }

      gsl_vector_free(c_core);
      gsl_matrix_free(G_core);
      gsl_matrix_free(G_core_orig);

      return s;
    }
}

static int
mfield_robust_print_stat(const char *str, const double value, const gsl_rstat_workspace *rstat_p)
{
  const size_t n = gsl_rstat_n(rstat_p);

  if (n > 0)
    {
      const double mean = gsl_rstat_mean(rstat_p);
      fprintf(stderr, "\t %18s = %.2f [nT], Robust weight mean = %.4f\n", str, value, mean);
    }

  return 0;
}

/*
mfield_robust_weights()
  Compute robust weights given a vector of residuals

Inputs: f   - vector of unweighted residuals, length nres
        wts - (output) robust weights, length nres
        w   - workspace
*/

static int
mfield_robust_weights(const gsl_vector * f, gsl_vector * wts, mfield_workspace * w)
{
  int s = 0;
  const double tune = w->robust_workspace_p->tune;
  const double qdlat_cutoff = w->params.qdlat_fit_cutoff; /* cutoff latitude for high/low statistics */
  size_t i, j;
  size_t idx = 0;
  gsl_rstat_workspace **rstat_x = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_y = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_z = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_f = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_dx_ns = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_dy_ns = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_low_dz_ns = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_high_dz_ns = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_dx_ew = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_dy_ew = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_low_dz_ew = malloc(w->nsat * sizeof(gsl_rstat_workspace *));
  gsl_rstat_workspace **rstat_high_dz_ew = malloc(w->nsat * sizeof(gsl_rstat_workspace *));

  gsl_rstat_workspace *rstat_obs_x = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_obs_y = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_obs_z = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_dxdt = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_dydt = gsl_rstat_alloc();
  gsl_rstat_workspace *rstat_dzdt = gsl_rstat_alloc();
  double mean_OBS_X, mean_OBS_Y, mean_OBS_Z;
  double sigma_OBS_X, sigma_OBS_Y, sigma_OBS_Z;
  double mean_DXDT, mean_DYDT, mean_DZDT;
  double sigma_DXDT, sigma_DYDT, sigma_DZDT;

  for (i = 0; i < w->nsat; ++i)
    {
      rstat_x[i] = gsl_rstat_alloc();
      rstat_y[i] = gsl_rstat_alloc();
      rstat_z[i] = gsl_rstat_alloc();
      rstat_f[i] = gsl_rstat_alloc();
      rstat_dx_ns[i] = gsl_rstat_alloc();
      rstat_dy_ns[i] = gsl_rstat_alloc();
      rstat_low_dz_ns[i] = gsl_rstat_alloc();
      rstat_high_dz_ns[i] = gsl_rstat_alloc();
      rstat_dx_ew[i] = gsl_rstat_alloc();
      rstat_dy_ew[i] = gsl_rstat_alloc();
      rstat_low_dz_ew[i] = gsl_rstat_alloc();
      rstat_high_dz_ew[i] = gsl_rstat_alloc();
    }

  fprintf(stderr, "\n");

  /*
   * first loop through the residuals and compute statistics for each residual type
   * (X,Y,Z,F,DX,DY,DZ,dX/dt,dY/dt,dZ/dt)
   */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          /* check if data point is discarded due to time interval */
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_ExistX(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
                    gsl_rstat_add(fi, rstat_obs_x);
                  else
                    gsl_rstat_add(fi, rstat_x[i]);
                }
            }

          if (MAGDATA_ExistY(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
                    gsl_rstat_add(fi, rstat_obs_y);
                  else
                    gsl_rstat_add(fi, rstat_y[i]);
                }
            }

          if (MAGDATA_ExistZ(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);

              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
                    gsl_rstat_add(fi, rstat_obs_z);
                  else
                    gsl_rstat_add(fi, rstat_z[i]);
                }
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx++);
              gsl_rstat_add(fi, rstat_f[i]);
            }

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              if (MAGDATA_ExistDXDT(mptr->flags[j]))
                {
                  double fi = gsl_vector_get(f, idx++);
                  gsl_rstat_add(fi, rstat_dxdt);
                }

              if (MAGDATA_ExistDYDT(mptr->flags[j]))
                {
                  double fi = gsl_vector_get(f, idx++);
                  gsl_rstat_add(fi, rstat_dydt);
                }

              if (MAGDATA_ExistDZDT(mptr->flags[j]))
                {
                  double fi = gsl_vector_get(f, idx++);
                  gsl_rstat_add(fi, rstat_dzdt);
                }
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

  assert(idx == w->nres);

  mean_OBS_X = gsl_rstat_mean(rstat_obs_x);
  mean_OBS_Y = gsl_rstat_mean(rstat_obs_y);
  mean_OBS_Z = gsl_rstat_mean(rstat_obs_z);
  sigma_OBS_X = gsl_rstat_sd(rstat_obs_x);
  sigma_OBS_Y = gsl_rstat_sd(rstat_obs_y);
  sigma_OBS_Z = gsl_rstat_sd(rstat_obs_z);
  gsl_rstat_reset(rstat_obs_x);
  gsl_rstat_reset(rstat_obs_y);
  gsl_rstat_reset(rstat_obs_z);

  mean_DXDT = gsl_rstat_mean(rstat_dxdt);
  mean_DYDT = gsl_rstat_mean(rstat_dydt);
  mean_DZDT = gsl_rstat_mean(rstat_dzdt);
  sigma_DXDT = gsl_rstat_sd(rstat_dxdt);
  sigma_DYDT = gsl_rstat_sd(rstat_dydt);
  sigma_DZDT = gsl_rstat_sd(rstat_dzdt);
  gsl_rstat_reset(rstat_dxdt);
  gsl_rstat_reset(rstat_dydt);
  gsl_rstat_reset(rstat_dzdt);

  /* loop through again and compute robust weights and the mean of the weights */

  idx = 0;
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      const double alpha = 1.0; /* constant to multiply sigma so that mean(weights) = 0.95 */
      double sigma_X = alpha * gsl_rstat_sd(rstat_x[i]);
      double sigma_Y = alpha * gsl_rstat_sd(rstat_y[i]);
      double sigma_Z = alpha * gsl_rstat_sd(rstat_z[i]);
      double sigma_F = alpha * gsl_rstat_sd(rstat_f[i]);
      double sigma_DX_NS = alpha * gsl_rstat_sd(rstat_dx_ns[i]);
      double sigma_DY_NS = alpha * gsl_rstat_sd(rstat_dy_ns[i]);
      double sigma_low_DZ_NS = alpha * gsl_rstat_sd(rstat_low_dz_ns[i]);
      double sigma_high_DZ_NS = alpha * gsl_rstat_sd(rstat_high_dz_ns[i]);
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
              double wi;

              if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
                {
                  wi = bisquare(fi / (tune * sigma_OBS_X));
                  gsl_rstat_add(wi, rstat_obs_x);
                }
              else
                {
                  wi = bisquare(fi / (tune * sigma_X));
                  gsl_rstat_add(wi, rstat_x[i]);
                }

              gsl_vector_set(wts, idx++, wi);
            }

          if (MAGDATA_ExistY(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi;

              if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
                {
                  wi = bisquare(fi / (tune * sigma_OBS_Y));
                  gsl_rstat_add(wi, rstat_obs_y);
                }
              else
                {
                  wi = bisquare(fi / (tune * sigma_Y));
                  gsl_rstat_add(wi, rstat_y[i]);
                }

              gsl_vector_set(wts, idx++, wi);
            }

          if (MAGDATA_ExistZ(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi;

              if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
                {
                  wi = bisquare(fi / (tune * sigma_OBS_Z));
                  gsl_rstat_add(wi, rstat_obs_z);
                }
              else
                {
                  wi = bisquare(fi / (tune * sigma_Z));
                  gsl_rstat_add(wi, rstat_z[i]);
                }

              gsl_vector_set(wts, idx++, wi);
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = bisquare(fi / (tune * sigma_F));
              gsl_vector_set(wts, idx++, wi);
              gsl_rstat_add(wi, rstat_f[i]);
            }

          if (MAGDATA_FitMF(mptr->flags[j]))
            {
              if (MAGDATA_ExistDXDT(mptr->flags[j]))
                {
                  double fi = gsl_vector_get(f, idx);
                  double wi = bisquare(fi / (tune * sigma_DXDT));
                  gsl_vector_set(wts, idx++, wi);
                  gsl_rstat_add(wi, rstat_dxdt);
                }

              if (MAGDATA_ExistDYDT(mptr->flags[j]))
                {
                  double fi = gsl_vector_get(f, idx);
                  double wi = bisquare(fi / (tune * sigma_DYDT));
                  gsl_vector_set(wts, idx++, wi);
                  gsl_rstat_add(wi, rstat_dydt);
                }

              if (MAGDATA_ExistDZDT(mptr->flags[j]))
                {
                  double fi = gsl_vector_get(f, idx);
                  double wi = bisquare(fi / (tune * sigma_DZDT));
                  gsl_vector_set(wts, idx++, wi);
                  gsl_rstat_add(wi, rstat_dzdt);
                }
            }

          if (MAGDATA_ExistDX_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = bisquare(fi / (tune * sigma_DX_NS));
              gsl_vector_set(wts, idx++, wi);
              gsl_rstat_add(wi, rstat_dx_ns[i]);
            }

          if (MAGDATA_ExistDY_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = bisquare(fi / (tune * sigma_DY_NS));
              gsl_vector_set(wts, idx++, wi);
              gsl_rstat_add(wi, rstat_dy_ns[i]);
            }

          if (MAGDATA_ExistDZ_NS(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi;
              
              if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                {
                  wi = bisquare(fi / (tune * sigma_low_DZ_NS));
                  gsl_rstat_add(wi, rstat_low_dz_ns[i]);
                }
              else
                {
                  wi = bisquare(fi / (tune * sigma_high_DZ_NS));
                  gsl_rstat_add(wi, rstat_high_dz_ns[i]);
                }

              gsl_vector_set(wts, idx++, wi);
            }

          if (MAGDATA_ExistDX_EW(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = bisquare(fi / (tune * sigma_DX_EW));
              gsl_vector_set(wts, idx++, wi);
              gsl_rstat_add(wi, rstat_dx_ew[i]);
            }

          if (MAGDATA_ExistDY_EW(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = bisquare(fi / (tune * sigma_DY_EW));
              gsl_vector_set(wts, idx++, wi);
              gsl_rstat_add(wi, rstat_dy_ew[i]);
            }

          if (MAGDATA_ExistDZ_EW(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi;

              if (fabs(mptr->qdlat[j]) <= qdlat_cutoff)
                {
                  wi = bisquare(fi / (tune * sigma_low_DZ_EW));
                  gsl_rstat_add(wi, rstat_low_dz_ew[i]);
                }
              else
                {
                  wi = bisquare(fi / (tune * sigma_high_DZ_EW));
                  gsl_rstat_add(wi, rstat_high_dz_ew[i]);
                }

              gsl_vector_set(wts, idx++, wi);
            }
        }

      if (mptr->global_flags & MAGDATA_GLOBFLG_SATELLITE)
        {
          fprintf(stderr, "\t === SATELLITE %zu (robust sigma) ===\n", i);

          mfield_robust_print_stat("sigma X", sigma_X, rstat_x[i]);
          mfield_robust_print_stat("sigma Y", sigma_Y, rstat_y[i]);
          mfield_robust_print_stat("sigma Z", sigma_Z, rstat_z[i]);
          mfield_robust_print_stat("sigma F", sigma_F, rstat_f[i]);

          mfield_robust_print_stat("sigma DX_NS", sigma_DX_NS, rstat_dx_ns[i]);
          mfield_robust_print_stat("sigma DY_NS", sigma_DY_NS, rstat_dy_ns[i]);
          mfield_robust_print_stat("sigma low DZ_NS", sigma_low_DZ_NS, rstat_low_dz_ns[i]);
          mfield_robust_print_stat("sigma high DZ_NS", sigma_high_DZ_NS, rstat_high_dz_ns[i]);

          mfield_robust_print_stat("sigma DX_EW", sigma_DX_EW, rstat_dx_ew[i]);
          mfield_robust_print_stat("sigma DY_EW", sigma_DY_EW, rstat_dy_ew[i]);
          mfield_robust_print_stat("sigma low DZ_EW", sigma_low_DZ_EW, rstat_low_dz_ew[i]);
          mfield_robust_print_stat("sigma high DZ_EW", sigma_high_DZ_EW, rstat_high_dz_ew[i]);
        }
    }

  fprintf(stderr, "\t === OBSERVATORY (robust sigma) ===\n");

  mfield_robust_print_stat("mean X", mean_OBS_X, rstat_obs_x);
  mfield_robust_print_stat("mean Y", mean_OBS_Y, rstat_obs_y);
  mfield_robust_print_stat("mean Z", mean_OBS_Z, rstat_obs_z);
  mfield_robust_print_stat("sigma X", sigma_OBS_X, rstat_obs_x);
  mfield_robust_print_stat("sigma Y", sigma_OBS_Y, rstat_obs_y);
  mfield_robust_print_stat("sigma Z", sigma_OBS_Z, rstat_obs_z);

  fprintf(stderr, "\t === OBSERVATORY SV (robust sigma) ===\n");

  mfield_robust_print_stat("mean dX/dt", mean_DXDT, rstat_dxdt);
  mfield_robust_print_stat("mean dY/dt", mean_DYDT, rstat_dydt);
  mfield_robust_print_stat("mean dZ/dt", mean_DZDT, rstat_dzdt);
  mfield_robust_print_stat("sigma dX/dt", sigma_DXDT, rstat_dxdt);
  mfield_robust_print_stat("sigma dY/dt", sigma_DYDT, rstat_dydt);
  mfield_robust_print_stat("sigma dZ/dt", sigma_DZDT, rstat_dzdt);

  assert(idx == w->nres);

  for (i = 0; i < w->nsat; ++i)
    {
      gsl_rstat_free(rstat_x[i]);
      gsl_rstat_free(rstat_y[i]);
      gsl_rstat_free(rstat_z[i]);
      gsl_rstat_free(rstat_f[i]);
      gsl_rstat_free(rstat_dx_ns[i]);
      gsl_rstat_free(rstat_dy_ns[i]);
      gsl_rstat_free(rstat_low_dz_ns[i]);
      gsl_rstat_free(rstat_high_dz_ns[i]);
      gsl_rstat_free(rstat_dx_ew[i]);
      gsl_rstat_free(rstat_dy_ew[i]);
      gsl_rstat_free(rstat_low_dz_ew[i]);
      gsl_rstat_free(rstat_high_dz_ew[i]);
    }

  gsl_rstat_free(rstat_obs_x);
  gsl_rstat_free(rstat_obs_y);
  gsl_rstat_free(rstat_obs_z);
  gsl_rstat_free(rstat_dxdt);
  gsl_rstat_free(rstat_dydt);
  gsl_rstat_free(rstat_dzdt);

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

static void
mfield_nonlinear_callback(const size_t iter, void *params,
                          const gsl_multifit_nlinear_workspace *multifit_p)
{
  mfield_workspace *w = (mfield_workspace *) params;
  gsl_vector *x = gsl_multifit_nlinear_position(multifit_p);
  gsl_vector *f = gsl_multifit_nlinear_residual(multifit_p);
  double avratio = gsl_multifit_nlinear_avratio(multifit_p);
  double rcond;

  /* print out state every 5 iterations */
  if (iter % 5 != 0 && iter != 1)
    return;

  fprintf(stderr, "iteration %zu (method: %s/%s):\n",
          iter,
          gsl_multifit_nlinear_name(multifit_p),
          gsl_multifit_nlinear_trs_name(multifit_p));

  if (w->nnm_core > 0)
    {
      const double t0 = epoch2year(w->data_workspace_p->t0_data);
      const double t1 = epoch2year(w->data_workspace_p->t1_data);
      const double t = 0.5 * (t0 + t1);
      double g10, g11, h11;
      double dg10, dg11, dh11;

      g10 = mfield_get_gnm(t, 1, 0, 0, x, w);
      dg10 = mfield_get_gnm(t, 1, 0, 1, x, w);

      g11 = mfield_get_gnm(t, 1, 1, 0, x, w);
      dg11 = mfield_get_gnm(t, 1, 1, 1, x, w);

      h11 = mfield_get_gnm(t, 1, -1, 0, x, w);
      dh11 = mfield_get_gnm(t, 1, -1, 1, x, w);

      fprintf(stderr, "\t %-10s %12.4f %12.4f %12.4f [nT]\n",
              "dipole:", g10, g11, h11);
      fprintf(stderr, "\t %-10s %12.4f %12.4f %12.4f [nT]\n",
              "SV dipole:",
              dg10,
              dg11,
              dh11);
    }

  if (w->params.fit_euler)
    {
      size_t i;

      for (i = 0; i < w->nsat; ++i)
        {
          magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
  
          if (mptr->n == 0)
            continue;

          if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
            {
              double t0 = w->data_workspace_p->t0[i];
              gsl_bspline2_workspace *euler_spline_p = w->euler_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)];
              size_t euler_idx = w->euler_offset + w->offset_euler[i];
              size_t ncontrol = gsl_bspline2_ncontrol(euler_spline_p);
              gsl_vector_const_view tmp = gsl_vector_const_subvector(x, euler_idx, EULER_P * ncontrol);
              gsl_matrix_const_view control_pts = gsl_matrix_const_view_vector(&tmp.vector, EULER_P, ncontrol);
              double euler_data[EULER_P];
              gsl_vector_view euler_params = gsl_vector_view_array(euler_data, EULER_P);
              char buf[32];

              if (mptr->euler_flags & EULER_FLG_ZYX)
                sprintf(buf, "ZYX");
              else if (mptr->euler_flags & EULER_FLG_ZYZ)
                sprintf(buf, "ZYZ");
              else
                *buf = '\0';

              gsl_bspline2_vector_eval(t0, &control_pts.matrix, &euler_params.vector, euler_spline_p);

              fprintf(stderr, "\t euler %zu: %12.4f %12.4f %12.4f [deg] [%s]\n",
                      i,
                      euler_data[0] * 180.0 / M_PI,
                      euler_data[1] * 180.0 / M_PI,
                      euler_data[2] * 180.0 / M_PI,
                      buf);
            }
        }
    }

  if (w->params.fit_fluxcal)
    {
      size_t i;

      for (i = 0; i < w->nsat; ++i)
        {
          magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
  
          if (mptr->n == 0)
            continue;

          if (mptr->global_flags & MAGDATA_GLOBFLG_FLUXCAL)
            {
              double t0 = w->data_workspace_p->t0[i];
              gsl_bspline2_workspace *fluxcal_spline_p = w->fluxcal_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)];
              size_t fluxcal_idx = w->fluxcal_offset + w->offset_fluxcal[i];
              size_t ncontrol = gsl_bspline2_ncontrol(fluxcal_spline_p);
              gsl_vector_const_view tmp = gsl_vector_const_subvector(x, fluxcal_idx, FLUXCAL_P * ncontrol);
              gsl_matrix_const_view control_pts = gsl_matrix_const_view_vector(&tmp.vector, FLUXCAL_P, ncontrol);
              double cal_data[FLUXCAL_P];
              gsl_vector_view cal_params = gsl_vector_view_array(cal_data, FLUXCAL_P);

              gsl_bspline2_vector_eval(t0, &control_pts.matrix, &cal_params.vector, fluxcal_spline_p);

              fprintf(stderr, "\t fluxcal %zu: S = %12.4f %12.4f %12.4f\n",
                      i,
                      cal_data[FLUXCAL_IDX_SX],
                      cal_data[FLUXCAL_IDX_SY],
                      cal_data[FLUXCAL_IDX_SZ]);
              fprintf(stderr, "\t            O = %12.4f %12.4f %12.4f [nT]\n",
                      cal_data[FLUXCAL_IDX_OX],
                      cal_data[FLUXCAL_IDX_OY],
                      cal_data[FLUXCAL_IDX_OZ]);
              fprintf(stderr, "\t            U = %12.4f %12.4f %12.4f [deg]\n",
                      cal_data[FLUXCAL_IDX_U1] * 180.0 / M_PI,
                      cal_data[FLUXCAL_IDX_U2] * 180.0 / M_PI,
                      cal_data[FLUXCAL_IDX_U3] * 180.0 / M_PI);
            }
        }
    }

  if (w->params.fit_cbias)
    {
      size_t i;

      for (i = 0; i < w->nsat; ++i)
        {
          magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

          if (!(mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY))
            continue;

          if (!strcasecmp(mptr->name, "MBO0"))
            {
              fprintf(stderr, "\t %-10s %12.4f %12.4f %12.4f [nT]\n",
                      "MBO bias:",
                      gsl_vector_get(x, w->bias_offset + w->bias_idx[i]),
                      gsl_vector_get(x, w->bias_offset + w->bias_idx[i] + 1),
                      gsl_vector_get(x, w->bias_offset + w->bias_idx[i] + 2));
            }
        }
    }

  fprintf(stderr, "\t |a|/|v|:    %12g\n", avratio);
  fprintf(stderr, "\t ||f(x)||:   %12g\n", gsl_blas_dnrm2(f));

  gsl_multifit_nlinear_rcond(&rcond, multifit_p);
  fprintf(stderr, "\t cond(J(x)): %12g\n", 1.0 / rcond);
}

static void
mfield_nonlinear_callback2(const size_t iter, void *params,
                           const gsl_multilarge_nlinear_workspace *multilarge_p)
{
  mfield_workspace *w = (mfield_workspace *) params;
  gsl_vector *x = gsl_multilarge_nlinear_position(multilarge_p);
  gsl_vector *dx = gsl_multilarge_nlinear_step(multilarge_p);
  gsl_vector *f = gsl_multilarge_nlinear_residual(multilarge_p);
  double avratio = gsl_multilarge_nlinear_avratio(multilarge_p);
  double rcond, normf, max_dx = 0.0;
  size_t i;

  /* print out state every 5 iterations */
  if (iter % 5 != 0 && iter != 1)
    return;

  fprintf(stderr, "iteration %zu (method: %s/%s):\n",
          iter,
          gsl_multilarge_nlinear_name(multilarge_p),
          gsl_multilarge_nlinear_trs_name(multilarge_p));

  if (w->nnm_core > 0)
    {
      const double t0 = epoch2year(w->data_workspace_p->t0_data);
      const double t1 = epoch2year(w->data_workspace_p->t1_data);
      const double t = 0.5 * (t0 + t1);
      double g10, g11, h11;
      double dg10, dg11, dh11;

      g10 = mfield_get_gnm(t, 1, 0, 0, x, w);
      dg10 = mfield_get_gnm(t, 1, 0, 1, x, w);

      g11 = mfield_get_gnm(t, 1, 1, 0, x, w);
      dg11 = mfield_get_gnm(t, 1, 1, 1, x, w);

      h11 = mfield_get_gnm(t, 1, -1, 0, x, w);
      dh11 = mfield_get_gnm(t, 1, -1, 1, x, w);

      fprintf(stderr, "\t %-10s %12.4f %12.4f %12.4f [nT]\n",
              "dipole:", g10, g11, h11);
      fprintf(stderr, "\t %-10s %12.4f %12.4f %12.4f [nT]\n",
              "SV dipole:",
              dg10,
              dg11,
              dh11);
    }

  if (w->params.fit_euler)
    {
      size_t i;

      for (i = 0; i < w->nsat; ++i)
        {
          magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
  
          if (mptr->n == 0)
            continue;

          if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
            {
              double t0 = w->data_workspace_p->t0[i];
              gsl_bspline2_workspace *euler_spline_p = w->euler_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)];
              size_t euler_idx = w->euler_offset + w->offset_euler[i];
              size_t ncontrol = gsl_bspline2_ncontrol(euler_spline_p);
              gsl_vector_const_view tmp = gsl_vector_const_subvector(x, euler_idx, EULER_P * ncontrol);
              gsl_matrix_const_view control_pts = gsl_matrix_const_view_vector(&tmp.vector, EULER_P, ncontrol);
              double euler_data[EULER_P];
              gsl_vector_view euler_params = gsl_vector_view_array(euler_data, EULER_P);
              char buf[32];

              if (mptr->euler_flags & EULER_FLG_ZYX)
                sprintf(buf, "ZYX");
              else if (mptr->euler_flags & EULER_FLG_ZYZ)
                sprintf(buf, "ZYZ");
              else
                *buf = '\0';

              gsl_bspline2_vector_eval(t0, &control_pts.matrix, &euler_params.vector, euler_spline_p);

              fprintf(stderr, "\t euler %zu: %12.4f %12.4f %12.4f [deg] [%s]\n",
                      i,
                      euler_data[0] * 180.0 / M_PI,
                      euler_data[1] * 180.0 / M_PI,
                      euler_data[2] * 180.0 / M_PI,
                      buf);
            }
        }
    }

  if (w->params.fit_fluxcal)
    {
      size_t i;

      for (i = 0; i < w->nsat; ++i)
        {
          magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
  
          if (mptr->n == 0)
            continue;

          if (mptr->global_flags & MAGDATA_GLOBFLG_FLUXCAL)
            {
              double t0 = w->data_workspace_p->t0[i];
              gsl_bspline2_workspace *fluxcal_spline_p = w->fluxcal_spline_workspace_p[CIDX2(i, w->nsat, 0, w->max_threads)];
              size_t fluxcal_idx = w->fluxcal_offset + w->offset_fluxcal[i];
              size_t ncontrol = gsl_bspline2_ncontrol(fluxcal_spline_p);
              gsl_vector_const_view tmp = gsl_vector_const_subvector(x, fluxcal_idx, FLUXCAL_P * ncontrol);
              gsl_matrix_const_view control_pts = gsl_matrix_const_view_vector(&tmp.vector, FLUXCAL_P, ncontrol);
              double cal_data[FLUXCAL_P];
              gsl_vector_view cal_params = gsl_vector_view_array(cal_data, FLUXCAL_P);

              gsl_bspline2_vector_eval(t0, &control_pts.matrix, &cal_params.vector, fluxcal_spline_p);

              fprintf(stderr, "\t fluxcal %zu: S = %12.4f %12.4f %12.4f\n",
                      i,
                      cal_data[FLUXCAL_IDX_SX],
                      cal_data[FLUXCAL_IDX_SY],
                      cal_data[FLUXCAL_IDX_SZ]);
              fprintf(stderr, "\t            O = %12.4f %12.4f %12.4f [nT]\n",
                      cal_data[FLUXCAL_IDX_OX],
                      cal_data[FLUXCAL_IDX_OY],
                      cal_data[FLUXCAL_IDX_OZ]);
              fprintf(stderr, "\t            U = %12.4f %12.4f %12.4f [deg]\n",
                      cal_data[FLUXCAL_IDX_U1] * 180.0 / M_PI,
                      cal_data[FLUXCAL_IDX_U2] * 180.0 / M_PI,
                      cal_data[FLUXCAL_IDX_U3] * 180.0 / M_PI);
            }
        }
    }

  /* compute max | dx_i / x_i | */
  for (i = 0; i < w->p; ++i)
    {
      double dxi = gsl_vector_get(dx, i);
      double xi = gsl_vector_get(x, i);

      if (fabs(xi) < GSL_DBL_EPSILON)
        continue;

      max_dx = GSL_MAX(max_dx, fabs(dxi / xi));
    }

  normf = gsl_blas_dnrm2(f);

  fprintf(stderr, "\t max |dx/x|:  %12g\n", max_dx);
  fprintf(stderr, "\t |a|/|v|:     %12g\n", avratio);
  fprintf(stderr, "\t |f(x)|:      %12g\n", normf);

  gsl_multilarge_nlinear_rcond(&rcond, multilarge_p);
  fprintf(stderr, "\t cond(J(x)):  %12g\n", 1.0 / rcond);
}

static int
gsl_linalg_symband_unpack(const gsl_matrix * AB, gsl_matrix * A)
{
  const size_t N = AB->size1;

  if (A->size1 != A->size2)
    {
      GSL_ERROR ("output matrix must be square", GSL_ENOTSQR);
    }
  else if (A->size1 != N)
    {
      GSL_ERROR ("matrix dimensions do not match", GSL_EBADLEN);
    }
  else
    {
      const size_t ndiag = AB->size2;
      size_t i;

      gsl_matrix_set_zero(A);

      for (i = 0; i < ndiag; ++i)
        {
          gsl_vector_const_view f = gsl_matrix_const_subcolumn(AB, i, 0, N - i);
          gsl_vector_view g = gsl_matrix_subdiagonal(A, i);
          gsl_vector_memcpy(&g.vector, &f.vector);
        }

      gsl_matrix_transpose_tricpy(CblasLower, CblasUnit, A, A);

      return GSL_SUCCESS;
    }
}
