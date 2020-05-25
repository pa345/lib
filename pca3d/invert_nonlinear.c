/* least squares solver to use,
 *
 * FDF_SOLVER = 0   GSL multifit solver / LM
 * FDF_SOLVER = 1   GSL multilarge solver / LM
 */
#define FDF_SOLVER     1

#include <gsl/gsl_integration.h>
#include <mainlib/ml_spatwt.h>

static int invert_init_nonlinear(invert_workspace *w);
static int invert_nonlinear_regularize_init(invert_workspace *w);
static void invert_nonlinear_callback(const size_t iter, void *params,
                                      const gsl_multifit_nlinear_workspace *multifit_p);
static int invert_robust_weights(const gsl_vector * f, gsl_vector * wts, invert_workspace * w);
static double huber(const double x);
static double bisquare(const double x);

#include "invert_multifit.c"
#include "invert_multilarge.c"

/*
invert_calc_nonlinear()
  Solve linear least squares system, using previously stored
satellite data

Inputs: c    - (input/output)
               on input, initial guess for coefficient vector
               on output, final coefficients
               units of nT, nT/year, nT/year^2
        w    - workspace

Notes:
1) On input, w->data_workspace_p must be filled with satellite data
2a) invert_init() must be called first to initialize various parameters,
    including weights
2b) this includes invert_init_nonlinear(), called from invert_init()
3) on output, coefficients are stored in w->c with following units:

Return: GSL_SUCCESS if converged, GSL_CONTINUE if not
*/

int
invert_calc_nonlinear(gsl_vector *c, invert_workspace *w)
{
  int s = 0;
  const invert_parameters *params = &(w->params);

  /* compute robust weights with coefficients from previous iteration */
  if (w->niter > 0)
    {
      gsl_vector_view wv = gsl_vector_subvector(w->wts_final, 0, w->nres);

      /* compute residuals f = Y_model - y_data with previous coefficients */
      invert_calc_f(c, w, w->fvec);

      /* compute robust weights */
      fprintf(stderr, "invert_calc_nonlinear: computing robust weights...");
      invert_robust_weights(w->fvec, w->wts_robust, w);
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

  s = invert_calc_nonlinear_multifit(c, w);

#elif FDF_SOLVER == 1 /* GSL multilarge */

  s = invert_calc_nonlinear_multilarge(c, w);

#endif

  gsl_vector_memcpy(c, w->c);

  w->niter++;

  return s;
}

gsl_vector *
invert_residual(const gsl_vector *c, invert_workspace *w)
{
#if FDF_SOLVER == 0
  gsl_vector *f = gsl_multifit_nlinear_residual(w->multifit_nlinear_p);
#elif FDF_SOLVER == 1
  gsl_vector *f = gsl_multilarge_nlinear_residual(w->nlinear_workspace_p);
#endif

  invert_calc_f(c, w, f);

  return f;
}

/*
invert_init_nonlinear()
  This function is called from invert_init() to count
total number of residuals and allocate nonlinear least
squares workspaces

Notes:
1) weight_calc() must be called prior to this function to
compute spatial weights
*/

static int
invert_init_nonlinear(invert_workspace *w)
{
  int s = 0;
  const invert_parameters *params = &(w->params);
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
      magdata *mptr = invert_data_ptr(i, w->data_workspace_p);

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
  if (nres_B[3] + nres_dB_ns[3] + nres_dB_ew[3] + nres_dBdt[3] == 0)
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

  fprintf(stderr, "invert_init_nonlinear: %zu distinct data points\n", ndata);
  fprintf(stderr, "invert_init_nonlinear: %zu scalar residuals\n", nres_B[3]);
  fprintf(stderr, "invert_init_nonlinear: %zu X residuals\n", nres_B[0]);
  fprintf(stderr, "invert_init_nonlinear: %zu Y residuals\n", nres_B[1]);
  fprintf(stderr, "invert_init_nonlinear: %zu Z residuals\n", nres_B[2]);
  fprintf(stderr, "invert_init_nonlinear: %zu d/dt X residuals\n", nres_dBdt[0]);
  fprintf(stderr, "invert_init_nonlinear: %zu d/dt Y residuals\n", nres_dBdt[1]);
  fprintf(stderr, "invert_init_nonlinear: %zu d/dt Z residuals\n", nres_dBdt[2]);
  fprintf(stderr, "invert_init_nonlinear: %zu dX_ns residuals\n", nres_dB_ns[0]);
  fprintf(stderr, "invert_init_nonlinear: %zu dY_ns residuals\n", nres_dB_ns[1]);
  fprintf(stderr, "invert_init_nonlinear: %zu dZ_ns residuals\n", nres_dB_ns[2]);
  fprintf(stderr, "invert_init_nonlinear: %zu dF_ns residuals\n", nres_dB_ns[3]);
  fprintf(stderr, "invert_init_nonlinear: %zu dX_ew residuals\n", nres_dB_ew[0]);
  fprintf(stderr, "invert_init_nonlinear: %zu dY_ew residuals\n", nres_dB_ew[1]);
  fprintf(stderr, "invert_init_nonlinear: %zu dZ_ew residuals\n", nres_dB_ew[2]);
  fprintf(stderr, "invert_init_nonlinear: %zu dF_ew residuals\n", nres_dB_ew[3]);
  fprintf(stderr, "invert_init_nonlinear: %zu data residuals\n", w->nres);
  fprintf(stderr, "invert_init_nonlinear: %zu total residuals (including regularization terms)\n", w->nres_tot);
  fprintf(stderr, "invert_init_nonlinear: %zu total parameters\n", p);

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

  w->multilarge_linear_p = gsl_multilarge_linear_alloc(gsl_multilarge_linear_tsqr, p);

#endif

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
#if 1 /*XXX*/
  {
    double weightfac[INVERT_DATA_IDX_END];

    weightfac[INVERT_DATA_IDX_X] = w->params.weight_X;
    weightfac[INVERT_DATA_IDX_Y] = w->params.weight_Y;
    weightfac[INVERT_DATA_IDX_Z] = w->params.weight_Z;
    weightfac[INVERT_DATA_IDX_F] = w->params.weight_F;

    invert_data_weights(w->wts_spatial, weightfac, w->data_workspace_p);
  }
#else
  {
    size_t idx = 0;
    size_t j;

    fprintf(stderr, "invert_init_nonlinear: calculating spatial weights...");

    for (i = 0; i < w->nsat; ++i)
      {
        magdata *mptr = invert_data_ptr(i, w->data_workspace_p);
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
#endif

  /* precompute regularization matrix */
  if (params->regularize && !params->synth_data)
    {
      fprintf(stderr, "invert_init_nonlinear: calculating regularization matrix...");
      invert_nonlinear_regularize_init(w);
      fprintf(stderr, "done\n");
    }
  else
    {
      gsl_spmatrix_set_zero(w->L);
      gsl_spmatrix_set_zero(w->Lambda);
    }

  return s;
} /* invert_init_nonlinear() */

/*
invert_nonlinear_regularize_init()
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
invert_nonlinear_regularize_init(invert_workspace *w)
{
  int s = 0;

  gsl_spmatrix_set_zero(w->L);
  gsl_spmatrix_set_zero(w->Lambda);

  return s;
}

static int
invert_robust_print_stat(const char *str, const double value, const gsl_rstat_workspace *rstat_p)
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
invert_robust_weights()
  Compute robust weights given a vector of residuals

Inputs: f   - vector of unweighted residuals, length nres
        wts - (output) robust weights, length nres
        w   - workspace
*/

static int
invert_robust_weights(const gsl_vector * f, gsl_vector * wts, invert_workspace * w)
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
      magdata *mptr = invert_data_ptr(i, w->data_workspace_p);

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
      magdata *mptr = invert_data_ptr(i, w->data_workspace_p);
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

          invert_robust_print_stat("sigma X", sigma_X, rstat_x[i]);
          invert_robust_print_stat("sigma Y", sigma_Y, rstat_y[i]);
          invert_robust_print_stat("sigma Z", sigma_Z, rstat_z[i]);
          invert_robust_print_stat("sigma F", sigma_F, rstat_f[i]);

          invert_robust_print_stat("sigma DX_NS", sigma_DX_NS, rstat_dx_ns[i]);
          invert_robust_print_stat("sigma DY_NS", sigma_DY_NS, rstat_dy_ns[i]);
          invert_robust_print_stat("sigma low DZ_NS", sigma_low_DZ_NS, rstat_low_dz_ns[i]);
          invert_robust_print_stat("sigma high DZ_NS", sigma_high_DZ_NS, rstat_high_dz_ns[i]);

          invert_robust_print_stat("sigma DX_EW", sigma_DX_EW, rstat_dx_ew[i]);
          invert_robust_print_stat("sigma DY_EW", sigma_DY_EW, rstat_dy_ew[i]);
          invert_robust_print_stat("sigma low DZ_EW", sigma_low_DZ_EW, rstat_low_dz_ew[i]);
          invert_robust_print_stat("sigma high DZ_EW", sigma_high_DZ_EW, rstat_high_dz_ew[i]);
        }
    }

  fprintf(stderr, "\t === OBSERVATORY (robust sigma) ===\n");

  invert_robust_print_stat("mean X", mean_OBS_X, rstat_obs_x);
  invert_robust_print_stat("mean Y", mean_OBS_Y, rstat_obs_y);
  invert_robust_print_stat("mean Z", mean_OBS_Z, rstat_obs_z);
  invert_robust_print_stat("sigma X", sigma_OBS_X, rstat_obs_x);
  invert_robust_print_stat("sigma Y", sigma_OBS_Y, rstat_obs_y);
  invert_robust_print_stat("sigma Z", sigma_OBS_Z, rstat_obs_z);

  fprintf(stderr, "\t === OBSERVATORY SV (robust sigma) ===\n");

  invert_robust_print_stat("mean dX/dt", mean_DXDT, rstat_dxdt);
  invert_robust_print_stat("mean dY/dt", mean_DYDT, rstat_dydt);
  invert_robust_print_stat("mean dZ/dt", mean_DZDT, rstat_dzdt);
  invert_robust_print_stat("sigma dX/dt", sigma_DXDT, rstat_dxdt);
  invert_robust_print_stat("sigma dY/dt", sigma_DYDT, rstat_dydt);
  invert_robust_print_stat("sigma dZ/dt", sigma_DZDT, rstat_dzdt);

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
invert_nonlinear_callback(const size_t iter, void *params,
                          const gsl_multifit_nlinear_workspace *multifit_p)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(multifit_p);
  double avratio = gsl_multifit_nlinear_avratio(multifit_p);
  double rcond;

  (void) params;

  /* print out state every 5 iterations */
  if (iter % 5 != 0 && iter != 1)
    return;

  fprintf(stderr, "iteration %zu (method: %s/%s):\n",
          iter,
          gsl_multifit_nlinear_name(multifit_p),
          gsl_multifit_nlinear_trs_name(multifit_p));

  fprintf(stderr, "\t |a|/|v|:    %12g\n", avratio);
  fprintf(stderr, "\t ||f(x)||:   %12g\n", gsl_blas_dnrm2(f));

  gsl_multifit_nlinear_rcond(&rcond, multifit_p);
  fprintf(stderr, "\t cond(J(x)): %12g\n", 1.0 / rcond);
}
