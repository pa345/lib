#define OLD_FDF     1
#define DEBUG       1

typedef struct
{
  mfield_workspace *w;
} mfield_nonlinear_params;

static int mfield_init_nonlinear(mfield_workspace *w);
static int mfield_calc_nonlinear_multilarge(const gsl_vector *c, mfield_workspace *w);
static int mfield_calc_Wf(const gsl_vector *x, void *params, gsl_vector *f);
static int mfield_calc_df2(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                           void *params, gsl_vector * v, gsl_matrix * JTJ);
static inline int mfield_jacobian_JTu(const double t, const size_t flags, const double weight,
                                      const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, const size_t extidx,
                                      const double dB_ext, const size_t euler_idx, const double B_nec_alpha,
                                      const double B_nec_beta, const double B_nec_gamma,
                                      gsl_vector *JTy, const mfield_workspace *w);
static inline int mfield_jacobian_grad_JTu(const double t, const double t_grad, const size_t flags, const double weight,
                                           const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, gsl_vector * dB_int_grad,
                                           gsl_vector *JTu, const mfield_workspace *w);
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
static int mfield_nonlinear_vector_precompute(const gsl_vector *weights, mfield_workspace *w);
static int mfield_vector_green(const double t, const double weight, const gsl_vector *g,
                               gsl_vector *G, mfield_workspace *w);
static int mfield_vector_green_grad(const double t, const double t_grad, const double weight, const gsl_vector *g,
                                    const gsl_vector *g_grad, gsl_vector *G, mfield_workspace *w);
static int mfield_nonlinear_model_ext(const double r, const double theta,
                                      const double phi, const gsl_vector *g,
                                      double dB[3], const mfield_workspace *w);
static int mfield_nonlinear_histogram(const gsl_vector *c,
                                      mfield_workspace *w);
static int mfield_nonlinear_regularize(gsl_vector *diag,
                                       mfield_workspace *w);
static void mfield_nonlinear_callback(const size_t iter, void *params,
                                      const gsl_multifit_nlinear_workspace *multifit_p);
static void mfield_nonlinear_callback2(const size_t iter, void *params,
                                       const gsl_multilarge_nlinear_workspace *multifit_p);
static int mfield_robust_weights(const gsl_vector * f, gsl_vector * wts, mfield_workspace * w);
static double huber(const double x);
static double bisquare(const double x);
static int mfield_nonlinear_alloc_multilarge(const gsl_multilarge_nlinear_trs * trs, mfield_workspace * w);

#include "mfield_multifit.c"


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
7) call mfield_coeffs() to convert coefficients to physical
   time units
*/

int
mfield_calc_nonlinear(gsl_vector *c, mfield_workspace *w)
{
  int s = 0;
  const mfield_parameters *params = &(w->params);
  const size_t max_iter = 50;     /* maximum iterations */
  const double xtol = 1.0e-4;
  const double gtol = 1.5e-3;
  const double ftol = 1.0e-6;
  int info;
  const size_t p = w->p;          /* number of coefficients */
  const size_t n = w->nres;       /* number of residuals */
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

  /* convert input vector from physical to dimensionless time units */
  mfield_coeffs(-1, c, c, w);

#if OLD_FDF
  /*
   * build and print residual histograms with previous coefficients
   * and previous wts_final vector
   */
  mfield_nonlinear_histogram(c, w);
#endif

  /* compute robust weights with coefficients from previous iteration */
  if (w->niter > 0)
    {
      /* compute residuals f = Y_model - y_data with previous coefficients */
      mfield_calc_f(c, w, w->fvec);

      /* compute robust weights */
      fprintf(stderr, "mfield_calc_nonlinear: computing robust weights...");
      mfield_robust_weights(w->fvec, w->wts_robust, w);
      fprintf(stderr, "done\n");

      /* compute final weights = wts_robust .* wts_spatial */
      gsl_vector_memcpy(w->wts_final, w->wts_robust);
      gsl_vector_mul(w->wts_final, w->wts_spatial);
    }
  else
    {
      gsl_vector_memcpy(w->wts_final, w->wts_spatial);
    }

  if (!params->use_weights)
    gsl_vector_set_all(w->wts_final, 1.0);

#if !OLD_FDF

  s = mfield_calc_nonlinear_multilarge(c, w);

#else

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
    mfield_coeffs(1, w->c, c, w);
  }
#endif

  /* convert coefficients to physical units */
  mfield_coeffs(1, w->c, c, w);

  printv_octave(c, "cfinal");

  w->niter++;

  return s;
} /* mfield_calc_nonlinear() */

/*
mfield_calc_nonlinear_multilarge()
  Calculate a solution to current inverse problem using multilarge

Inputs: c - coefficient vector
        w - workspace

Notes:
1) w->wts_final must be initialized prior to calling this function
2) On output, w->c contains the solution coefficients in dimensionless units
*/

static int
mfield_calc_nonlinear_multilarge(const gsl_vector *c, mfield_workspace *w)
{
  int s = 0;
  const mfield_parameters *params = &(w->params);
  const size_t max_iter = 50;     /* maximum iterations */
  const double xtol = 1.0e-6;
  const double gtol = 1.0e-6;
  const double ftol = 1.0e-6;
  int info;
  const size_t p = w->p;          /* number of coefficients */
  size_t n;                       /* number of residuals */
  gsl_multilarge_nlinear_fdf fdf;
  struct timeval tv0, tv1;
  double res0;                    /* initial residual */
  gsl_vector *f;

  n = w->nres;
  if (params->regularize == 1 && !params->synth_data)
    n += p;

  fdf.f = mfield_calc_Wf;
  fdf.df = mfield_calc_df2;
  fdf.fvv = mfield_calc_Wfvv;
  fdf.n = n;
  fdf.p = p;
  fdf.params = w;

  fprintf(stderr, "mfield_calc_nonlinear: precomputing vector J_int^T W J_int...");
  gettimeofday(&tv0, NULL);
  mfield_nonlinear_vector_precompute(w->wts_final, w);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (w->lls_solution == 1)
    {
      /*
       * There are no scalar residuals or Euler angles in the inverse problem,
       * so it is linear and we can use a LLS method
       */

      gsl_matrix *JTJ = w->JTJ_vec;
      gsl_vector *JTf = w->nlinear_workspace_p->g;
      size_t N = JTJ->size1;
      gsl_vector_view diag = gsl_matrix_diagonal(JTJ);
      gsl_permutation *perm = gsl_permutation_alloc(N);
      gsl_matrix *L = gsl_matrix_alloc(N, N); /* Cholesky factor */
      gsl_vector *work = gsl_vector_alloc(3 * N);
      double rcond;

      /* compute J^T f where f is the right hand side (residual vector when c = 0) */
      fprintf(stderr, "mfield_calc_nonlinear: computing RHS of linear system...");
      gsl_vector_set_zero(w->c);
      mfield_calc_Wf(w->c, w, w->wfvec);
      mfield_calc_df2(CblasTrans, w->c, w->wfvec, w, JTf, NULL);
      fprintf(stderr, "done\n");

      /* regularize JTJ matrix */
      gsl_vector_add(&diag.vector, w->LTL);

      fprintf(stderr, "mfield_calc_nonlinear: solving linear normal equations system...");
      gettimeofday(&tv0, NULL);

      lapack_cholesky_solve(JTJ, JTf, w->c, &rcond, L);

      gsl_vector_scale(w->c, -1.0);

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, cond(A) = %g)\n", time_diff(tv0, tv1), 1.0 / rcond);

      {
        const char *error_file = "error.txt";
        gsl_vector_const_view d = gsl_matrix_const_diagonal(L);
        FILE *fp;
        size_t n;

        /* compute (J^T J)^{-1} from Cholesky factor */
        fprintf(stderr, "mfield_calc_nonlinear: computing (J^T J)^{-1}...");
        gettimeofday(&tv0, NULL);

        lapack_cholesky_invert(L);

        gettimeofday(&tv1, NULL);
        fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

        fprintf(stderr, "mfield_calc_nonlinear: printing parameter uncertainties to %s...", error_file);

        fp = fopen(error_file, "w");

        n = 1;
        fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
        fprintf(fp, "# Field %zu: spherical harmonic order m\n", n++);
        fprintf(fp, "# Field %zu: uncertainty in g(n,m) (dimensionless)\n", n++);
        fprintf(fp, "# Field %zu: g(n,m) (nT)\n", n++);

        for (n = 1; n <= w->nmax_max; ++n)
          {
            int m, ni = (int) n;

            for (m = -ni; m <= ni; ++m)
              {
                size_t cidx = mfield_coeff_nmidx(n, m);
                double gnm = gsl_vector_get(w->c, cidx);
                double err_gnm = gsl_vector_get(&d.vector, cidx);

                fprintf(fp, "%5d %5zu %20.4e %20.4e\n", m, n, err_gnm, gnm);
              }

            fprintf(fp, "\n");
          }

        fclose(fp);

        fprintf(stderr, "done\n");
      }

      gsl_permutation_free(perm);
      gsl_vector_free(work);
      gsl_matrix_free(L);
    }
  else
    {
      fprintf(stderr, "mfield_calc_nonlinear: initializing multilarge...");
      gettimeofday(&tv0, NULL);
      gsl_multilarge_nlinear_init(c, &fdf, w->nlinear_workspace_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      /* compute initial residual */
      f = gsl_multilarge_nlinear_residual(w->nlinear_workspace_p);
      res0 = gsl_blas_dnrm2(f);

      fprintf(stderr, "mfield_calc_nonlinear: computing nonlinear least squares solution...");
      gettimeofday(&tv0, NULL);
      s = gsl_multilarge_nlinear_driver(max_iter, xtol, gtol, ftol,
                                        mfield_nonlinear_callback2, (void *) w,
                                        &info, w->nlinear_workspace_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      /* store final coefficients in dimensionless units in w->c */
      {
        gsl_vector *x_final = gsl_multilarge_nlinear_position(w->nlinear_workspace_p);
        gsl_vector_memcpy(w->c, x_final);
      }

      if (s == GSL_SUCCESS)
        {
          fprintf(stderr, "mfield_calc_nonlinear: NITER  = %zu\n",
                  gsl_multilarge_nlinear_niter(w->nlinear_workspace_p));
          fprintf(stderr, "mfield_calc_nonlinear: NFEV   = %zu\n", fdf.nevalf);
          fprintf(stderr, "mfield_calc_nonlinear: NJUEV  = %zu\n", fdf.nevaldfu);
          fprintf(stderr, "mfield_calc_nonlinear: NJTJEV = %zu\n", fdf.nevaldf2);
          fprintf(stderr, "mfield_calc_nonlinear: NAEV   = %zu\n", fdf.nevalfvv);
          fprintf(stderr, "mfield_calc_nonlinear: reason for stopping: %d\n", info);
          fprintf(stderr, "mfield_calc_nonlinear: initial |f(x)|: %.12e\n", res0);
          fprintf(stderr, "mfield_calc_nonlinear: final   |f(x)|: %.12e\n",
                  gsl_blas_dnrm2(f));
        }
      else
        {
          fprintf(stderr, "mfield_calc_nonlinear: failed: %s\n",
                  gsl_strerror(s));

          /* LM solver stalled or is converging slowly, switch to dogleg for next iteration */
          mfield_nonlinear_alloc_multilarge(gsl_multilarge_nlinear_trs_ddogleg, w);
        }
    } /* w->lls_solution == 0 */

  return s;
}

gsl_vector *
mfield_residual(const gsl_vector *c, mfield_workspace *w)
{
#if OLD_FDF
  gsl_vector *f = gsl_multifit_nlinear_residual(w->multifit_nlinear_p);
#else
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

#if OLD_FDF

  /* allocate fit workspace */
  {
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_parameters fdf_params =
      gsl_multifit_nlinear_default_parameters();

    fdf_params.solver = gsl_multifit_nlinear_solver_cholesky;
    fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
    fdf_params.trs = gsl_multifit_nlinear_trs_ddogleg;
    fdf_params.fdtype = GSL_MULTIFIT_NLINEAR_CTRDIFF;
    w->multifit_nlinear_p = gsl_multifit_nlinear_alloc(T, &fdf_params, w->nres_tot, p);
  }

#else

  /* allocate fit workspace - start with Levenberg-Marquardt solver */
  mfield_nonlinear_alloc_multilarge(gsl_multilarge_nlinear_trs_lm, w);

#endif

  w->wts_spatial = gsl_vector_alloc(nres);
  w->wts_robust = gsl_vector_alloc(nres);
  w->wts_final = gsl_vector_alloc(nres);

  gsl_vector_set_all(w->wts_robust, 1.0);

  /* to save memory, allocate p-by-p workspace since a full n-by-p isn't needed */
  w->robust_workspace_p = gsl_multifit_robust_alloc(gsl_multifit_robust_bisquare, p, p);

#if 0
  /* suggested by Nils to use 1.5 tuning constant */
  gsl_multifit_robust_tune(1.5, w->robust_workspace_p);
#endif

  w->fvec = gsl_vector_alloc(w->nres_tot);
  w->wfvec = gsl_vector_alloc(w->nres_tot);

  gsl_vector_set_all(w->wts_final, 1.0);

  /* calculate spatial weights */
  {
    size_t idx = 0;
    size_t j;

    fprintf(stderr, "mfield_init_nonlinear: calculating spatial weights...");

    for (i = 0; i < w->nsat; ++i)
      {
        magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

        for (j = 0; j < mptr->n; ++j)
          {
            double wt; /* spatial weight */

            if (MAGDATA_Discarded(mptr->flags[j]))
              continue;

            track_weight_get(mptr->phi[j], mptr->theta[j], &wt, w->weight_workspace_p);

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

      /* compute diag(L) */
      mfield_nonlinear_regularize(w->lambda_diag, w);

      /* compute L^T L */
      gsl_vector_memcpy(w->LTL, w->lambda_diag);
      gsl_vector_mul(w->LTL, w->lambda_diag);

      fprintf(stderr, "done\n");
    }
  else
    {
      gsl_vector_set_all(w->lambda_diag, 0.0);
      gsl_vector_set_all(w->LTL, 0.0);
    }

  return s;
} /* mfield_init_nonlinear() */

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
#if DEBUG
  fprintf(stderr, "mfield_calc_df2: entering function...\n");

  fprintf(stderr, "mfield_calc_df2: TransJ = %s...\n",
          TransJ == CblasTrans ? "trans" : "notrans");
#endif

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
          green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data,
                         w->green_array_p[thread_id]);

          /* calculate internal Green's functions for gradient point (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
            {
              green_calc_int(mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j],
                             vx_grad.vector.data, vy_grad.vector.data, vz_grad.vector.data,
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
  printv_octave(u, "u");
  printv_octave(v, "JTu");
  printsym_octave(JTJ, "JTJ");
  exit(1);
#endif

  gettimeofday(&tv1, NULL);
#if DEBUG
  fprintf(stderr, "mfield_calc_df2: leaving function... (%g seconds)\n", time_diff(tv0, tv1));
#endif

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
        extidx      - index of external field coefficient in [0,next-1]
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
                    const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, const size_t extidx,
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
          gsl_vector_view g_mf = gsl_vector_subvector(dB_int, 0, w->nnm_mf);

          v = gsl_vector_subvector(JTu, 0, w->nnm_mf);
          gsl_blas_daxpy(sWy, &g_mf.vector, &v.vector);
        }

      if (w->nnm_sv > 0)
        {
          gsl_vector_view g_sv = gsl_vector_subvector(dB_int, 0, w->nnm_sv);

          v = gsl_vector_subvector(JTu, w->sv_offset, w->nnm_sv);
          gsl_blas_daxpy(t * sWy, &g_sv.vector, &v.vector);
        }

      if (w->nnm_sa > 0)
        {
          gsl_vector_view g_sa = gsl_vector_subvector(dB_int, 0, w->nnm_sa);

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
                         const gsl_vector * u, const size_t ridx, gsl_vector * dB_int, gsl_vector * dB_int_grad,
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
          gsl_vector_view g_mf = gsl_vector_subvector(dB_int, 0, w->nnm_mf);
          gsl_vector_view dg_mf = gsl_vector_subvector(dB_int_grad, 0, w->nnm_mf);

          v = gsl_vector_subvector(JTu, 0, w->nnm_mf);
          gsl_blas_daxpy(sWy, &dg_mf.vector, &v.vector);
          gsl_blas_daxpy(-sWy, &g_mf.vector, &v.vector);
        }

      if (w->nnm_sv > 0)
        {
          gsl_vector_view g_sv = gsl_vector_subvector(dB_int, 0, w->nnm_sv);
          gsl_vector_view dg_sv = gsl_vector_subvector(dB_int_grad, 0, w->nnm_sv);

          v = gsl_vector_subvector(JTu, w->sv_offset, w->nnm_sv);
          gsl_blas_daxpy(t_grad * sWy, &dg_sv.vector, &v.vector);
          gsl_blas_daxpy(-t * sWy, &g_sv.vector, &v.vector);
        }

      if (w->nnm_sa > 0)
        {
          gsl_vector_view g_sa = gsl_vector_subvector(dB_int, 0, w->nnm_sa);
          gsl_vector_view dg_sa = gsl_vector_subvector(dB_int_grad, 0, w->nnm_sa);

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
        extidx      - index of external field coefficient in [0,next-1]
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
        extidx      - index of external field coefficient in [0,next-1]
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
mfield_nonlinear_vector_precompute()
  Precompute J_int^T W J_int for vector measurements, since
this submatrix is independent of the model parameters and
only needs to be computed once per iteration. This function
uses OpenMP to speed up the calculation

Inputs: weights - weight vector
        w       - workspace

Notes:
1) w->JTJ_vec is updated with J_int^T W J_int for vector
residuals
*/

static int
mfield_nonlinear_vector_precompute(const gsl_vector *weights, mfield_workspace *w)
{
  int s = GSL_SUCCESS;
  size_t i, j;
  size_t *omp_nrows; /* number of rows processed by each thread */
  size_t nres_vec = w->nres_vec + w->nres_vec_grad;

  gsl_matrix_set_zero(w->JTJ_vec);

  /* check for quick return */
  if (nres_vec == 0)
    return GSL_SUCCESS;

  omp_nrows = malloc(w->max_threads * sizeof(size_t));

  /*
   * omp_rowidx[thread_id] contains the number of currently filled rows
   * of omp_J[thread_id]. When omp_J[thread_id] is filled, it is folded
   * into the JTJ_vec matrix and then omp_rowidx[thread_id] is reset to 0
   */
  for (i = 0; i < w->max_threads; ++i)
    {
      omp_nrows[i] = 0;
      w->omp_rowidx[i] = 0;
    }

  fprintf(stderr, "\n");

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */
          double t = mptr->ts[j];
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];

          gsl_vector_view vx = gsl_matrix_row(w->omp_dX, thread_id);
          gsl_vector_view vy = gsl_matrix_row(w->omp_dY, thread_id);
          gsl_vector_view vz = gsl_matrix_row(w->omp_dZ, thread_id);

          gsl_vector_view vx_grad = gsl_matrix_row(w->omp_dX_grad, thread_id);
          gsl_vector_view vy_grad = gsl_matrix_row(w->omp_dY_grad, thread_id);
          gsl_vector_view vz_grad = gsl_matrix_row(w->omp_dZ_grad, thread_id);

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* calculate internal Green's functions */
          green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data,
                         w->green_array_p[thread_id]);

          /* calculate internal Green's functions for gradient point (N/S or E/W) */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW))
            {
              green_calc_int(mptr->r_ns[j], mptr->theta_ns[j], mptr->phi_ns[j],
                             vx_grad.vector.data, vy_grad.vector.data, vz_grad.vector.data,
                             w->green_array_p[thread_id]);
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green(t, wj, &vx.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green(t, wj, &vy.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green(t, wj, &vz.vector, &v.vector, w);
                }
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green_grad(t, mptr->ts_ns[j], wj, &vx.vector, &vx_grad.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green_grad(t, mptr->ts_ns[j], wj, &vy.vector, &vy_grad.vector, &v.vector, w);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              double wj = gsl_vector_get(weights, ridx++);
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view v = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
                  mfield_vector_green_grad(t, mptr->ts_ns[j], wj, &vz.vector, &vz_grad.vector, &v.vector, w);
                }
            }

          /*
           * check if omp_J[thread_id] is full and should be folded into JTJ; the
           * 15 is just some slop to prevent trying to fill rows past the matrix buffer
           * in the loop above
           */
          if (w->omp_rowidx[thread_id] >= w->omp_J[thread_id]->size1 - 15)
            {
              /* fold current matrix block into JTJ_vec, one thread at a time */
              gsl_matrix_view m = gsl_matrix_submatrix(w->omp_J[thread_id], 0, 0, w->omp_rowidx[thread_id], w->p_int);

              /* keep cumulative total of rows processed by this thread for progress bar */
              omp_nrows[thread_id] += w->omp_rowidx[thread_id];
              w->omp_rowidx[thread_id] = 0;

#pragma omp critical
              {
                gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &m.matrix, 1.0, w->JTJ_vec);
              }

              if (thread_id == 0)
                {
                  double progress = 0.0;
                  size_t k;

                  for (k = 0; k < w->max_threads; ++k)
                    progress += (double) omp_nrows[k];

                  progress /= (double) nres_vec;

                  fprintf(stderr, "\t");
                  progress_bar(stderr, progress, 70);
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
        }
    }

  free(omp_nrows);

#if 0
  printsym_octave(w->JTJ_vec, "JTJ_vec");
  exit(1);
#endif

  return s;
}

/*
mfield_vector_green()
  Function to compute sqrt(w) [ J_mf J_sv J_sa ] for a given set of
vector Green's functions

Inputs: t      - scaled timestamp
        weight - weight for this measurement
        g      - vector Green's functions J_mf, size nnm_max
        G      - (output) combined vector G = sqrt(w) [ g ; t*g ; 0.5*t*t*g ],
                 size w->p_int
        w      - workspace
*/

static int
mfield_vector_green(const double t, const double weight, const gsl_vector *g,
                    gsl_vector *G, mfield_workspace *w)
{
  const double sqrt_weight = sqrt(weight);
  size_t i;

  /* form G */
  for (i = 0; i < w->nnm_max; ++i)
    {
      double gi = sqrt_weight * gsl_vector_get(g, i);

      mfield_set_mf(G, i, gi, w);
      mfield_set_sv(G, i, t * gi, w);
      mfield_set_sa(G, i, 0.5 * t * t * gi, w);
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
mfield_nonlinear_histogram()
  Print residual histogram

Inputs: c - scaled/dimensionless coefficient vector
        w - workspace

Notes:
1) w->wts_final must be initialized prior to calling this function
*/

static int
mfield_nonlinear_histogram(const gsl_vector *c, mfield_workspace *w)
{
  int s = 0;
  FILE *fp;
  char filename[2048];

  /* reset histograms */
  gsl_histogram_reset(w->hf);
  gsl_histogram_reset(w->hz);

  /* loop through data and construct residual histograms */
  mfield_calc_f(c, w, NULL);

  /* scale histograms */
  gsl_histogram_scale(w->hf, 1.0 / gsl_histogram_sum(w->hf));
  gsl_histogram_scale(w->hz, 1.0 / gsl_histogram_sum(w->hz));

  /* print histograms to file */

  sprintf(filename, "reshistF.nlin.iter%zu.dat", w->niter);
  fprintf(stderr, "mfield_nonlinear_histogram: writing %s...", filename);
  fp = fopen(filename, "w");
  mfield_print_histogram(fp, w->hf);
  fclose(fp);
  fprintf(stderr, "done\n");

  sprintf(filename, "reshistZ.nlin.iter%zu.dat", w->niter);
  fprintf(stderr, "mfield_nonlinear_histogram: writing %s...", filename);
  fp = fopen(filename, "w");
  mfield_print_histogram(fp, w->hz);
  fclose(fp);
  fprintf(stderr, "done\n");

  return s;
} /* mfield_nonlinear_histogram() */

/*
mfield_nonlinear_regularize()
  Construct diag = diag(L) for regularized fit by minimizing
[d/dt B_r]^2 and/or [d^2/dt^2 B_r]^2 at the CMB

< | d^2/dt^2 B_r |^2 > = \int | d^2/dt^2 B_r |^2 dOmega = || L g ||^2

with

L_{nm} = (n+1)/sqrt(2n + 1) (a/c)^{n+2}
*/

static int
mfield_nonlinear_regularize(gsl_vector *diag, mfield_workspace *w)
{
  int s = 0;
  const size_t nmax = GSL_MAX(w->nmax_sv, w->nmax_sa);
  const double c = 3485.0;       /* Earth core radius */
  const double a = w->params.R;  /* Earth surface radius */
  const double ratio = a / c;
  size_t n;

  gsl_vector_set_zero(diag);

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;
      int m;
      double term = (n + 1.0) / sqrt(2.0*n + 1.0) * pow(ratio, n + 2.0);

      for (m = -ni; m <= ni; ++m)
        {
          size_t cidx = mfield_coeff_nmidx(n, m);
          mfield_set_sv(diag, cidx, w->lambda_sv * term, w);
          mfield_set_sa(diag, cidx, w->lambda_sa * term, w);
        }
    }

  return s;
} /* mfield_nonlinear_regularize() */

static int
mfield_robust_print_stat(const char *str, const double sigma, const gsl_rstat_workspace *rstat_p)
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
   * (X,Y,Z,F,DX,DY,DZ)
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
              double wi = bisquare(fi / (tune * sigma_X));
              gsl_vector_set(wts, idx++, wi);
              gsl_rstat_add(wi, rstat_x[i]);
            }

          if (MAGDATA_ExistY(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = bisquare(fi / (tune * sigma_Y));
              gsl_vector_set(wts, idx++, wi);
              gsl_rstat_add(wi, rstat_y[i]);
            }

          if (MAGDATA_ExistZ(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = bisquare(fi / (tune * sigma_Z));
              gsl_vector_set(wts, idx++, wi);
              gsl_rstat_add(wi, rstat_z[i]);
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              double fi = gsl_vector_get(f, idx);
              double wi = bisquare(fi / (tune * sigma_F));
              gsl_vector_set(wts, idx++, wi);
              gsl_rstat_add(wi, rstat_f[i]);
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

  fprintf(stderr, "iteration %zu:\n", iter);

  if (w->nnm_mf > 0)
    {
      fprintf(stderr, "\t %-10s %12.4f %12.4f %12.4f [nT]\n",
              "dipole:",
              mfield_get_mf(x, mfield_coeff_nmidx(1, 0), w),
              mfield_get_mf(x, mfield_coeff_nmidx(1, 1), w),
              mfield_get_mf(x, mfield_coeff_nmidx(1, -1), w));
    }

  if (w->nnm_sv > 0)
    {
      fprintf(stderr, "\t %-10s %12.4f %12.4f %12.4f [nT]\n",
              "SV dipole:",
              mfield_get_sv(x, mfield_coeff_nmidx(1, 0), w),
              mfield_get_sv(x, mfield_coeff_nmidx(1, 1), w),
              mfield_get_sv(x, mfield_coeff_nmidx(1, -1), w));
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
              size_t euler_idx = mfield_euler_idx(i, t0, w);

              fprintf(stderr, "\t euler : %12.4f %12.4f %12.4f [deg]\n",
                      gsl_vector_get(x, euler_idx) * 180.0 / M_PI,
                      gsl_vector_get(x, euler_idx + 1) * 180.0 / M_PI,
                      gsl_vector_get(x, euler_idx + 2) * 180.0 / M_PI);
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

  fprintf(stderr, "\t dipole: %12.4f %12.4f %12.4f [nT]\n",
          gsl_vector_get(x, mfield_coeff_nmidx(1, 0)),
          gsl_vector_get(x, mfield_coeff_nmidx(1, 1)),
          gsl_vector_get(x, mfield_coeff_nmidx(1, -1)));

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
              size_t euler_idx = mfield_euler_idx(i, t0, w);
              char buf[32];

              if (mptr->euler_flags & EULER_FLG_ZYX)
                sprintf(buf, "ZYX");
              else if (mptr->euler_flags & EULER_FLG_ZYZ)
                sprintf(buf, "ZYZ");
              else
                *buf = '\0';

              fprintf(stderr, "\t euler%zu: %12.4f %12.4f %12.4f [deg] [%s]\n",
                      i,
                      gsl_vector_get(x, euler_idx) * 180.0 / M_PI,
                      gsl_vector_get(x, euler_idx + 1) * 180.0 / M_PI,
                      gsl_vector_get(x, euler_idx + 2) * 180.0 / M_PI,
                      buf);
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

/* allocate w->nlinear_workspace_p with given TRS method; this function makes
 * it easy to switch TRS methods in the middle of an iteration */
static int
mfield_nonlinear_alloc_multilarge(const gsl_multilarge_nlinear_trs * trs,
                                  mfield_workspace * w)
{
  const gsl_multilarge_nlinear_type *T = gsl_multilarge_nlinear_trust;
  gsl_multilarge_nlinear_parameters fdf_params =
    gsl_multilarge_nlinear_default_parameters();

  if (w->nlinear_workspace_p)
    gsl_multilarge_nlinear_free(w->nlinear_workspace_p);

  fdf_params.trs = trs;
  fdf_params.scale = gsl_multilarge_nlinear_scale_levenberg;

  w->nlinear_workspace_p = gsl_multilarge_nlinear_alloc(T, &fdf_params, w->nres_tot, w->p);

  return 0;
}
