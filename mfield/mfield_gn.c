/*
 * mfield_gn.c
 *
 * Gauss-Newton solver for nonlinear least squares
 */

static int mfield_nonlinear_alloc_gn(mfield_workspace * w);
static int mfield_calc_nonlinear_gn(const gsl_vector *c, mfield_workspace *w);

/* allocate w->nlinear_workspace_p for Gauss-Newton solver */
static int
mfield_nonlinear_alloc_gn(mfield_workspace * w)
{
  const gsl_multilarge_nlinear_type *T = gsl_multilarge_nlinear_trust;
  gsl_multilarge_nlinear_parameters fdf_params =
    gsl_multilarge_nlinear_default_parameters();

  if (w->nlinear_workspace_p)
    gsl_multilarge_nlinear_free(w->nlinear_workspace_p);

  fdf_params.trs = gsl_multilarge_nlinear_trs_lm;
  fdf_params.scale = gsl_multilarge_nlinear_scale_levenberg;

  w->nlinear_workspace_p = gsl_multilarge_nlinear_alloc(T, &fdf_params, w->nres_tot, w->p);

  return 0;
}

/*
mfield_calc_nonlinear_gn()
  Calculate a solution to current inverse problem using Gauss-Newton

Inputs: c - coefficient vector
        w - workspace

Notes:
1) w->wts_final must be initialized prior to calling this function
2) On output, w->c contains the solution coefficients in dimensionless units
*/

static int
mfield_calc_nonlinear_gn(const gsl_vector *c, mfield_workspace *w)
{
  int s = 0;
  const mfield_parameters *params = &(w->params);
  int info;
  const size_t p = w->p;          /* number of coefficients */
  size_t n;                       /* number of residuals */
  gsl_multilarge_nlinear_fdf fdf;
  struct timeval tv0, tv1;

  n = w->nres;
  if (params->regularize == 1 && !params->synth_data)
    n += p;

  fdf.f = mfield_calc_Wf;
  fdf.df = mfield_calc_df3;
  fdf.fvv = mfield_calc_Wfvv;
  fdf.n = n;
  fdf.p = p;
  fdf.params = w;

  fprintf(stderr, "mfield_calc_nonlinear_gn: precomputing vector J_int^T W J_int...");
  gettimeofday(&tv0, NULL);
  mfield_nonlinear_precompute(w->sqrt_wts_final, w);
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
      double rcond;

      /* compute J^T f where f is the right hand side (residual vector when c = 0) */
      fprintf(stderr, "mfield_calc_nonlinear_gn: computing RHS of linear system...");
      gsl_vector_set_zero(w->c);
      mfield_calc_Wf(w->c, w, w->wfvec);
      mfield_calc_df3(CblasTrans, w->c, w->wfvec, w, JTf, NULL);
      fprintf(stderr, "done\n");

      /* regularize JTJ matrix */
      gsl_spmatrix_add_to_dense(JTJ, w->Lambda);

      fprintf(stderr, "mfield_calc_nonlinear_gn: solving linear normal equations system...");
      gettimeofday(&tv0, NULL);

      info = lapack_cholesky_solve(JTJ, JTf, w->c, &rcond, w->choleskyL);

      gsl_vector_scale(w->c, -1.0);

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, cond(A) = %g, status = %d)\n", time_diff(tv0, tv1), 1.0 / rcond, info);
    }
  else
    {
      const double xtol = 1.0e-6;
      const double gtol = 1.0e-6;
      const double ftol = 1.0e-6;
      static size_t iter = 0;
      const double alpha = 1.0;
      gsl_matrix *JTJ = w->nlinear_workspace_p->JTJ;
      gsl_vector *JTf = w->nlinear_workspace_p->g;
      gsl_vector *dx = w->nlinear_workspace_p->dx;
      gsl_vector *x = w->nlinear_workspace_p->x;
      double rcond;

      /* compute f, JTf, and JTJ */
      fprintf(stderr, "mfield_calc_nonlinear_gn: initializing GN...");
      gettimeofday(&tv0, NULL);
      gsl_multilarge_nlinear_init(c, &fdf, w->nlinear_workspace_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      fprintf(stderr, "mfield_calc_nonlinear_gn: solving normal equations system...");
      gettimeofday(&tv0, NULL);

      info = lapack_cholesky_solve(JTJ, JTf, dx, &rcond, NULL);

      /* c_{i+1} = c_i - alpha*dx */
      gsl_vector_memcpy(w->c, dx);
      gsl_vector_scale(w->c, -alpha);
      gsl_vector_add(w->c, c);

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, cond(A) = %g, status = %d)\n", time_diff(tv0, tv1), 1.0 / rcond, info);

      mfield_nonlinear_callback2(++iter, w, w->nlinear_workspace_p);

      gsl_vector_memcpy(x, w->c);
      s = gsl_multilarge_nlinear_test(xtol, gtol, ftol, &info, w->nlinear_workspace_p);

      if (s == GSL_SUCCESS)
        {
          fprintf(stderr, "CONVERGED\n");
        }
      else
        fprintf(stderr, "NOT CONVERGED\n");
    } /* w->lls_solution == 0 */

  return s;
}
