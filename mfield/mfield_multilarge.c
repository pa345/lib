/*
 * mfield_multilarge.c
 *
 * Contains routines for fitting magnetic field module using gsl_multilarge_nlinear framework
 */

#include <bspline2/gsl_bspline2.h>

static int mfield_calc_nonlinear_multilarge(const gsl_vector *c, mfield_workspace *w);
static int mfield_nonlinear_alloc_multilarge(const gsl_multilarge_nlinear_trs * trs, mfield_workspace * w);
static int mfield_nonlinear_precompute_core(const gsl_vector *sqrt_weights, gsl_matrix * JTJ_core, mfield_workspace *w);
static int mfield_nonlinear_precompute_crust(const gsl_vector *sqrt_weights, gsl_matrix * JTJ_crust, mfield_workspace *w);
static int mfield_nonlinear_precompute_mixed(const gsl_vector *sqrt_weights, gsl_matrix * JTJ_mixed, mfield_workspace *w);
static int mfield_nonlinear_precompute(const gsl_vector *sqrt_weights, mfield_workspace *w);
static int mfield_calc_df3(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                           void *params, gsl_vector * v, gsl_matrix * JTJ);
static int mfield_nonlinear_scalar_core(const gsl_vector * x, const gsl_vector * sqrt_weights, gsl_matrix * JTJ_core, mfield_workspace *w);
static int mfield_nonlinear_scalar_crust(const gsl_vector * x, const gsl_vector * sqrt_weights, gsl_matrix * JTJ_core, mfield_workspace *w);
static inline int mfield_jacobian_J1Tu(const gsl_vector * N, const size_t istart, const double sqrt_wj, const double uj,
                                       const gsl_vector * dB_int, gsl_vector *JTu, const mfield_workspace *w);
static inline int mfield_jacobian_J1Tu_scalar(const gsl_vector * N, const size_t istart, const double uj,
                                              const gsl_vector * dB, gsl_vector *JTu, const mfield_workspace *w);
static inline int mfield_jacobian_J1u(const double t, const double sqrt_wt, const gsl_vector * u,
                                      const size_t ridx, const gsl_vector * dB_int, gsl_vector *Ju,
                                      const mfield_workspace *w);
static inline int mfield_jacobian_F(CBLAS_TRANSPOSE_t TransJ, const size_t istart, const gsl_vector * u,
                                    const size_t ridx, const gsl_spmatrix * J2, const gsl_vector *J_int,
                                    gsl_vector *v, gsl_matrix * J2TJ1, const mfield_workspace *w);
static int mfield_jacobian_J2TJ1(const gsl_vector * N, const size_t istart, const double sqrt_wj, const size_t ridx,
                                 const gsl_vector * dB_int, const gsl_spmatrix * J2, gsl_matrix * J2TJ1,
                                 const mfield_workspace * w);

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
  const double xtol = 1.0e-6;
  const double gtol = 1.0e-6;
  const double ftol = 1.0e-6;
  int info;
  const size_t p = w->p;          /* number of coefficients */
  size_t n;                       /* number of residuals */
  size_t max_iter = 5;            /* maximum "inner" iterations */
  gsl_multilarge_nlinear_fdf fdf;
  struct timeval tv0, tv1;
  double res0;                    /* initial residual */
  gsl_vector *f;

  /*
   * On the first iteration (niter = 0), there are no robust weights yet, use more inner iterations.
   * On the second iteration (niter = 1), it is the first set of robust weights, so also use more inner iterations.
   * For the third iteration and onward, they are usually minor refinements, so we can use less inner iterations.
   */
  if (w->niter <= 1)
    max_iter = 30;
  else
    max_iter = 5;

  n = w->nres;
  if (params->regularize == 1 && !params->synth_data)
    n += p;

  fdf.f = mfield_calc_Wf;
  fdf.df = mfield_calc_df3;
  fdf.fvv = mfield_calc_Wfvv;
  fdf.n = n;
  fdf.p = p;
  fdf.params = w;

  fprintf(stderr, "mfield_calc_nonlinear: precomputing vector J_int^T W J_int...");
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
      fprintf(stderr, "mfield_calc_nonlinear: computing RHS of linear system...");
      gsl_vector_set_zero(w->c);
      mfield_calc_Wf(w->c, w, w->wfvec);
      mfield_calc_df3(CblasTrans, w->c, w->wfvec, w, JTf, NULL);
      fprintf(stderr, "done\n");

      /* regularize JTJ matrix */
      gsl_spmatrix_add_to_dense(JTJ, w->Lambda);

      fprintf(stderr, "mfield_calc_nonlinear: solving linear normal equations system...");
      gettimeofday(&tv0, NULL);

      lapack_cholesky_solve(JTJ, JTf, w->c, &rcond, w->choleskyL);

      gsl_vector_scale(w->c, -1.0);

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, cond(A) = %g)\n", time_diff(tv0, tv1), 1.0 / rcond);

#if 0
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
#endif
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
        }
    } /* w->lls_solution == 0 */

  return s;
}

/*
mfield_nonlinear_precompute_core()
  Compute J_core^T W J_core for vector residuals

Inputs: sqrt_weights - sqrt(wts), length nres
        JTJ_core     - (output) output matrix, p_core-by-p_core
        w            - workspace

Notes:
1) Lower triangle of JTJ_core is filled
2) JTJ_core must be initialized to zero by calling function
*/

static int
mfield_nonlinear_precompute_core(const gsl_vector *sqrt_weights, gsl_matrix * JTJ_core,
                                 mfield_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax_core;
  const size_t mmax = nmax;
  const size_t ncontrol = gsl_bspline2_ncontrol(w->gauss_spline_workspace_p[0]);
  const size_t order = gsl_bspline2_order(w->gauss_spline_workspace_p[0]);
  size_t ii, jj, i, j;
  struct timeval tv0, tv1;

  /* compute J_core^T W J_core
   *
   * This matrix can be divided into nnm_core-by-nnm_core blocks:
   *
   * A_{ij} = \sum_{k=1}^{nres} N_i(t_k) N_j(t_k) B(r_k) B^T(r_k)
   *
   * We have outer loops over i,j, followed by threaded loop over the residuals,
   * accumulate the internal Greens functions into a large block matrix T,
   * which is then folded into A_{ij} with DSYRK when it is full.
   */

  fprintf(stderr, "\n");
  fprintf(stderr, "mfield_nonlinear_precompute_core: computing J_core^T W J_core...");
  fprintf(stderr, "\n");
  gettimeofday(&tv0, NULL);

  for (ii = 0; ii < ncontrol; ++ii)
    {
      for (jj = 0; jj <= ii; ++jj)
        {
          gsl_matrix_view JTJ_block = gsl_matrix_submatrix(JTJ_core, ii * w->nnm_core, jj * w->nnm_core,
                                                           w->nnm_core, w->nnm_core);

          for (i = 0; i < w->max_threads; ++i)
            w->omp_rowidx[i] = 0;

          for (i = 0; i < w->nsat; ++i)
            {
              magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
              for (j = 0; j < mptr->n; ++j)
                {
                  int thread_id, flag;
                  gsl_vector *N_gauss;
                  size_t ridx, left, istart;
                  double t, Nii, Njj, alpha;
                  gsl_vector *VX, *VY, *VZ;
                  gsl_vector_view vx, vy, vz;

                  if (MAGDATA_Discarded(mptr->flags[j]))
                    continue;

                  if (!MAGDATA_FitMF(mptr->flags[j]))
                    continue;

                  if (!(mptr->flags[j] & (MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z)))
                    continue;

                  thread_id = omp_get_thread_num();
                  ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */
                  N_gauss = w->gauss_spline_workspace_p[thread_id]->B;
                  VX = VY = VZ = NULL;

                  t = epoch2year(mptr->t[j]);
                  left = gsl_bspline2_find_interval(t, &flag, w->gauss_spline_workspace_p[thread_id]);
                  if (flag != 0)
                    continue;

                  if (ii < (left + 1) - order || ii > left)
                    continue; /* N_{ii}(t) = 0 */
                  if (jj < (left + 1) - order || jj > left)
                    continue; /* N_{ii}(t) = 0 */

                  if (mptr->flags[j] & MAGDATA_FLG_X)
                    {
                      vx = gsl_matrix_subrow(w->omp_T[thread_id], w->omp_rowidx[thread_id]++, 0, w->nnm_core);
                      VX = &vx.vector;
                    }

                  if (mptr->flags[j] & MAGDATA_FLG_Y)
                    {
                      vy = gsl_matrix_subrow(w->omp_T[thread_id], w->omp_rowidx[thread_id]++, 0, w->nnm_core);
                      VY = &vy.vector;
                    }

                  if (mptr->flags[j] & MAGDATA_FLG_Z)
                    {
                      vz = gsl_matrix_subrow(w->omp_T[thread_id], w->omp_rowidx[thread_id]++, 0, w->nnm_core);
                      VZ = &vz.vector;
                    }

                  /* calculate internal Green's functions */
                  green_calc_intopt(1, nmax, mmax, mptr->r[j], mptr->theta[j], mptr->phi[j], VX, VY, VZ, w->green_array_p[thread_id]);

                  /* calculate non-zero B-splines for the Gauss coefficients for this timestamp */
                  gsl_bspline2_eval_basis_nonzero(t, N_gauss, &istart, w->gauss_spline_workspace_p[thread_id]);

                  Nii = gsl_vector_get(N_gauss, ii - istart);
                  Njj = gsl_vector_get(N_gauss, jj - istart);
                  alpha = sqrt(Nii * Njj);

                  if (mptr->flags[j] & MAGDATA_FLG_X)
                    {
                      double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
                      gsl_blas_dscal(alpha * sqrt_wj, VX);
                    }

                  if (mptr->flags[j] & MAGDATA_FLG_Y)
                    {
                      double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
                      gsl_blas_dscal(alpha * sqrt_wj, VY);
                    }

                  if (mptr->flags[j] & MAGDATA_FLG_Z)
                    {
                      double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
                      gsl_blas_dscal(alpha * sqrt_wj, VZ);
                    }

                  /*
                   * check if omp_T[thread_id] is full and should be folded into JTJ; the
                   * 15 is just some slop to prevent trying to fill columns past the matrix buffer
                   * in the loop above
                   */
                  if (w->omp_rowidx[thread_id] >= w->omp_T[thread_id]->size1 - 15)
                    {
                      /* fold current matrix block into JTJ_vec, one thread at a time */
                      gsl_matrix_view m = gsl_matrix_submatrix(w->omp_T[thread_id], 0, 0, w->omp_rowidx[thread_id], w->nnm_core);

#pragma omp critical
                      {
                        gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &m.matrix, 1.0, &JTJ_block.matrix);
                      }

                      w->omp_rowidx[thread_id] = 0;
                    }
                } /* for (j = 0; j < mptr->n; ++j) */
            } /* for (i = 0; i < w->nsat; ++i) */

          for (i = 0; i < w->max_threads; ++i)
            {
              if (w->omp_rowidx[i] > 0)
                {
                  /* accumulate final Green's functions into this (ii,jj) block */
                  gsl_matrix_view m = gsl_matrix_submatrix(w->omp_T[i], 0, 0, w->omp_rowidx[i], w->nnm_core);
                  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &m.matrix, 1.0, &JTJ_block.matrix);
                }
            }

          /* at this point the (ii, jj) block of JTJ is calculated, but if it is not
           * a diagonal block, we need to copy the lower triangle to the upper */
          if (ii != jj)
            {
              gsl_matrix_transpose_tricpy('L', 0, &JTJ_block.matrix, &JTJ_block.matrix);
            }

        } /* for (jj = 0; jj <= ii; ++jj) */

      progress_bar(stderr, (double) ii / (double) ncontrol, 70);
    } /* for (ii = 0; ii < ncontrol; ++ii) */

  progress_bar(stderr, 1.0, 70);
  fprintf(stderr, "\n");

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  return s;
}

/*
mfield_nonlinear_precompute_crust()
  Compute J_crust^T W J_crust for vector residuals

Inputs: sqrt_weights - sqrt(wts), length nres
        JTJ_crust    - (output) output matrix, p_crust-by-p_crust
        w            - workspace

Notes:
1) Lower triangle of JTJ_crust is filled
2) JTJ_crust must be initialized to zero by calling function
*/

static int
mfield_nonlinear_precompute_crust(const gsl_vector *sqrt_weights, gsl_matrix * JTJ_crust,
                                  mfield_workspace *w)
{
  int s = 0;
  const size_t nmin = w->nmax_core + 1;
  const size_t nmax = w->nmax;
  const size_t mmax = nmax;
  size_t i, j;
  struct timeval tv0, tv1;

  if (nmin > nmax)
    return s; /* nothing to do */

  /* compute J_crust^T W J_crust
   *
   * Threaded loop over the residuals, accumulate the internal Greens
   * functions into a large block matrix T, which is then folded into
   * JTJ_crust with DSYRK when it is full.
   */

  fprintf(stderr, "mfield_nonlinear_precompute_crust: computing J_crust^T W J_crust...");
  fprintf(stderr, "\n");
  gettimeofday(&tv0, NULL);

  for (i = 0; i < w->max_threads; ++i)
    w->omp_rowidx[i] = 0;

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id;
          size_t ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */
          gsl_vector *VX, *VY, *VZ;
          gsl_vector_view vx, vy, vz;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (!MAGDATA_FitMF(mptr->flags[j]))
            continue;

          if (!(mptr->flags[j] & (MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z)))
            continue;

          thread_id = omp_get_thread_num();
          VX = VY = VZ = NULL;

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              vx = gsl_matrix_subrow(w->omp_T[thread_id], w->omp_rowidx[thread_id]++, 0, w->nnm_crust);
              VX = &vx.vector;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              vy = gsl_matrix_subrow(w->omp_T[thread_id], w->omp_rowidx[thread_id]++, 0, w->nnm_crust);
              VY = &vy.vector;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              vz = gsl_matrix_subrow(w->omp_T[thread_id], w->omp_rowidx[thread_id]++, 0, w->nnm_crust);
              VZ = &vz.vector;
            }

          /* calculate internal Green's functions for crustal field only */
          green_calc_intopt(nmin, nmax, mmax, mptr->r[j], mptr->theta[j], mptr->phi[j], VX, VY, VZ, w->green_array_p[thread_id]);

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              gsl_blas_dscal(sqrt_wj, VX);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              gsl_blas_dscal(sqrt_wj, VY);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              gsl_blas_dscal(sqrt_wj, VZ);
            }

          /*
           * check if omp_T[thread_id] is full and should be folded into JTJ; the
           * 15 is just some slop to prevent trying to fill columns past the matrix buffer
           * in the loop above
           */
          if (w->omp_rowidx[thread_id] >= w->omp_T[thread_id]->size1 - 15)
            {
              /* fold current matrix block into JTJ_vec, one thread at a time */
              gsl_matrix_view m = gsl_matrix_submatrix(w->omp_T[thread_id], 0, 0, w->omp_rowidx[thread_id], w->nnm_crust);

#pragma omp critical
              {
                gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &m.matrix, 1.0, JTJ_crust);
              }

              w->omp_rowidx[thread_id] = 0;
            }
        } /* for (j = 0; j < mptr->n; ++j) */

      progress_bar(stderr, (i + 1.0) / (double) w->nsat, 70);
    } /* for (i = 0; i < w->nsat; ++i) */

  for (i = 0; i < w->max_threads; ++i)
    {
      if (w->omp_rowidx[i] > 0)
        {
          /* accumulate final Green's functions */
          gsl_matrix_view m = gsl_matrix_submatrix(w->omp_T[i], 0, 0, w->omp_rowidx[i], w->nnm_crust);
          gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &m.matrix, 1.0, JTJ_crust);
        }
    }

  progress_bar(stderr, 1.0, 70);
  fprintf(stderr, "\n");

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  return s;
}

static int
mfield_nonlinear_precompute_crust2(const gsl_vector *sqrt_weights, gsl_matrix * JTJ_crust,
                                  mfield_workspace *w)
{
  int s = 0;
  const size_t nmin = w->nmax_core + 1;
  const size_t nmax = w->nmax;
  const size_t mmax = nmax;
  size_t i, j;
  struct timeval tv0, tv1;

  if (nmin > nmax)
    return s; /* nothing to do */

  /* compute J_crust^T W J_crust
   *
   * Threaded loop over the residuals, accumulate the internal Greens
   * functions into a large block matrix T, which is then folded into
   * JTJ_crust with DSYRK when it is full.
   */

  fprintf(stderr, "mfield_nonlinear_precompute_crust2: computing J_crust^T W J_crust...");
  fprintf(stderr, "\n");
  gettimeofday(&tv0, NULL);

  for (i = 0; i < w->max_threads; ++i)
    w->omp_colidx[i] = 0;

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id;
          size_t ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */
          gsl_vector *VX, *VY, *VZ;
          gsl_vector_view vx, vy, vz;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (!MAGDATA_FitMF(mptr->flags[j]))
            continue;

          if (!(mptr->flags[j] & (MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z)))
            continue;

          thread_id = omp_get_thread_num();
          VX = VY = VZ = NULL;

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              vx = gsl_matrix_subcolumn(w->omp_S[thread_id], w->omp_colidx[thread_id]++, 0, w->nnm_crust);
              VX = &vx.vector;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              vy = gsl_matrix_subcolumn(w->omp_S[thread_id], w->omp_colidx[thread_id]++, 0, w->nnm_crust);
              VY = &vy.vector;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              vz = gsl_matrix_subcolumn(w->omp_S[thread_id], w->omp_colidx[thread_id]++, 0, w->nnm_crust);
              VZ = &vz.vector;
            }

          /* calculate internal Green's functions for crustal field only */
          green_calc_intopt(nmin, nmax, mmax, mptr->r[j], mptr->theta[j], mptr->phi[j], VX, VY, VZ, w->green_array_p[thread_id]);

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              gsl_blas_dscal(sqrt_wj, VX);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              gsl_blas_dscal(sqrt_wj, VY);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
              gsl_blas_dscal(sqrt_wj, VZ);
            }

          /*
           * check if omp_S[thread_id] is full and should be folded into JTJ; the
           * 15 is just some slop to prevent trying to fill columns past the matrix buffer
           * in the loop above
           */
          if (w->omp_colidx[thread_id] >= w->omp_S[thread_id]->size2 - 15)
            {
              /* fold current matrix block into JTJ_vec, one thread at a time */
              gsl_matrix_view m = gsl_matrix_submatrix(w->omp_S[thread_id], 0, 0, w->nnm_crust, w->omp_colidx[thread_id]);

#pragma omp critical
              {
                gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &m.matrix, 1.0, JTJ_crust);
              }

              w->omp_colidx[thread_id] = 0;
            }
        } /* for (j = 0; j < mptr->n; ++j) */

      progress_bar(stderr, (i + 1.0) / (double) w->nsat, 70);
    } /* for (i = 0; i < w->nsat; ++i) */

  for (i = 0; i < w->max_threads; ++i)
    {
      if (w->omp_colidx[i] > 0)
        {
          /* accumulate final Green's functions */
          gsl_matrix_view m = gsl_matrix_submatrix(w->omp_S[i], 0, 0, w->nnm_crust, w->omp_colidx[i]);
          gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &m.matrix, 1.0, JTJ_crust);
        }
    }

  progress_bar(stderr, 1.0, 70);
  fprintf(stderr, "\n");

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  return s;
}

/*
mfield_nonlinear_precompute_mixed()
  Compute J_crust^T W J_core for vector residuals

Inputs: sqrt_weights - sqrt(wts), length nres
        JTJ_mixed    - (output) output matrix, p_crust-by-p_core
        w            - workspace

Notes:
1) Lower triangle of JTJ_mixed is filled
2) JTJ_mixed must be initialized to zero by calling function
*/

static int
mfield_nonlinear_precompute_mixed(const gsl_vector *sqrt_weights, gsl_matrix * JTJ_mixed,
                                  mfield_workspace *w)
{
  int s = 0;
  const size_t ncontrol = gsl_bspline2_ncontrol(w->gauss_spline_workspace_p[0]);
  const size_t order = gsl_bspline2_order(w->gauss_spline_workspace_p[0]);
  size_t ii, i, j;
  struct timeval tv0, tv1;

  /* compute J_crust^T W J_core
   *
   * This matrix can be divided into nnm_crust-by-nnm_core blocks:
   *
   * A_i = \sum_{k=1}^{nres} N_i(t_k) B_crust(r_k) B_core^T(r_k)
   *
   * We have outer loop over i, followed by threaded loop over the residuals,
   * accumulate the internal Greens functions into large block matrices T and U,
   * which are then folded into A_i with DSYRK when it is full.
   */

  fprintf(stderr, "mfield_nonlinear_precompute_mixed: computing J_crust^T W J_core...");
  fprintf(stderr, "\n");
  gettimeofday(&tv0, NULL);

  for (ii = 0; ii < ncontrol; ++ii)
    {
      gsl_matrix_view JTJ_block = gsl_matrix_submatrix(JTJ_mixed, 0, ii * w->nnm_core,
                                                       w->nnm_crust, w->nnm_core);

      for (i = 0; i < w->max_threads; ++i)
        w->omp_rowidx[i] = 0;

      for (i = 0; i < w->nsat; ++i)
        {
          magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
          for (j = 0; j < mptr->n; ++j)
            {
              int thread_id, flag;
              gsl_vector *N_gauss;
              size_t ridx, left, istart;
              double t, Nii, alpha;
              gsl_vector *VX, *VY, *VZ;
              gsl_vector_view vx, vy, vz;

              if (MAGDATA_Discarded(mptr->flags[j]))
                continue;

              if (!MAGDATA_FitMF(mptr->flags[j]))
                continue;

              if (!(mptr->flags[j] & (MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z)))
                continue;

              thread_id = omp_get_thread_num();
              ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */
              t = epoch2year(mptr->t[j]);
              N_gauss = w->gauss_spline_workspace_p[thread_id]->B;
              VX = VY = VZ = NULL;

              left = gsl_bspline2_find_interval(t, &flag, w->gauss_spline_workspace_p[thread_id]);
              if (flag != 0)
                continue;

              if (ii < (left + 1) - order || ii > left)
                continue; /* N_{ii}(t) = 0 */

              if (mptr->flags[j] & MAGDATA_FLG_X)
                {
                  vx = gsl_matrix_row(w->omp_T[thread_id], w->omp_rowidx[thread_id]++);
                  VX = &vx.vector;
                }

              if (mptr->flags[j] & MAGDATA_FLG_Y)
                {
                  vy = gsl_matrix_row(w->omp_T[thread_id], w->omp_rowidx[thread_id]++);
                  VY = &vy.vector;
                }

              if (mptr->flags[j] & MAGDATA_FLG_Z)
                {
                  vz = gsl_matrix_row(w->omp_T[thread_id], w->omp_rowidx[thread_id]++);
                  VZ = &vz.vector;
                }

              /* calculate internal Green's functions for both core and crust */
              green_calc_int2(mptr->r[j], mptr->theta[j], mptr->phi[j], VX, VY, VZ, w->green_array_p[thread_id]);

              /* calculate non-zero B-splines for the Gauss coefficients for this timestamp */
              gsl_bspline2_eval_basis_nonzero(t, N_gauss, &istart, w->gauss_spline_workspace_p[thread_id]);

              Nii = gsl_vector_get(N_gauss, ii - istart);
              alpha = sqrt(Nii);

              if (mptr->flags[j] & MAGDATA_FLG_X)
                {
                  double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
                  gsl_blas_dscal(alpha * sqrt_wj, VX);
                }

              if (mptr->flags[j] & MAGDATA_FLG_Y)
                {
                  double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
                  gsl_blas_dscal(alpha * sqrt_wj, VY);
                }

              if (mptr->flags[j] & MAGDATA_FLG_Z)
                {
                  double sqrt_wj = gsl_vector_get(sqrt_weights, ridx++);
                  gsl_blas_dscal(alpha * sqrt_wj, VZ);
                }

              /*
               * check if omp_T[thread_id] is full and should be folded into JTJ; the
               * 15 is just some slop to prevent trying to fill columns past the matrix buffer
               * in the loop above
               */
              if (w->omp_rowidx[thread_id] >= w->omp_T[thread_id]->size1 - 15)
                {
                  /* fold current matrix blocks into JTJ_block, one thread at a time */
                  gsl_matrix_view m_core = gsl_matrix_submatrix(w->omp_T[thread_id], 0, 0, w->omp_rowidx[thread_id], w->nnm_core);
                  gsl_matrix_view m_crust = gsl_matrix_submatrix(w->omp_T[thread_id], 0, w->nnm_core, w->omp_rowidx[thread_id], w->nnm_crust);

#pragma omp critical
                  {
                    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &m_crust.matrix, &m_core.matrix, 1.0, &JTJ_block.matrix);
                  }

                  w->omp_rowidx[thread_id] = 0;
                }
            } /* for (j = 0; j < mptr->n; ++j) */
        } /* for (i = 0; i < w->nsat; ++i) */

      for (i = 0; i < w->max_threads; ++i)
        {
          if (w->omp_rowidx[i] > 0)
            {
              /* accumulate final Green's functions into this ii block */
              gsl_matrix_view m_core = gsl_matrix_submatrix(w->omp_T[i], 0, 0, w->omp_rowidx[i], w->nnm_core);
              gsl_matrix_view m_crust = gsl_matrix_submatrix(w->omp_T[i], 0, w->nnm_core, w->omp_rowidx[i], w->nnm_crust);
              gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &m_crust.matrix, &m_core.matrix, 1.0, &JTJ_block.matrix);
            }
        }

      progress_bar(stderr, (double) ii / (double) ncontrol, 70);
    } /* for (ii = 0; ii < ncontrol; ++ii) */

  progress_bar(stderr, 1.0, 70);
  fprintf(stderr, "\n");

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  return s;
}

/*
mfield_nonlinear_precompute()
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
mfield_nonlinear_precompute(const gsl_vector *sqrt_weights, mfield_workspace *w)
{
  int s = GSL_SUCCESS;
  size_t nres_vec = 0;
  size_t i, j;

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

#if 1/*XXX*/
  if (w->p_core > 0)
    {
      gsl_matrix_view JTJ_core = gsl_matrix_submatrix(w->JTJ_vec, 0, 0, w->p_core, w->p_core);
      mfield_nonlinear_precompute_core(sqrt_weights, &JTJ_core.matrix, w);
    }
#endif

  if (w->p_crust > 0)
    {
      gsl_matrix_view JTJ_crust = gsl_matrix_submatrix(w->JTJ_vec, w->p_core, w->p_core, w->p_crust, w->p_crust);
#if 0
      mfield_nonlinear_precompute_crust(sqrt_weights, &JTJ_crust.matrix, w);
#else
      mfield_nonlinear_precompute_crust2(sqrt_weights, &JTJ_crust.matrix, w);
#endif
    }

  if (w->p_core > 0 && w->p_crust > 0)
    {
      gsl_matrix_view JTJ_mixed = gsl_matrix_submatrix(w->JTJ_vec, w->p_core, 0, w->p_crust, w->p_core);
      mfield_nonlinear_precompute_mixed(sqrt_weights, &JTJ_mixed.matrix, w);
    }

#if 0
  printsym_octave(w->JTJ_vec, "JTJ_vec");
  exit(1);
#endif

  return s;
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
  const mfield_parameters *mparams = &(w->params);
  size_t i, j;
  gsl_matrix_view JTJ_int; /* internal field portion of J^T J */
  gsl_matrix_view vJ2TJ1;   /* J_2^T J_1, p_sparse-by-p_int */
  gsl_matrix *J2TJ1 = NULL;
  size_t ndata_completed = 0;
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
      {
        gsl_vector_view tmp = gsl_vector_subvector(w->sqrt_wts_final, 0, w->nres);
        gsl_spmatrix_scale_rows(w->J2_csr, &tmp.vector);
      }
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
          gsl_vector_view vtmp = gsl_vector_subvector(v, w->p_int, w->p_sparse);
          gsl_vector_const_view utmp = gsl_vector_const_subvector(u, 0, w->nres);
          gsl_spblas_dusmv(CblasTrans, 1.0, w->J2_csr, &utmp.vector, 0.0, &vtmp.vector);
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

              /*XXX*/
              gsl_vector_set_zero(&Jv.vector);

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

#if MFIELD_DEBUG

#pragma omp atomic
          ndata_completed++;

#pragma omp critical
          if (ndata_completed % 1000 == 0)
            {
              mfield_debug("\t");
              progress_bar(stderr, (double) ndata_completed / (double) w->ndata, 70);
            }
#endif /* MFIELD_DEBUG */
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  /* accumulate any last rows of internal field Green's functions */
  for (i = 0; i < w->max_threads; ++i)
    {
      if (JTJ && w->omp_rowidx[i] > 0)
        {
          gsl_matrix_view Jm = gsl_matrix_submatrix(w->omp_J[i], 0, 0, w->omp_rowidx[i], w->p_int);
          gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &Jm.matrix, 1.0, &JTJ_int.matrix);
        }
    }

#if MFIELD_DEBUG
  mfield_debug("\t");
  progress_bar(stderr, (double) ndata_completed / (double) w->ndata, 70);
  mfield_debug("\n");
#endif

  if (mparams->regularize && !mparams->synth_data)
    {
      if (JTJ)
        {
          /* add Lambda to J^T J */
          gsl_spmatrix_add_to_dense(JTJ, w->Lambda);
        }

      if (TransJ == CblasTrans)
        {
          /* add L u(n+1:end) to v */
          gsl_vector_const_view tmp = gsl_vector_const_subvector(u, w->nres, w->p);
          gsl_spblas_dusmv(CblasNoTrans, 1.0, w->L, &tmp.vector, 1.0, v);
        }
    }

#if 0
  printv_octave(w->wts_final, "wts");

  if (u)
    printv_octave(u, "u");

  if (v)
    printv_octave(v, "JTu");

  if (JTJ)
    printsym_octave(JTJ, "JTJ");

  if (w->p_sparse > 0)
    printsp_octave(w->J2, "J2");

  exit(1);
#endif

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_df3: leaving function... (%g seconds)\n", time_diff(tv0, tv1));

  return GSL_SUCCESS;
}

/*
mfield_nonlinear_scalar_core()
  Compute J_core^T W J_core for scalar residuals

Inputs: sqrt_weights - sqrt(wts), length nres
        JTJ_core     - (output) output matrix, p_core-by-p_core
        w            - workspace

Notes:
1) Lower triangle of JTJ_core is filled
2) JTJ_core must be initialized to zero by calling function
*/

static int
mfield_nonlinear_scalar_core(const gsl_vector * x, const gsl_vector * sqrt_weights, gsl_matrix * JTJ_core,
                             mfield_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax_core;
  const size_t mmax = nmax;
  const size_t ncontrol = gsl_bspline2_ncontrol(w->gauss_spline_workspace_p[0]);
  const size_t order = gsl_bspline2_order(w->gauss_spline_workspace_p[0]);
  size_t ii, jj, i, j;
  struct timeval tv0, tv1;

  /* compute J_core^T W J_core
   *
   * This matrix can be divided into nnm_core-by-nnm_core blocks:
   *
   * A_{ij} = \sum_{k=1}^{nres} N_i(t_k) N_j(t_k) B(r_k) B^T(r_k)
   *
   * We have outer loops over i,j, followed by threaded loop over the residuals,
   * accumulate the internal Greens functions into a large block matrix T,
   * which is then folded into A_{ij} with DSYRK when it is full.
   */

  fprintf(stderr, "mfield_nonlinear_scalar_core: computing J_core^T W J_core...");
  fprintf(stderr, "\n");
  gettimeofday(&tv0, NULL);

  for (ii = 0; ii < ncontrol; ++ii)
    {
      for (jj = 0; jj <= ii; ++jj)
        {
          gsl_matrix_view JTJ_block = gsl_matrix_submatrix(JTJ_core, ii * w->nnm_core, jj * w->nnm_core,
                                                           w->nnm_core, w->nnm_core);

          for (i = 0; i < w->max_threads; ++i)
            w->omp_colidx[i] = 0;

          for (i = 0; i < w->nsat; ++i)
            {
              magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
              for (j = 0; j < mptr->n; ++j)
                {
                  int thread_id, flag;
                  gsl_vector *N_gauss;
                  size_t ridx, left, istart;
                  double t, Nii, Njj, alpha;
                  gsl_vector_view vx, vy, vz;

                  if (MAGDATA_Discarded(mptr->flags[j]))
                    continue;

                  if (!MAGDATA_FitMF(mptr->flags[j]))
                    continue;

                  if (!MAGDATA_ExistScalar(mptr->flags[j]))
                    continue;

                  thread_id = omp_get_thread_num();
                  ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */
                  N_gauss = w->gauss_spline_workspace_p[thread_id]->B;
                  vx = gsl_matrix_subrow(w->omp_dB[thread_id], 0, 0, w->nnm_core);
                  vy = gsl_matrix_subrow(w->omp_dB[thread_id], 1, 0, w->nnm_core);
                  vz = gsl_matrix_subrow(w->omp_dB[thread_id], 2, 0, w->nnm_core);
                  t = epoch2year(mptr->t[j]);

                  left = gsl_bspline2_find_interval(t, &flag, w->gauss_spline_workspace_p[thread_id]);
                  if (flag != 0)
                    continue;

                  if (ii < (left + 1) - order || ii > left)
                    continue; /* N_{ii}(t) = 0 */
                  if (jj < (left + 1) - order || jj > left)
                    continue; /* N_{ii}(t) = 0 */

                  /* calculate internal Green's functions for core field */
                  green_calc_intopt(1, nmax, mmax, mptr->r[j], mptr->theta[j], mptr->phi[j],
                                    &vx.vector, &vy.vector, &vz.vector, w->green_array_p[thread_id]);

                  /* calculate non-zero B-splines for the Gauss coefficients for this timestamp */
                  gsl_bspline2_eval_basis_nonzero(t, N_gauss, &istart, w->gauss_spline_workspace_p[thread_id]);

                  Nii = gsl_vector_get(N_gauss, ii - istart);
                  Njj = gsl_vector_get(N_gauss, jj - istart);
                  alpha = sqrt(Nii * Njj);

                  if (mptr->flags[j] & MAGDATA_FLG_X)
                    {
                      ++ridx;
                    }

                  if (mptr->flags[j] & MAGDATA_FLG_Y)
                    {
                      ++ridx;
                    }

                  if (mptr->flags[j] & MAGDATA_FLG_Z)
                    {
                      ++ridx;
                    }

                  if (MAGDATA_ExistScalar(mptr->flags[j]))
                    {
                      double sqrt_wj = gsl_vector_get(sqrt_weights, ridx);
                      double b[3], F;
                      gsl_vector_view bv = gsl_vector_view_array(b, 3);
                      gsl_vector_view v = gsl_matrix_subcolumn(w->omp_T[thread_id], w->omp_colidx[thread_id]++, 0, w->nnm_core);
                      gsl_matrix_view m = gsl_matrix_submatrix(w->omp_dB[thread_id], 0, 0, 3, w->nnm_core);

                      /* compute internal field model */
                      b[0] = mfield_nonlinear_model_int(N_gauss, istart, &vx.vector, x, w);
                      b[1] = mfield_nonlinear_model_int(N_gauss, istart, &vy.vector, x, w);
                      b[2] = mfield_nonlinear_model_int(N_gauss, istart, &vz.vector, x, w);

                      /* add apriori model of external (and possibly crustal) field */
                      b[0] += mptr->Bx_model[j];
                      b[1] += mptr->By_model[j];
                      b[2] += mptr->Bz_model[j];

                      F = gsl_hypot3(b[0], b[1], b[2]);

                      b[0] /= F;
                      b[1] /= F;
                      b[2] /= F;

                      gsl_blas_dgemv(CblasTrans, sqrt_wj * alpha, &m.matrix, &bv.vector, 0.0, &v.vector);

                      ++ridx;
                    }

                  /*
                   * check if omp_T[thread_id] is full and should be folded into JTJ; the
                   * 15 is just some slop to prevent trying to fill columns past the matrix buffer
                   * in the loop above
                   */
                  if (w->omp_colidx[thread_id] >= w->omp_T[thread_id]->size2 - 15)
                    {
                      /* fold current matrix block into JTJ_vec, one thread at a time */
                      gsl_matrix_view m = gsl_matrix_submatrix(w->omp_T[thread_id], 0, 0, w->nnm_core, w->omp_colidx[thread_id]);

#pragma omp critical
                      {
                        gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &m.matrix, 1.0, &JTJ_block.matrix);
                      }

                      w->omp_colidx[thread_id] = 0;
                    }
                } /* for (j = 0; j < mptr->n; ++j) */
            } /* for (i = 0; i < w->nsat; ++i) */

          for (i = 0; i < w->max_threads; ++i)
            {
              if (w->omp_colidx[i] > 0)
                {
                  /* accumulate final Green's functions into this (ii,jj) block */
                  gsl_matrix_view m = gsl_matrix_submatrix(w->omp_T[i], 0, 0, w->nnm_core, w->omp_colidx[i]);
                  gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &m.matrix, 1.0, &JTJ_block.matrix);
                }
            }

          /* at this point the (ii, jj) block of JTJ is calculated, but if it is not
           * a diagonal block, we need to copy the lower triangle to the upper */
          if (ii != jj)
            {
              gsl_matrix_transpose_tricpy('L', 0, &JTJ_block.matrix, &JTJ_block.matrix);
            }

        } /* for (jj = 0; jj <= ii; ++jj) */

      progress_bar(stderr, (double) ii / (double) ncontrol, 70);
    } /* for (ii = 0; ii < ncontrol; ++ii) */

  progress_bar(stderr, 1.0, 70);
  fprintf(stderr, "\n");

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  return s;
}

/*
mfield_nonlinear_scalar_crust()
  Compute J_crust^T W J_crust for scalar residuals

Inputs: sqrt_weights - sqrt(wts), length nres
        JTJ_crust    - (output) output matrix, p_crust-by-p_crust
        w            - workspace

Notes:
1) Lower triangle of JTJ_crust is filled
2) JTJ_crust must be initialized to zero by calling function
*/

static int
mfield_nonlinear_scalar_crust(const gsl_vector * x, const gsl_vector * sqrt_weights, gsl_matrix * JTJ_crust,
                              mfield_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = nmax;
  size_t i, j;
  struct timeval tv0, tv1;

  /* compute J_crust^T W J_crust
   *
   * JTJ_crust = \sum_{k=1}^{nres} B^T(r_k) b(r_k,t_k) b^T(r_k,t_k) B(r_k)
   */

  fprintf(stderr, "mfield_nonlinear_scalar_crust: computing J_crust^T W J_crust...");
  fprintf(stderr, "\n");
  gettimeofday(&tv0, NULL);

  for (i = 0; i < w->max_threads; ++i)
    w->omp_colidx[i] = 0;

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          size_t ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (!MAGDATA_FitMF(mptr->flags[j]))
            continue;

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              ++ridx;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              ++ridx;
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]))
            {
              double t = epoch2year(mptr->t[j]);
              double sqrt_wj = gsl_vector_get(sqrt_weights, ridx);
              gsl_vector_view vx = gsl_matrix_row(w->omp_dB[thread_id], 0);
              gsl_vector_view vy = gsl_matrix_row(w->omp_dB[thread_id], 1);
              gsl_vector_view vz = gsl_matrix_row(w->omp_dB[thread_id], 2);
              double b[3], F;
              gsl_vector_view bv = gsl_vector_view_array(b, 3);
              gsl_vector_view v = gsl_matrix_subcolumn(w->omp_T[thread_id], w->omp_colidx[thread_id]++, 0, w->nnm_crust);
              gsl_matrix_view m = gsl_matrix_submatrix(w->omp_dB[thread_id], 0, w->nnm_core, 3, w->nnm_crust);
              gsl_vector *N_gauss = w->gauss_spline_workspace_p[thread_id]->B;
              size_t istart;

              /* calculate non-zero B-splines for the Gauss coefficients for this timestamp */
              gsl_bspline2_eval_basis_nonzero(t, N_gauss, &istart, w->gauss_spline_workspace_p[thread_id]);

              /* calculate internal Green's functions for core and crustal field */
              green_calc_intopt(1, nmax, mmax, mptr->r[j], mptr->theta[j], mptr->phi[j],
                                &vx.vector, &vy.vector, &vz.vector, w->green_array_p[thread_id]);

              /* compute internal field model */
              b[0] = mfield_nonlinear_model_int(N_gauss, istart, &vx.vector, x, w);
              b[1] = mfield_nonlinear_model_int(N_gauss, istart, &vy.vector, x, w);
              b[2] = mfield_nonlinear_model_int(N_gauss, istart, &vz.vector, x, w);

              /* add apriori model of external (and possibly crustal) field */
              b[0] += mptr->Bx_model[j];
              b[1] += mptr->By_model[j];
              b[2] += mptr->Bz_model[j];

              F = gsl_hypot3(b[0], b[1], b[2]);

              b[0] /= F;
              b[1] /= F;
              b[2] /= F;

              gsl_blas_dgemv(CblasTrans, sqrt_wj, &m.matrix, &bv.vector, 0.0, &v.vector);

              ++ridx;
            }

          /*
           * check if omp_T[thread_id] is full and should be folded into JTJ; the
           * 15 is just some slop to prevent trying to fill columns past the matrix buffer
           * in the loop above
           */
          if (w->omp_colidx[thread_id] >= w->omp_T[thread_id]->size2 - 15)
            {
              /* fold current matrix block into JTJ_vec, one thread at a time */
              gsl_matrix_view m = gsl_matrix_submatrix(w->omp_T[thread_id], 0, 0, w->nnm_crust, w->omp_colidx[thread_id]);

#pragma omp critical
              {
                gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &m.matrix, 1.0, JTJ_crust);
              }

              w->omp_colidx[thread_id] = 0;
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  for (i = 0; i < w->max_threads; ++i)
    {
      if (w->omp_colidx[i] > 0)
        {
          /* accumulate final Green's functions into this (ii,jj) block */
          gsl_matrix_view m = gsl_matrix_submatrix(w->omp_T[i], 0, 0, w->nnm_crust, w->omp_colidx[i]);
          gsl_blas_dsyrk(CblasLower, CblasNoTrans, 1.0, &m.matrix, 1.0, JTJ_crust);
        }
    }

  progress_bar(stderr, 1.0, 70);
  fprintf(stderr, "\n");

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  return s;
}

/*
mfield_nonlinear_scalar_mixed()
  Compute J_crust^T W J_core for scalar residuals

Inputs: sqrt_weights - sqrt(wts), length nres
        JTJ_mixed    - (output) output matrix, p_crust-by-p_core
        w            - workspace

Notes:
1) Lower triangle of JTJ_mixed is filled
2) JTJ_mixed must be initialized to zero by calling function
*/

static int
mfield_nonlinear_scalar_mixed(const gsl_vector * x, const gsl_vector * sqrt_weights, gsl_matrix * JTJ_mixed,
                              mfield_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = nmax;
  const size_t ncontrol = gsl_bspline2_ncontrol(w->gauss_spline_workspace_p[0]);
  const size_t order = gsl_bspline2_order(w->gauss_spline_workspace_p[0]);
  size_t ii, i, j;
  struct timeval tv0, tv1;

  /* compute J_crust^T W J_core
   *
   * This matrix can be divided into nnm_crust-by-nnm_core blocks:
   *
   * A_i = \sum_{k=1}^{nres} N_i(t_k) B^T_crust(r_k) b(r_k,t_k) b^T(r_k,t_k) B_core(r_k)
   *
   * We have outer loop over i, followed by threaded loop over the residuals,
   * accumulate the internal Greens functions into a large block matrix T,
   * which is then folded into A_i with DGEMM when it is full.
   */

  fprintf(stderr, "mfield_nonlinear_scalar_mixed: computing J_crust^T W J_core...");
  fprintf(stderr, "\n");
  gettimeofday(&tv0, NULL);

  for (ii = 0; ii < ncontrol; ++ii)
    {
      gsl_matrix_view JTJ_block = gsl_matrix_submatrix(JTJ_mixed, 0, ii * w->nnm_core, w->nnm_crust, w->nnm_core);

      for (i = 0; i < w->max_threads; ++i)
        w->omp_colidx[i] = 0;

      for (i = 0; i < w->nsat; ++i)
        {
          magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
          for (j = 0; j < mptr->n; ++j)
            {
              int thread_id, flag;
              double t;
              size_t ridx, left;

              if (MAGDATA_Discarded(mptr->flags[j]))
                continue;

              if (!MAGDATA_FitMF(mptr->flags[j]))
                continue;

              if (!MAGDATA_ExistScalar(mptr->flags[j]))
                continue;

              thread_id = omp_get_thread_num();
              ridx = mptr->index[j]; /* residual index for this data point in [0:nres-1] */
              t = epoch2year(mptr->t[j]);

              left = gsl_bspline2_find_interval(t, &flag, w->gauss_spline_workspace_p[thread_id]);
              if (flag != 0)
                continue;

              if (ii < (left + 1) - order || ii > left)
                continue; /* N_{ii}(t) = 0 */

              if (mptr->flags[j] & MAGDATA_FLG_X)
                {
                  ++ridx;
                }

              if (mptr->flags[j] & MAGDATA_FLG_Y)
                {
                  ++ridx;
                }

              if (mptr->flags[j] & MAGDATA_FLG_Z)
                {
                  ++ridx;
                }

              if (MAGDATA_ExistScalar(mptr->flags[j]))
                {
                  double sqrt_wj = gsl_vector_get(sqrt_weights, ridx);
                  gsl_vector_view vx = gsl_matrix_row(w->omp_dB[thread_id], 0);
                  gsl_vector_view vy = gsl_matrix_row(w->omp_dB[thread_id], 1);
                  gsl_vector_view vz = gsl_matrix_row(w->omp_dB[thread_id], 2);
                  double b[3], F;
                  gsl_vector_view bv = gsl_vector_view_array(b, 3);
                  gsl_vector_view v = gsl_matrix_subcolumn(w->omp_T[thread_id], w->omp_colidx[thread_id]++, 0, w->nnm_tot);
                  gsl_vector *N_gauss = w->gauss_spline_workspace_p[thread_id]->B;
                  double Nii, alpha;
                  size_t istart;

                  /* calculate internal Green's functions for core and crustal field */
                  green_calc_intopt(1, nmax, mmax, mptr->r[j], mptr->theta[j], mptr->phi[j],
                                    &vx.vector, &vy.vector, &vz.vector, w->green_array_p[thread_id]);

                  /* calculate non-zero B-splines for the Gauss coefficients for this timestamp */
                  gsl_bspline2_eval_basis_nonzero(t, N_gauss, &istart, w->gauss_spline_workspace_p[thread_id]);

                  Nii = gsl_vector_get(N_gauss, ii - istart);
                  alpha = sqrt(Nii);

                  /* compute internal field model */
                  b[0] = mfield_nonlinear_model_int(N_gauss, istart, &vx.vector, x, w);
                  b[1] = mfield_nonlinear_model_int(N_gauss, istart, &vy.vector, x, w);
                  b[2] = mfield_nonlinear_model_int(N_gauss, istart, &vz.vector, x, w);

                  /* add apriori model of external (and possibly crustal) field */
                  b[0] += mptr->Bx_model[j];
                  b[1] += mptr->By_model[j];
                  b[2] += mptr->Bz_model[j];

                  F = gsl_hypot3(b[0], b[1], b[2]);

                  b[0] /= F;
                  b[1] /= F;
                  b[2] /= F;

                  gsl_blas_dgemv(CblasTrans, sqrt_wj * alpha, w->omp_dB[thread_id], &bv.vector, 0.0, &v.vector);

                  ++ridx;
                }

              /*
               * check if omp_T[thread_id] is full and should be folded into JTJ; the
               * 15 is just some slop to prevent trying to fill columns past the matrix buffer
               * in the loop above
               */
              if (w->omp_colidx[thread_id] >= w->omp_T[thread_id]->size2 - 15)
                {
                  /* fold current matrix block into JTJ_vec, one thread at a time */
                  gsl_matrix_view m_core = gsl_matrix_submatrix(w->omp_T[thread_id], 0, 0, w->nnm_core, w->omp_colidx[thread_id]);
                  gsl_matrix_view m_crust = gsl_matrix_submatrix(w->omp_T[thread_id], w->nnm_core, 0, w->nnm_crust, w->omp_colidx[thread_id]);

#pragma omp critical
                  {
                    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &m_crust.matrix, &m_core.matrix, 1.0, &JTJ_block.matrix);
                  }

                  w->omp_colidx[thread_id] = 0;
                }
            } /* for (j = 0; j < mptr->n; ++j) */
        } /* for (i = 0; i < w->nsat; ++i) */

      for (i = 0; i < w->max_threads; ++i)
        {
          if (w->omp_colidx[i] > 0)
            {
              /* accumulate final Green's functions into this (ii,jj) block */
              gsl_matrix_view m_core = gsl_matrix_submatrix(w->omp_T[i], 0, 0, w->nnm_core, w->omp_colidx[i]);
              gsl_matrix_view m_crust = gsl_matrix_submatrix(w->omp_T[i], w->nnm_core, 0, w->nnm_crust, w->omp_colidx[i]);
              gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &m_crust.matrix, &m_core.matrix, 1.0, &JTJ_block.matrix);
            }
        }

      progress_bar(stderr, (double) ii / (double) ncontrol, 70);
    } /* for (ii = 0; ii < ncontrol; ++ii) */

  progress_bar(stderr, 1.0, 70);
  fprintf(stderr, "\n");

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  return s;
}

static int
mfield_calc_df4(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u,
                void *params, gsl_vector * v, gsl_matrix * JTJ)
{
  mfield_workspace *w = (mfield_workspace *) params;
  size_t i, j;
  gsl_matrix_view JTJ_int; /* internal field portion of J^T J */
  gsl_matrix_view vJ2TJ1;   /* J_2^T J_1, p_sparse-by-p_int */
  gsl_matrix *J2TJ1 = NULL;
  size_t ndata_completed = 0;
  struct timeval tv0, tv1, tv2, tv3;

  gettimeofday(&tv0, NULL);
  mfield_debug("mfield_calc_df4: entering function...\n");
  mfield_debug("mfield_calc_df4: TransJ = %s...\n",
               TransJ == CblasTrans ? "trans" : "notrans");

  /* initialize outputs to 0 */
  if (v)
    gsl_vector_set_zero(v);

  if (JTJ)
    {
      gsl_matrix_view vJTJ;

      gsl_matrix_set_zero(JTJ);

      /* copy previously computed vector internal field portion of J^T J
       * (doesn't depend on x) */
      JTJ_int = gsl_matrix_submatrix(JTJ, 0, 0, w->p_int, w->p_int);
      gsl_matrix_tricpy('L', 1, &JTJ_int.matrix, w->JTJ_vec);

      vJTJ = gsl_matrix_submatrix(JTJ, 0, 0, w->p_core, w->p_core);
      mfield_nonlinear_scalar_core(x, w->sqrt_wts_final, &vJTJ.matrix, w);

      vJTJ = gsl_matrix_submatrix(JTJ, w->p_core, w->p_core, w->p_crust, w->p_crust);
      mfield_nonlinear_scalar_crust(x, w->sqrt_wts_final, &vJTJ.matrix, w);

      vJTJ = gsl_matrix_submatrix(JTJ, w->p_core, 0, w->p_crust, w->p_core);
      mfield_nonlinear_scalar_mixed(x, w->sqrt_wts_final, &vJTJ.matrix, w);
    }

  if (w->p_sparse > 0)
    {
      /* calculate sparse part of Jacobian J_2 */
      mfield_calc_J2(x, w->J2, w);

      if (w->J2_csr == NULL)
        w->J2_csr = gsl_spmatrix_alloc_nzmax(w->J2->size1, w->J2->size2, w->J2->nz, GSL_SPMATRIX_CSR);

      mfield_debug("mfield_calc_df4: compressing J_2 to CSR...");
      gettimeofday(&tv2, NULL);
      gsl_spmatrix_csr(w->J2_csr, w->J2);
      gettimeofday(&tv3, NULL);
      mfield_debug("done (%g seconds)\n", time_diff(tv2, tv3));

      mfield_debug("mfield_calc_df4: computing sqrt(W) J_2...");
      gettimeofday(&tv2, NULL);
      gsl_spmatrix_scale_rows(w->J2_csr, w->sqrt_wts_final);
      gettimeofday(&tv3, NULL);
      mfield_debug("done (%g seconds)\n", time_diff(tv2, tv3));

      if (JTJ)
        {
          /* store J_2^T J_2 in lower right portion of JTJ */
          gsl_matrix_view J2TJ2 = gsl_matrix_submatrix(JTJ, w->p_int, w->p_int, w->p_sparse, w->p_sparse);

          mfield_debug("mfield_calc_df4: computing J_2^T W J_2...");
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
          double t = mptr->ts[j];       /* use scaled time */
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          gsl_vector *N_gauss = w->gauss_spline_workspace_p[thread_id]->B;
          size_t istart_gauss;

          gsl_vector_view vx = gsl_matrix_row(w->omp_dB[thread_id], 0);
          gsl_vector_view vy = gsl_matrix_row(w->omp_dB[thread_id], 1);
          gsl_vector_view vz = gsl_matrix_row(w->omp_dB[thread_id], 2);

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
              double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx);
              double uj = gsl_vector_get(u, ridx);
              double b[3], F;
              gsl_vector_view bv = gsl_vector_view_array(b, 3);
              gsl_vector_view tmp = gsl_matrix_row(w->omp_dF, thread_id);

              /* compute internal field model */
              b[0] = mfield_nonlinear_model_int(N_gauss, istart_gauss, &vx.vector, x, w);
              b[1] = mfield_nonlinear_model_int(N_gauss, istart_gauss, &vy.vector, x, w);
              b[2] = mfield_nonlinear_model_int(N_gauss, istart_gauss, &vz.vector, x, w);

              /* add apriori model of external (and possibly crustal) field */
              b[0] += mptr->Bx_model[j];
              b[1] += mptr->By_model[j];
              b[2] += mptr->Bz_model[j];

              F = gsl_hypot3(b[0], b[1], b[2]);

              b[0] /= F;
              b[1] /= F;
              b[2] /= F;

              gsl_blas_dgemv(CblasTrans, sqrt_wj, w->omp_dB[thread_id], &bv.vector, 0.0, &tmp.vector);

#pragma omp critical
              {
                mfield_jacobian_J1Tu_scalar(N_gauss, istart_gauss, uj, &tmp.vector, v, w);
              }

#if 0
#pragma omp critical
              {
                mfield_jacobian_F(TransJ, istart_gauss, u, ridx, w->J2_csr, &Jv.vector, v, J2TJ1, w);
              }
#endif

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

#if MFIELD_DEBUG

#pragma omp atomic
          ndata_completed++;

#pragma omp critical
          if (ndata_completed % 1000 == 0)
            {
              progress_bar(stderr, (double) ndata_completed / (double) w->ndata, 70);
            }
#endif /* MFIELD_DEBUG */
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  /* accumulate any last rows of internal field Green's functions */
  for (i = 0; i < w->max_threads; ++i)
    {
      if (JTJ && w->omp_rowidx[i] > 0)
        {
          gsl_matrix_view Jm = gsl_matrix_submatrix(w->omp_J[i], 0, 0, w->omp_rowidx[i], w->p_int);
          gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &Jm.matrix, 1.0, &JTJ_int.matrix);
        }
    }

#if MFIELD_DEBUG
  progress_bar(stderr, (double) ndata_completed / (double) w->ndata, 70);
  mfield_debug("\n");
#endif

  /* regularize by adding Lambda to J^T J */
  if (JTJ)
    {
      gsl_spmatrix_add_to_dense(JTJ, w->Lambda);
    }

#if 0
  printv_octave(w->wts_final, "wts");

  if (u)
    printv_octave(u, "u");

  if (v)
    printv_octave(v, "JTu");

  if (JTJ)
    printsym_octave(JTJ, "JTJ");

  exit(1);
#endif

  gettimeofday(&tv1, NULL);
  mfield_debug("mfield_calc_df4: leaving function... (%g seconds)\n", time_diff(tv0, tv1));

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
mfield_jacobian_J1Tu_scalar()
  Update the J_1^T sqrt(W) u vector with a new row of the Jacobian matrix,
corresponding to a scalar residual.

J^T sqrt(W) u = [   J_1^T sqrt(W) u  ]
                [ J_2^(x)T sqrt(W) u ]

Inputs: N       - B-spline functions for t(j)
        istart  - index of first non-zero B-spline in [0, ncontrol-1]
        uj      - input vector element u(j)
        dB      - B^T(r_k) b(r_k, t_k) vector, length nnm_tot
        JTu     - (output) J_1^T sqrt(W) u vector, length p_int
        w       - workspace
*/

static inline int
mfield_jacobian_J1Tu_scalar(const gsl_vector * N, const size_t istart, const double uj,
                            const gsl_vector * dB, gsl_vector *JTu, const mfield_workspace *w)
{
  gsl_vector_const_view dB_core = gsl_vector_const_subvector(dB, 0, w->nnm_core);
  size_t k;

  /* update J^T u */
  for (k = 0; k < w->params.gauss_spline_order; ++k)
    {
      double Nk = gsl_vector_get(N, k);
      gsl_vector_view v = gsl_vector_subvector(JTu, (k + istart) * w->nnm_core, w->nnm_core);
      gsl_blas_daxpy(-Nk * uj, &dB_core.vector, &v.vector);
    }

  if (w->nnm_crust > 0)
    {
      gsl_vector_const_view dB_crust = gsl_vector_const_subvector(dB, w->nnm_core, w->nnm_crust);
      gsl_vector_view v = gsl_vector_subvector(JTu, w->p_core, w->nnm_crust);
      gsl_blas_daxpy(-uj, &dB_crust.vector, &v.vector);
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
