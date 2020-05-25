#include <mainlib/ml_oct.h>

static int invert_calc_nonlinear_multilarge(const gsl_vector * c, invert_workspace * w);
static int invert_accumulate(const gsl_vector *x, const gsl_vector * sqrt_wts, invert_workspace * w);

/*
invert_calc_nonlinear_multilarge()
  Calculate a solution to inverse problem using multilarge

Inputs: c - coefficient vector
        w - workspace

Notes:
1) w->wts_final must be initialized prior to calling this function
2) On output, w->c contains the solution vector
*/

static int
invert_calc_nonlinear_multilarge(const gsl_vector * c, invert_workspace * w)
{
  int s = 0;
  struct timeval tv0, tv1;

  if (w->lls_solution == 1)
    {
      double rnorm, snorm, rcond;

      invert_accumulate(c, w->sqrt_wts_final, w);

      fprintf(stderr, "invert_calc_nonlinear_multilarge: computing condition number of LS matrix...");
      gettimeofday(&tv0, NULL);
      gsl_multilarge_linear_rcond(&rcond, w->multilarge_linear_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, condition number = %.12e)\n", time_diff(tv0, tv1), 1.0 / rcond);

      fprintf(stderr, "invert_calc_nonlinear_multilarge: solving LS system...");
      gettimeofday(&tv0, NULL);
      gsl_multilarge_linear_solve(0.0, w->c, &rnorm, &snorm, w->multilarge_linear_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, rnorm = %e, snorm = %e)\n", time_diff(tv0, tv1), rnorm, snorm);

      printv_octave(w->c, "c");
    }
  else
    {
    }

  return s;
}

/*
invert_accumulate()
  Compute Jacobian matrix J(x) using OpenMP to
calculate Green's functions quickly.

Inputs: x        - parameter vector, length p
        sqrt_wts - sqrt(W) vector
        w        - workspace
*/

static int
invert_accumulate(const gsl_vector *x, const gsl_vector * sqrt_wts, invert_workspace * w)
{
  int s = GSL_SUCCESS;
  const invert_parameters *mparams = &(w->params);
  size_t i, j;
  struct timeval tv0, tv1;
  size_t ndata_completed = 0;

  invert_debug("invert_accumulate: entering function...\n");
  gettimeofday(&tv0, NULL);

  gsl_multilarge_linear_reset(w->multilarge_linear_p);

  /*
   * omp_rowidx[thread_id] contains the number of currently filled rows
   * of omp_J[thread_id]. When omp_J[thread_id] is filled, it is folded
   * into the LS system and then omp_rowidx[thread_id] is reset to 0
   */
  for (i = 0; i < w->max_threads; ++i)
    w->omp_rowidx[i] = 0;

  invert_debug("invert_accumulate: building LS system...\n");

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = invert_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          gsl_vector_view vx, vy, vz;
          gsl_vector *VX = NULL, *VY = NULL, *VZ = NULL;
          double *fx = NULL, *fy = NULL, *fz = NULL;
          double B_obs[3],   /* observed NEC vector */
                 B_prior[3], /* prior model */
                 B_model[3]; /* fitted model */

          /* use supplied NEC vector */
          B_obs[0] = mptr->Bx_nec[j];
          B_obs[1] = mptr->By_nec[j];
          B_obs[2] = mptr->Bz_nec[j];

          B_prior[0] = mptr->Bx_model[j];
          B_prior[1] = mptr->By_model[j];
          B_prior[2] = mptr->Bz_model[j];

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              fx = gsl_vector_ptr(w->omp_f[thread_id], w->omp_rowidx[thread_id]);
              vx = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
              VX = &vx.vector;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              fy = gsl_vector_ptr(w->omp_f[thread_id], w->omp_rowidx[thread_id]);
              vy = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
              VY = &vy.vector;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              fz = gsl_vector_ptr(w->omp_f[thread_id], w->omp_rowidx[thread_id]);
              vz = gsl_matrix_row(w->omp_J[thread_id], w->omp_rowidx[thread_id]++);
              VZ = &vz.vector;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z))
            {
              invert_jacobian_vector(thread_id, mptr->t[j], r, theta, phi, VX, VY, VZ, w);
            }

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx++);
              *fx = sqrt_wj * (B_obs[0] - B_prior[0]);
              gsl_blas_dscal(sqrt_wj, VX);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx++);
              *fy = sqrt_wj * (B_obs[1] - B_prior[1]);
              gsl_blas_dscal(sqrt_wj, VY);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              double sqrt_wj = gsl_vector_get(w->sqrt_wts_final, ridx++);
              *fz = sqrt_wj * (B_obs[2] - B_prior[2]);
              gsl_blas_dscal(sqrt_wj, VZ);
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) &&
              MAGDATA_FitMF(mptr->flags[j]))
            {
              ++ridx;
            }

          if (w->omp_rowidx[thread_id] > w->omp_J[thread_id]->size1 - 4)
            {
              gsl_matrix_view m = gsl_matrix_submatrix(w->omp_J[thread_id], 0, 0, w->omp_rowidx[thread_id], w->p);
              gsl_vector_view f = gsl_vector_subvector(w->omp_f[thread_id], 0, w->omp_rowidx[thread_id]);

#pragma omp critical
              {
                gsl_multilarge_linear_accumulate(&m.matrix, &f.vector, w->multilarge_linear_p);
              }

              w->omp_rowidx[thread_id] = 0;
            }

#pragma omp atomic
          ndata_completed++;

#pragma omp critical
          if (ndata_completed % 1000 == 0)
            {
              progress_bar(stderr, (double) ndata_completed / (double) w->ndata, 70);
            }
        } /* for (j = 0; j < mptr->n; ++j) */
    } /* for (i = 0; i < w->nsat; ++i) */

  /* accumulate final partial blocks */
  for (i = 0; i < w->max_threads; ++i)
    {
      if (w->omp_rowidx[i] > 0)
        {
          gsl_matrix_view m = gsl_matrix_submatrix(w->omp_J[i], 0, 0, w->omp_rowidx[i], w->p);
          gsl_vector_view f = gsl_vector_subvector(w->omp_f[i], 0, w->omp_rowidx[i]);
          gsl_multilarge_linear_accumulate(&m.matrix, &f.vector, w->multilarge_linear_p);
        }
    }

  progress_bar(stderr, 1.0, 70);

  invert_debug("done\n");

  gettimeofday(&tv1, NULL);
  invert_debug("invert_accumulate: leaving function (%g seconds)\n", time_diff(tv0, tv1));

#if 0
  if (mparams->regularize && !mparams->synth_data)
    {
      /* copy L^T into lower portion of J */
      gsl_matrix_view m = gsl_matrix_submatrix(J, w->nres, 0, w->p, w->p);

      for (i = 0; i < w->L->nz; ++i)
        gsl_matrix_set(&m.matrix, w->L->p[i], w->L->i[i], w->L->data[i]);
    }
#endif

#if 0
  {
    gsl_matrix_view m = gsl_matrix_submatrix(w->omp_J[0], 0, 0, w->omp_rowidx[0], w->p);
    gsl_vector_view f = gsl_vector_subvector(w->omp_f[0], 0, w->omp_rowidx[0]);
    printv_octave(x, "x");
    printv_octave(&f.vector, "f2");
    print_octave(&m.matrix, "J2");
    exit(1);
  }
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

        gsl_matrix_free(JTJ);
        exit(1);
      }
  }
#endif

  return s;
}
