/*
 * invert_multifit.c
 *
 * Contains routines for fitting magnetic field module using gsl_multifit_nlinear framework
 */

#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <omp.h>

#include <mainlib/ml_oct.h>

static int invert_calc_nonlinear_multifit(const gsl_vector * c, invert_workspace * w);
static int invert_calc_f(const gsl_vector *x, void *params, gsl_vector *f);
static int invert_calc_Wf(const gsl_vector *x, void *params, gsl_vector *f);
static int invert_calc_df(const gsl_vector *x, void *params, gsl_matrix *J);
static int invert_nonlinear_model(const gsl_vector * x, const double t, const double r, const double theta, const double phi,
                                  const int thread_id, double B_model[3], invert_workspace *w);
static int invert_jacobian_vector(const int thread_id, const double t, const double r, const double theta, const double phi,
                                  gsl_vector * X, gsl_vector * Y, gsl_vector * Z,
                                  const invert_workspace *w);
static int invert_nonlinear_model_J(const gsl_vector * x, const double t, const double r,
                                    const double theta, const double phi,
                                    const int thread_id, double J_model[3], invert_workspace *w);
static int invert_jacobian_vector_J(const int thread_id, const double t, const double r, const double theta, const double phi,
                                    gsl_vector * X, gsl_vector * Y, gsl_vector * Z,
                                    const invert_workspace *w);

/*
invert_calc_nonlinear_multifit()
  Calculate a solution to inverse problem using multifit

Inputs: c - coefficient vector
        w - workspace

Notes:
1) w->wts_final must be initialized prior to calling this function
2) On output, w->c contains the solution vector
*/

static int
invert_calc_nonlinear_multifit(const gsl_vector * c, invert_workspace * w)
{
  int s = 0;
  const size_t p = w->p;          /* number of coefficients */
  const size_t n = w->nres_tot;   /* number of residuals */
  struct timeval tv0, tv1;

  if (w->lls_solution == 1)
    {
      gsl_matrix * J = w->multifit_nlinear_p->J;
      gsl_vector * f = w->multifit_nlinear_p->f;
      gsl_matrix * T = gsl_matrix_alloc(p, p);
      gsl_vector * x = gsl_vector_alloc(n);
      gsl_vector * work = gsl_vector_alloc(p);
      gsl_vector * work3 = gsl_vector_alloc(3 * p);
      gsl_vector_view x1 = gsl_vector_subvector(x, 0, p);
      gsl_vector_view x2 = gsl_vector_subvector(x, p, n - p);
      double rcond;
      size_t i;

      fprintf(stderr, "invert_calc_nonlinear_multifit: computing RHS of linear system (%zu-by-1)...", n);
      gettimeofday(&tv0, NULL);
      invert_calc_Wf(w->c, w, f);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, || f || = %.12e)\n", time_diff(tv0, tv1), gsl_blas_dnrm2(f));

      fprintf(stderr, "invert_calc_nonlinear_multifit: computing design matrix of linear system (%zu-by-%zu)...", n, p);
      gettimeofday(&tv0, NULL);
      invert_calc_df(w->c, w, J);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      fprintf(stderr, "invert_calc_nonlinear_multifit: adding weights to Jacobian matrix...");
      gettimeofday(&tv0, NULL);
      for (i = 0; i < w->nres; ++i)
        {
          double swi = gsl_vector_get(w->sqrt_wts_final, i);
          gsl_vector_view v = gsl_matrix_row(J, i);

          if (swi != 1.0)
            gsl_vector_scale(&v.vector, swi);
        }
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      fprintf(stderr, "invert_calc_nonlinear_multifit: computing QR decomposition of design matrix...");
      gettimeofday(&tv0, NULL);
      gsl_linalg_QR_decomp_r(J, T);
      gettimeofday(&tv1, NULL);
      gsl_linalg_QR_rcond(J, &rcond, work3);
      fprintf(stderr, "done (%g seconds, condition number = %.12e)\n", time_diff(tv0, tv1), 1.0 / rcond);

      fprintf(stderr, "invert_calc_nonlinear_multifit: solving least squares system...");
      gettimeofday(&tv0, NULL);
      gsl_linalg_QR_lssolve_r(J, T, f, x, work);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, residual norm = %.12e)\n",
              time_diff(tv0, tv1), gsl_blas_dnrm2(&x2.vector));

      gsl_vector_memcpy(w->c, &x1.vector);

      printv_octave(w->c, "c");

      gsl_matrix_free(T);
      gsl_vector_free(x);
      gsl_vector_free(work);
      gsl_vector_free(work3);
    }
  else
    {
#if 0
      const size_t max_iter = 50;     /* maximum iterations */
      const double xtol = 1.0e-5;
      const double gtol = 1.0e-6;
      const double ftol = 1.0e-6;
      gsl_vector *f;
      int info;
      gsl_multifit_nlinear_fdf fdf;
      double res0;                    /* initial residual */

      fdf.f = invert_calc_f;
      fdf.df = invert_calc_df;
      fdf.fvv = NULL;
      fdf.n = n;
      fdf.p = p;
      fdf.params = w;

      fprintf(stderr, "invert_calc_nonlinear_multifit: initializing multifit...");
      gettimeofday(&tv0, NULL);
      gsl_multifit_nlinear_winit(c, w->wts_final, &fdf, w->multifit_nlinear_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      /* compute initial residual */
      f = gsl_multifit_nlinear_residual(w->multifit_nlinear_p);
      res0 = gsl_blas_dnrm2(f);

      fprintf(stderr, "invert_calc_nonlinear_multifit: computing nonlinear least squares solution...");
      gettimeofday(&tv0, NULL);
      s = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                                      invert_nonlinear_callback, (void *) w,
                                      &info, w->multifit_nlinear_p);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      if (s == GSL_SUCCESS)
        {
          fprintf(stderr, "invert_calc_nonlinear_multifit: NITER = %zu\n",
                  gsl_multifit_nlinear_niter(w->multifit_nlinear_p));
          fprintf(stderr, "invert_calc_nonlinear_multifit: NFEV  = %zu\n", fdf.nevalf);
          fprintf(stderr, "invert_calc_nonlinear_multifit: NJEV  = %zu\n", fdf.nevaldf);
          fprintf(stderr, "invert_calc_nonlinear_multifit: NAEV  = %zu\n", fdf.nevalfvv);
          fprintf(stderr, "invert_calc_nonlinear_multifit: reason for stopping: %d\n", info);
          fprintf(stderr, "invert_calc_nonlinear_multifit: initial |f(x)|: %.12e\n", res0);
          fprintf(stderr, "invert_calc_nonlinear_multifit: final   |f(x)|: %.12e\n",
                  gsl_blas_dnrm2(f));
        }
      else
        {
          fprintf(stderr, "invert_calc_nonlinear_multifit: multifit failed: %s\n",
                  gsl_strerror(s));
        }

      /* store final coefficients in physical units */
      {
        gsl_vector *x_final = gsl_multifit_nlinear_position(w->multifit_nlinear_p);
        gsl_vector_memcpy(w->c, x_final);
      }
#endif
    }

  return s;
}

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
  double xnorm = gsl_blas_dnrm2(x);
  int xzero = xnorm == 0.0 ? 1 : 0;

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
          double B_model[3];            /* model vector */
          double B_obs[3];              /* observation vector NEC frame */
          double B_prior[3];            /* a priori model NEC frame */
          double dBdt_model[3];         /* dB/dt (SV of internal model) */
          double dBdt_obs[3];           /* SV observation vector (NEC frame) */
          double F_obs;                 /* scalar field measurement */

          /* compute vector model for this residual */
          if (xzero)
            B_model[0] = B_model[1] = B_model[2] = 0.0;
          else
            invert_nonlinear_model(x, mptr->t[j], mptr->r[j], mptr->theta[j], mptr->phi[j], thread_id, B_model, w);

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

          B_prior[0] = mptr->Bx_model[j];
          B_prior[1] = mptr->By_model[j];
          B_prior[2] = mptr->Bz_model[j];

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              gsl_vector_set(f, ridx++, B_obs[0] - B_prior[0] - B_model[0]);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              gsl_vector_set(f, ridx++, B_obs[1] - B_prior[1] - B_model[1]);
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              gsl_vector_set(f, ridx++, B_obs[2] - B_prior[2] - B_model[2]);
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

#if 0/*XXX*/
  if (mparams->regularize && !mparams->synth_data)
    {
      /* store L^T*x in bottom of f for regularization */
      gsl_vector_view v = gsl_vector_subvector(f, w->nres, w->p);
      gsl_spblas_dusmv(CblasTrans, 1.0, w->L, x, 0.0, &v.vector);
    }
#endif

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
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          size_t ridx = mptr->index[j]; /* residual index for this data point */
          gsl_vector_view vx, vy, vz;
          gsl_vector *VX = NULL, *VY = NULL, *VZ = NULL;

          if (mptr->flags[j] & MAGDATA_FLG_X)
            {
              vx = gsl_matrix_row(J, ridx++);
              VX = &vx.vector;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            {
              vy = gsl_matrix_row(J, ridx++);
              VY = &vy.vector;
            }

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            {
              vz = gsl_matrix_row(J, ridx++);
              VZ = &vz.vector;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z))
            {
              invert_jacobian_vector(thread_id, mptr->t[j], r, theta, phi, VX, VY, VZ, w);
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
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DYDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p);
                }

              if (mptr->flags[j] & MAGDATA_FLG_DZDT)
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx++, 0, w->p);
                }
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p);
                }

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p);
                }

              ++ridx;
            }

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            {
              /* check if fitting MF to this data point */
              if (MAGDATA_FitMF(mptr->flags[j]))
                {
                  gsl_vector_view Jv = gsl_matrix_subrow(J, ridx, 0, w->p);
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
  printv_octave(x, "x");
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

        gsl_matrix_free(JTJ);
        exit(1);
      }
  }
#endif

  return s;
}

/*
invert_nonlinear_model()
  Compute total B model vector for a given residual

Inputs: x         - parameter vector
        t         - timestamp (CDF_EPOCH)
        r         - radius (km)
        theta     - colatitude (radians)
        phi       - longitude (radians)
        thread_id - OpenMP thread id
        B_model   - (output) B_model (X,Y,Z) in NEC (nT)
        w         - workspace
*/

static int
invert_nonlinear_model(const gsl_vector * x, const double t, const double r, const double theta, const double phi,
                       const int thread_id, double B_model[3], invert_workspace *w)
{
  int s = 0;
  gsl_matrix * JB = w->omp_B[thread_id];
  gsl_vector_view X = gsl_matrix_row(w->omp_B[thread_id], 0);
  gsl_vector_view Y = gsl_matrix_row(w->omp_B[thread_id], 1);
  gsl_vector_view Z = gsl_matrix_row(w->omp_B[thread_id], 2);
  gsl_vector_view out = gsl_vector_view_array(B_model, 3);

  /* compute 3 rows of Jacobian corresponding to X,Y,Z residuals */
  s = invert_jacobian_vector(thread_id, t, r, theta, phi, &X.vector, &Y.vector, &Z.vector, w);
  if (s)
    return s;

  /* B_model = JB * x */
  gsl_blas_dgemv(CblasNoTrans, 1.0, JB, x, 0.0, &out.vector);

  return s;
}

/*
invert_jacobian_vector()
  Construct a block of 3 rows of the Jacobian corresponding to vector residuals

Inputs: thread_id - thread id
        t         - timestamp (CDF_EPOCH)
        r         - radius (km)
        theta     - colatitude (radians)
        phi       - longitude (radians)
        X         - (output) X row of Jacobian, length p
        Y         - (output) X row of Jacobian, length p
        Z         - (output) X row of Jacobian, length p
        w         - workspace
*/

static int
invert_jacobian_vector(const int thread_id, const double t, const double r, const double theta, const double phi,
                       gsl_vector * X, gsl_vector * Y, gsl_vector * Z,
                       const invert_workspace *w)
{
  int s = 0;
  invert_tmode_workspace * tmode_p = w->tmode_workspace_p;
  invert_smode_workspace * smode_p = w->smode_workspace_p;
  const size_t nfreq = w->nfreq;
  size_t i, j, k;

  /* precompute magfield arrays for this (r,theta,phi) */
  invert_smode_precompute(thread_id, r, theta, phi, smode_p);

#if 0
  for (i = 0; i < nfreq; ++i)
    {
      for (j = 0; j < smode_p->nmodes[i]; ++j)
        {
          gsl_complex phi_ij[3] = { GSL_COMPLEX_ONE, GSL_COMPLEX_ONE, GSL_COMPLEX_ONE };/*XXX*/

          /*invert_smode_get(r, theta, phi, i, j, phi_ij, smode_p);*/

          for (k = 0; k < tmode_p->nmodes[i]; ++k)
            {
              /* compute temporal mode k of band i for this timestamp */
              complex double alpha_ik = CastComplex(invert_tmode_get(t, i, k, tmode_p));
              size_t idx = invert_coeff_idx(i, k, j, w);
              complex double z;

              if (X)
                {
                  z = alpha_ik * CastComplex(phi_ij[0]);
                  gsl_vector_set(X, idx, creal(z));
                  gsl_vector_set(X, w->p_complex + idx, -cimag(z));
                }

              if (Y)
                {
                  z = alpha_ik * CastComplex(phi_ij[1]);
                  gsl_vector_set(Y, idx, creal(z));
                  gsl_vector_set(Y, w->p_complex + idx, -cimag(z));
                }

              if (Z)
                {
                  z = alpha_ik * CastComplex(phi_ij[2]);
                  gsl_vector_set(Z, idx, creal(z));
                  gsl_vector_set(Z, w->p_complex + idx, -cimag(z));
                }
            }
        }
    }
#elif 1

  for (i = 0; i < nfreq; ++i)
    {
      for (j = 0; j < smode_p->nmodes[i]; ++j)
        {
          gsl_complex phi_ij[3];

          invert_smode_get(thread_id, r, theta, phi, i, j, phi_ij, smode_p);

          for (k = 0; k < tmode_p->nmodes[i]; ++k)
            {
              /* compute temporal mode k of band i for this timestamp */
              gsl_complex alpha_ik = invert_tmode_get(t, i, k, tmode_p);
              size_t idx = invert_coeff_idx(i, k, j, w);

              if (X)
                {
                  gsl_complex z = gsl_complex_mul(alpha_ik, phi_ij[0]);
                  gsl_vector_set(X, idx, GSL_REAL(z));
                  gsl_vector_set(X, w->p_complex + idx, -GSL_IMAG(z));
                }

              if (Y)
                {
                  gsl_complex z = gsl_complex_mul(alpha_ik, phi_ij[1]);
                  gsl_vector_set(Y, idx, GSL_REAL(z));
                  gsl_vector_set(Y, w->p_complex + idx, -GSL_IMAG(z));
                }

              if (Z)
                {
                  gsl_complex z = gsl_complex_mul(alpha_ik, phi_ij[2]);
                  gsl_vector_set(Z, idx, GSL_REAL(z));
                  gsl_vector_set(Z, w->p_complex + idx, -GSL_IMAG(z));
                }
            }
        }
    }
#endif

  return s;
}

/*
invert_nonlinear_model_J()
  Compute total J model vector for a given residual

Inputs: x         - parameter vector
        t         - timestamp (CDF_EPOCH)
        r         - radius (km)
        theta     - colatitude (radians)
        phi       - longitude (radians)
        thread_id - OpenMP thread id
        J_model   - (output) J_model (X,Y,Z) in NEC
        w         - workspace
*/

static int
invert_nonlinear_model_J(const gsl_vector * x, const double t, const double r,
                         const double theta, const double phi,
                         const int thread_id, double J_model[3], invert_workspace *w)
{
  int s = 0;
  gsl_matrix * JB = w->omp_B[thread_id];
  gsl_vector_view X = gsl_matrix_row(w->omp_B[thread_id], 0);
  gsl_vector_view Y = gsl_matrix_row(w->omp_B[thread_id], 1);
  gsl_vector_view Z = gsl_matrix_row(w->omp_B[thread_id], 2);
  gsl_vector_view out = gsl_vector_view_array(J_model, 3);

  /* compute 3 rows of Jacobian corresponding to X,Y,Z residuals */
  s = invert_jacobian_vector_J(thread_id, t, r, theta, phi, &X.vector, &Y.vector, &Z.vector, w);
  if (s)
    return s;

  /* J_model = JB * x */
  gsl_blas_dgemv(CblasNoTrans, 1.0, JB, x, 0.0, &out.vector);

  return s;
}

/*
invert_jacobian_vector_J()
  Construct a block of 3 rows of the Jacobian corresponding to vector residuals
for the current density

Inputs: thread_id - thread id
        t         - timestamp (CDF_EPOCH)
        r         - radius (km)
        theta     - colatitude (radians)
        phi       - longitude (radians)
        X         - (output) X row of Jacobian, length p
        Y         - (output) X row of Jacobian, length p
        Z         - (output) X row of Jacobian, length p
        w         - workspace
*/

static int
invert_jacobian_vector_J(const int thread_id, const double t, const double r, const double theta, const double phi,
                         gsl_vector * X, gsl_vector * Y, gsl_vector * Z,
                         const invert_workspace *w)
{
  int s = 0;
  invert_tmode_workspace * tmode_p = w->tmode_workspace_p;
  invert_smode_workspace * smode_p = w->smode_workspace_p;
  const size_t nfreq = w->nfreq;
  size_t i, j, k;

  /* precompute magfield arrays for this (r,theta,phi) */
  invert_smode_precompute_J(thread_id, r, theta, phi, smode_p);

  for (i = 0; i < nfreq; ++i)
    {
      for (j = 0; j < smode_p->nmodes[i]; ++j)
        {
          gsl_complex phi_ij[3];

          invert_smode_get_J(thread_id, r, theta, phi, i, j, phi_ij, smode_p);

          for (k = 0; k < tmode_p->nmodes[i]; ++k)
            {
              /* compute temporal mode j of band i for this timestamp */
              gsl_complex alpha_ik = invert_tmode_get(t, i, k, tmode_p);
              size_t idx = invert_coeff_idx(i, k, j, w);

              if (X)
                {
                  gsl_complex z = gsl_complex_mul(alpha_ik, phi_ij[0]);
                  gsl_vector_set(X, idx, GSL_REAL(z));
                  gsl_vector_set(X, w->p_complex + idx, -GSL_IMAG(z));
                }

              if (Y)
                {
                  gsl_complex z = gsl_complex_mul(alpha_ik, phi_ij[1]);
                  gsl_vector_set(Y, idx, GSL_REAL(z));
                  gsl_vector_set(Y, w->p_complex + idx, -GSL_IMAG(z));
                }

              if (Z)
                {
                  gsl_complex z = gsl_complex_mul(alpha_ik, phi_ij[2]);
                  gsl_vector_set(Z, idx, GSL_REAL(z));
                  gsl_vector_set(Z, w->p_complex + idx, -GSL_IMAG(z));
                }
            }
        }
    }

  return s;
}
