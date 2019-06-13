/*
 * invert.c
 *
 * This module contains code to perform magnetic main field modeling
 *
 * Calling sequences:
 * 1. invert_data_alloc     - allocate invert_data_workspace structure
 * 2. invert_data_copy      - copy all satellite data to internal
 *                            structure, ignoring flagged data
 * 3. invert_data_map       - print out map of spatial coverage of data
 * 4. invert_alloc          - allocate invert_workspace
 * 5. invert_init           - initialize various internal parameters
 *                            (time scaling, etc)
 * 6. invert_calc_linear    - build linear LS matrix and solve system
 * 7. invert_calc_nonlinear - solve nonlinear LS system; the initial
 *                            guess for the nonlinear procedure is the
 *                            linearized solution
 * 8. invert_free
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <stdarg.h>

#include <omp.h>

#include <satdata/satdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_spmatrix.h>

#include <spblas2/gsl_spblas2.h>
#include <spblas2/gsl_spcblas.h>

#include <common/common.h>
#include <common/oct.h>
#include <track/track_weight.h>
#include <bspline2/gsl_bspline2.h>

#include "lls.h"
#include "lapack_wrapper.h"
#include "invert.h"
#include "invert_tmode.h"

static int invert_compare_int(const void *a, const void *b);
static int invert_debug(const char *format, ...);

#include "invert_nonlinear.c"

/*
invert_alloc()
  Allocate a invert workspace

Inputs: params - model parameters
                 epoch   - model epoch t0 (years)
                 R       - reference radius (km)
                 nsat    - number of different satellites
                 flags   - flags
*/

invert_workspace *
invert_alloc(const invert_parameters *params)
{
  invert_workspace *w;
  const size_t ntheta = 100;
  const size_t nphi = 100;
  size_t plm_size;
  size_t i, j;

  w = calloc(1, sizeof(invert_workspace));
  if (!w)
    return 0;

  w->nsat = params->nsat;
  w->epoch = params->epoch;
  w->R = params->R;
  w->nmax = params->nmax;
  w->data_workspace_p = params->invert_data_p;
  w->max_threads = (size_t) omp_get_max_threads();

  w->params = *params;

  w->weight_workspace_p = track_weight_alloc(ntheta, nphi);
  w->spatwtMF_workspace_p = spatwt_alloc(8, 12);
  w->spatwtSV_workspace_p = spatwt_alloc(8, 12);

  /*
   * Add up all the contributions to the coefficient vector, which will
   * be partitioned as:
   *
   * c = [ MF | SV | SA | Euler | external ]
   */

  /* subtract 1 to exclude the (0,0) coefficient */
  w->nnm_tot = (w->nmax + 1) * (w->nmax + 1) - 1;

  /*XXX*/
  w->nnm_max = w->nnm_tot;

  /* compute total (internal) model coefficients */

  w->p_core = 0;
  if (params->fit_mf && w->data_workspace_p)
    {
    }

  w->p_int = GSL_MAX(w->p_core, 1);

  if (w->p_int == 0)
    {
      invert_free(w);
      GSL_ERROR_NULL("no internal model parameters to fit", GSL_EINVAL);
    }

  w->p = w->p_int;

  /* compute max nmax */
  w->nmax_max = 0;

  plm_size = gsl_sf_legendre_array_n(w->nmax);

  if (w->p == 0)
    {
      invert_free(w);
      GSL_ERROR_NULL("no parameters to fit", GSL_EINVAL);
    }

  w->cosmphi = malloc((w->nmax + 1) * sizeof(double));
  w->sinmphi = malloc((w->nmax + 1) * sizeof(double));

  w->Plm = malloc(plm_size * sizeof(double));
  w->dPlm = malloc(plm_size * sizeof(double));
  if (!w->Plm || !w->dPlm)
    {
      invert_free(w);
      return 0;
    }

  w->c = gsl_vector_calloc(w->p);
  w->c_copy = gsl_vector_alloc(w->p);

  /* covariance matrix of the internal field coefficients */
  w->covar = gsl_matrix_alloc(w->p, w->p);

  w->dX = malloc(w->nnm_tot * sizeof(double));
  w->dY = malloc(w->nnm_tot * sizeof(double));
  w->dZ = malloc(w->nnm_tot * sizeof(double));

  w->nobs_cnt = 0;

  /* these are computed later in invert_init() */
  w->t0_data = 0.0;

  w->niter = 0;

  w->JTJ_vec = gsl_matrix_alloc(w->p_int, w->p_int);
  w->choleskyL = gsl_matrix_alloc(w->p_int, w->p_int);

  w->eigen_workspace_p = gsl_eigen_symm_alloc(w->p);

  w->L = gsl_spmatrix_alloc(w->p, w->p);
  w->Lambda = gsl_spmatrix_alloc(w->p, w->p);

  w->omp_dX = gsl_matrix_alloc(w->max_threads, w->nnm_tot);
  w->omp_dY = gsl_matrix_alloc(w->max_threads, w->nnm_tot);
  w->omp_dZ = gsl_matrix_alloc(w->max_threads, w->nnm_tot);
  w->omp_dF = gsl_matrix_alloc(w->max_threads, w->nnm_tot);
  w->omp_dX_grad = gsl_matrix_alloc(w->max_threads, w->nnm_tot);
  w->omp_dY_grad = gsl_matrix_alloc(w->max_threads, w->nnm_tot);
  w->omp_dZ_grad = gsl_matrix_alloc(w->max_threads, w->nnm_tot);

  /* initialize green workspaces and omp_J */
  {
    /*
     * Number of components added to omp_J matrix:
     * X, Y, Z, F, DX_NS, DY_NS, DZ_NS, DF_NS, DX_EW, DY_EW, DZ_EW, DF_EW, dX/dt, dY/dt, dZ/dt
     */
    /*const size_t ncomp = 15;*/
    const size_t ncomp = 4;

    /*
     * maximum observations to accumulate at once in LS system, calculated to make
     * each omp_J matrix approximately of size 'MFIELD_MATRIX_SIZE'
     */
    w->data_block = MFIELD_MATRIX_SIZE / (ncomp * w->p_int * sizeof(double));
    w->data_block_tot = MFIELD_MATRIX_SIZE / (ncomp * w->nnm_tot * sizeof(double));

    w->green_array_p = malloc(w->max_threads * sizeof(green_workspace *));
    w->omp_J = malloc(w->max_threads * sizeof(gsl_matrix *));
    w->omp_dB = malloc(w->max_threads * sizeof(gsl_matrix *));
    w->omp_rowidx = malloc(w->max_threads * sizeof(size_t));
    w->omp_T = malloc(w->max_threads * sizeof(gsl_matrix *));
    w->omp_S = malloc(w->max_threads * sizeof(gsl_matrix *));
    w->omp_colidx = malloc(w->max_threads * sizeof(size_t));

    for (i = 0; i < w->max_threads; ++i)
      {
        w->green_array_p[i] = green_alloc(w->nmax, w->nmax, w->R);
        w->omp_J[i] = gsl_matrix_alloc(ncomp * w->data_block, w->p_int);
        w->omp_dB[i] = gsl_matrix_alloc(3, w->nnm_tot);
        w->omp_S[i] = gsl_matrix_alloc(w->nnm_tot, ncomp * w->data_block_tot);
        w->omp_T[i] = gsl_matrix_alloc(ncomp * w->data_block_tot, w->nnm_tot);
      }

    fprintf(stderr, "invert_alloc: data_block     = %zu\n", w->data_block);
    fprintf(stderr, "invert_alloc: data_block_tot = %zu\n", w->data_block_tot);
  }

  w->tmode_workspace_p = invert_tmode_read_binary(params->tmode_file);

  return w;
}

void
invert_free(invert_workspace *w)
{
  size_t i;

  if (w->cosmphi)
    free(w->cosmphi);

  if (w->sinmphi)
    free(w->sinmphi);

  if (w->Plm)
    free(w->Plm);

  if (w->dPlm)
    free(w->dPlm);

  if (w->c)
    gsl_vector_free(w->c);

  if (w->c_copy)
    gsl_vector_free(w->c_copy);

  if (w->covar)
    gsl_matrix_free(w->covar);

  if (w->dX)
    free(w->dX);

  if (w->dY)
    free(w->dY);

  if (w->dZ)
    free(w->dZ);

  if (w->bias_idx)
    free(w->bias_idx);

  if (w->L)
    gsl_spmatrix_free(w->L);

  if (w->Lambda)
    gsl_spmatrix_free(w->Lambda);

  if (w->wts_spatial)
    gsl_vector_free(w->wts_spatial);

  if (w->wts_robust)
    gsl_vector_free(w->wts_robust);

  if (w->wts_final)
    gsl_vector_free(w->wts_final);

  if (w->sqrt_wts_final)
    gsl_vector_free(w->sqrt_wts_final);

  if (w->multifit_nlinear_p)
    gsl_multifit_nlinear_free(w->multifit_nlinear_p);

  if (w->fvec)
    gsl_vector_free(w->fvec);

  if (w->wfvec)
    gsl_vector_free(w->wfvec);

  if (w->robust_workspace_p)
    gsl_multifit_robust_free(w->robust_workspace_p);

  if (w->weight_workspace_p)
    track_weight_free(w->weight_workspace_p);

  if (w->spatwtMF_workspace_p)
    spatwt_free(w->spatwtMF_workspace_p);

  if (w->spatwtSV_workspace_p)
    spatwt_free(w->spatwtSV_workspace_p);

  if (w->JTJ_vec)
    gsl_matrix_free(w->JTJ_vec);

  if (w->choleskyL)
    gsl_matrix_free(w->choleskyL);

  if (w->omp_dX)
    gsl_matrix_free(w->omp_dX);

  if (w->omp_dY)
    gsl_matrix_free(w->omp_dY);

  if (w->omp_dZ)
    gsl_matrix_free(w->omp_dZ);

  if (w->omp_dF)
    gsl_matrix_free(w->omp_dF);

  if (w->omp_dX_grad)
    gsl_matrix_free(w->omp_dX_grad);

  if (w->omp_dY_grad)
    gsl_matrix_free(w->omp_dY_grad);

  if (w->omp_dZ_grad)
    gsl_matrix_free(w->omp_dZ_grad);

  if (w->nlinear_workspace_p)
    gsl_multilarge_nlinear_free(w->nlinear_workspace_p);

  if (w->eigen_workspace_p)
    gsl_eigen_symm_free(w->eigen_workspace_p);

  for (i = 0; i < w->max_threads; ++i)
    {
      green_free(w->green_array_p[i]);
      gsl_matrix_free(w->omp_J[i]);
      gsl_matrix_free(w->omp_dB[i]);
      gsl_matrix_free(w->omp_T[i]);
      gsl_matrix_free(w->omp_S[i]);
    }

  free(w->green_array_p);
  free(w->omp_J);
  free(w->omp_dB);
  free(w->omp_T);
  free(w->omp_S);
  free(w->omp_rowidx);
  free(w->omp_colidx);

  if (w->tmode_workspace_p)
    invert_tmode_free(w->tmode_workspace_p);

  free(w);
}

/* initialize parameter structure */
int
invert_init_params(invert_parameters * params)
{
  params->epoch = -1.0;
  params->R = -1.0;
  params->nmax = 0;
  params->nsat = 0;
  params->max_iter = 0;
  params->fit_mf = 0;
  params->qdlat_fit_cutoff = -1.0;
  params->regularize = 0;
  params->use_weights = 0;
  params->weight_X = 0.0;
  params->weight_Y = 0.0;
  params->weight_Z = 0.0;
  params->weight_F = 0.0;
  params->weight_DXDT = 0.0;
  params->weight_DYDT = 0.0;
  params->weight_DZDT = 0.0;
  params->weight_DX = 0.0;
  params->weight_DY = 0.0;
  params->weight_DZ = 0.0;
  params->synth_data = 0;
  params->synth_noise = 0;
  params->synth_nmin = 0;

  return 0;
}

/*
invert_copy()
  Make a copy of an invert workspace that can be used
interchangeably with invert_eval(). This is useful for
OpenMP applications, since invert_eval() uses internal
arrays (dX,dY,dZ) so we need separate workspaces for each
thread for evaluation
*/

invert_workspace *
invert_copy(const invert_workspace *w)
{
  invert_workspace *w_copy;

  w_copy = invert_alloc(&(w->params));
  if (!w_copy)
    return 0;

  /* copy coefficients */
  gsl_vector_memcpy(w_copy->c, w->c);

  return w_copy;
} /* invert_copy() */

/*
invert_reset()
  Reset invert workspace to be ready for a new
LS iteration
*/

int
invert_reset(invert_workspace *w)
{
  int s = 0;

#if 0
  /* reset weight histogram */
  s += weight_reset(w->weight_workspace_p);
#endif

  w->nobs_cnt = 0;

  return s;
}

/*
invert_init()
  Initialize model

Notes:
1) On input, w->data_workspace_p must be initialized with satellite data

2) On output, w->t0_data is initialized to the time of the first available
data point (CDF_EPOCH)

3) Area/density weights are computed from previously constructed
histogram
*/

int
invert_init(invert_workspace *w)
{
  int s = 0;
  const invert_parameters *params = &(w->params);
  size_t i, j;

  invert_data_init(w->data_workspace_p);

  fprintf(stderr, "invert_init: data tmin = %.2f\n", satdata_epoch2year(w->data_workspace_p->t0_data));
  fprintf(stderr, "invert_init: data tmax = %.2f\n", satdata_epoch2year(w->data_workspace_p->t1_data));

  /* find time of first available data in CDF_EPOCH */
  w->t0_data = w->data_workspace_p->t0_data;

  /* initialize spatial weighting histogram and time scaling */
  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = invert_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          const double u = satdata_epoch2year(mptr->t[j]) - w->epoch;
          const double v = satdata_epoch2year(mptr->t_ns[j]) - w->epoch;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY)
            {
              if (mptr->flags[j] & MAGDATA_FLG_X)
                spatwt_add_data(mptr->theta[j], mptr->phi[j], w->spatwtMF_workspace_p);

              if (mptr->flags[j] & MAGDATA_FLG_Y)
                spatwt_add_data(mptr->theta[j], mptr->phi[j], w->spatwtMF_workspace_p);

              if (mptr->flags[j] & MAGDATA_FLG_Z)
                spatwt_add_data(mptr->theta[j], mptr->phi[j], w->spatwtMF_workspace_p);

              if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
                spatwt_add_data(mptr->theta[j], mptr->phi[j], w->spatwtMF_workspace_p);
            }
          else if (mptr->global_flags & MAGDATA_GLOBFLG_OBSERVATORY_SV)
            {
              if (mptr->flags[j] & MAGDATA_FLG_DXDT)
                spatwt_add_data(mptr->theta[j], mptr->phi[j], w->spatwtSV_workspace_p);

              if (mptr->flags[j] & MAGDATA_FLG_DYDT)
                spatwt_add_data(mptr->theta[j], mptr->phi[j], w->spatwtSV_workspace_p);

              if (mptr->flags[j] & MAGDATA_FLG_DZDT)
                spatwt_add_data(mptr->theta[j], mptr->phi[j], w->spatwtSV_workspace_p);
            }

#if 1
          if (mptr->flags[j] & MAGDATA_FLG_X)
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);

          if (mptr->flags[j] & MAGDATA_FLG_Y)
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);

          if (mptr->flags[j] & MAGDATA_FLG_Z)
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);

          if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);

          if (MAGDATA_ExistDF_NS(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);

          if (MAGDATA_ExistDF_EW(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);
#else
          track_weight_add_data(mptr->theta[j], mptr->phi[j], w->weight_workspace_p);
#endif
        }
    }

  /* compute data weights with histogram */
  track_weight_calc(w->weight_workspace_p);
  spatwt_calc(w->spatwtMF_workspace_p);
  spatwt_calc(w->spatwtSV_workspace_p);

  spatwt_print("spatwtMF.txt", w->spatwtMF_workspace_p);
  spatwt_print("spatwtSV.txt", w->spatwtSV_workspace_p);

  invert_init_nonlinear(w);

  return s;
} /* invert_init() */

/*
invert_eval()
  Evaluate magnetic field model at given point

Inputs: t     - timestamp (CDF_EPOCH)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - geocentric longitude (radians)
        B     - (output) magnetic field
                B[0] = B_x
                B[1] = B_y
                B[2] = B_z
                B[3] = |B|
        w     - workspace
*/

int
invert_eval(const double t, const double r, const double theta,
            const double phi, double B[4], invert_workspace *w)
{
#if 0/*XXX*/
  int s = invert_eval_g(t, r, theta, phi, w->c, B, w);
  return s;
#endif
}

int
invert_write(const char *filename, invert_workspace *w)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "invert_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fwrite(&(w->params), sizeof(invert_parameters), 1, fp);
  fwrite(&(w->t0_data), sizeof(double), 1, fp);

  /*
   * only write internal coefficients since when we later read
   * the file we won't be able to recalculate w->p_euler
   */
  fwrite(w->c->data, sizeof(double), w->p_int, fp);

  gsl_matrix_fwrite(fp, w->covar);

  fclose(fp);

  return s;
} /* invert_write() */

invert_workspace *
invert_read(const char *filename)
{
  invert_workspace *w;
  invert_parameters params;
  FILE *fp;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "invert_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  fread(&params, sizeof(invert_parameters), 1, fp);

  params.invert_data_p = NULL;

  w = invert_alloc(&params);

  fread(&(w->t0_data), sizeof(double), 1, fp);
  fread(w->c->data, sizeof(double), w->p_int, fp);
  gsl_matrix_fread(fp, w->covar);

  fclose(fp);

  return w;
} /* invert_read() */

/*
invert_write_ascii()
  Write ascii coefficient file

Inputs: filename - output file
        epoch    - epoch to evaluate B-splines (decimal years)
        c        - coefficient vector, length at least p_int
        w        - workspace
*/

int
invert_write_ascii(const char *filename, const double epoch,
                   const gsl_vector * c, invert_workspace *w)
{
  int s = 0;
  const invert_parameters *params = &(w->params);
  FILE *fp;
  size_t n;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "invert_write_ascii: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  /* print header information */
  fprintf(fp, "%% Magnetic field model coefficients\n");
  fprintf(fp, "%% nmax:  %zu\n", w->nmax);
  fprintf(fp, "%% epoch: %.4f\n", epoch);
  fprintf(fp, "%% radius: %.1f\n", w->R);

  fclose(fp);

  return s;
}

/*******************************************************
 *      INTERNAL ROUTINES                              *
 *******************************************************/

/*
invert_coeff_nmidx()
  This function returns a unique index in [0,w->p-1] corresponding
to a given (l,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Inputs: n - SH degree (> 0)
        m - SH order (-l <= m <= l)

Return: index in [0,nnm-1]
*/

size_t
invert_coeff_nmidx(const size_t n, const int m)
{
  size_t base = n * n; /* index of block for this n */
  int offset = m + n;  /* offset within block for this m */
  size_t nmidx;

  if (n == 0)
    {
      fprintf(stderr, "invert_coeff_nmidx: error: n = 0\n");
      return 0;
    }

  nmidx = base + offset;

  /* subtract 1 to exclude (0,0) coefficient */
  return nmidx - 1;
} /* invert_coeff_nmidx() */

static int
invert_compare_int(const void *a, const void *b)
{
  int ai = *(int *) a;
  int bi = *(int *) b;

  if (ai < bi)
    return -1;
  else if (ai == bi)
    return 0;
  else
    return 1;
}

static int
invert_debug(const char *format, ...)
{
  int s = 0;
#if MFIELD_DEBUG
  va_list args;

  va_start(args, format);

  vfprintf(stderr, format, args);

  va_end(args);

  fflush(stderr);
#endif

  return s;
}
