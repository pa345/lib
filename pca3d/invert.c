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

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_test.h>

#include <spblas2/gsl_spblas2.h>
#include <spblas2/gsl_spcblas.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_track_weight.h>
#include <mainlib/ml_lapack_wrapper.h>
#include <bspline2/gsl_bspline2.h>

#include "invert.h"
#include "invert_tmode.h"

static int invert_compare_int(const void *a, const void *b);

#include "invert_nonlinear.c"

/*
invert_alloc()
  Allocate a invert workspace

Inputs: params - model parameters
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
  size_t * nsmodes;
  size_t i;

  w = calloc(1, sizeof(invert_workspace));
  if (!w)
    return 0;

  w->nsat = params->nsat;
  w->R = params->R;
  w->nfreq = params->nfreq;
  w->data_workspace_p = params->invert_data_p;
  w->max_threads = (size_t) omp_get_max_threads();

  w->params = *params;

  w->weight_workspace_p = track_weight_alloc(ntheta, nphi);
  w->spatwtMF_workspace_p = spatwt_alloc(8, 12);
  w->spatwtSV_workspace_p = spatwt_alloc(8, 12);

  fprintf(stderr, "invert_alloc: reading temporal modes from %s...", params->tmode_file);
  w->tmode_workspace_p = invert_tmode_read_binary(params->tmode_file);
  w->tmode_workspace_p->nfreq = w->nfreq;
  fprintf(stderr, "done (%zu frequency bands)\n", w->tmode_workspace_p->nfreq);

  /* use same number of spatial modes per bin as the temporal modes */
  nsmodes = malloc(w->nfreq * sizeof(size_t));
#if 0 /*XXX*/
  for (i = 0; i < w->nfreq; ++i)
    nsmodes[i] = w->tmode_workspace_p->nmodes[i];
#elif 1
  for (i = 0; i < w->nfreq; ++i)
    nsmodes[i] = 15;

#if 0
  nsmodes[10] = 20; /* 1 cpd */
#endif
#else
  for (i = 0; i < w->nfreq; ++i)
    nsmodes[i] = GSL_MIN(w->tmode_workspace_p->nmodes[i], 6);
#endif

  fprintf(stderr, "invert_alloc: reading spatial modes...");
  w->smode_workspace_p = invert_smode_alloc(w->nfreq, nsmodes);
  w->smode_workspace_p->nfreq = w->nfreq;
  fprintf(stderr, "done (%zu frequency bands)\n", w->smode_workspace_p->nfreq);

  /* sanity check that tmodes and smodes have same frequency bands */
  for (i = 0; i < w->nfreq; ++i)
    {
      double tfreq = w->tmode_workspace_p->freqs[i];
      double sfreq = w->smode_workspace_p->freqs[i];

      gsl_test_rel(tfreq, sfreq, GSL_DBL_EPSILON, "frequency check ifreq=%zu", i);
    }

  w->tmode_idx = malloc(w->tmode_workspace_p->nfreq * sizeof(size_t));
  w->smode_idx = malloc(w->smode_workspace_p->nfreq * sizeof(size_t));
  w->mode_idx = malloc(w->smode_workspace_p->nfreq * sizeof(size_t));

  /*
   * coefficient vector will be partitioned as:
   *
   * c = [ beta_r | beta_i ]
   */

  /* count number of temporal modes in each frequency band */
  w->ntmodes = 0;
  for (i = 0; i < w->tmode_workspace_p->nfreq; ++i)
    {
      w->tmode_idx[i] = w->ntmodes;
      w->ntmodes += w->tmode_workspace_p->nmodes[i];
    }

  /* count number of spatial modes in each frequency band */
  w->nsmodes = 0;
  for (i = 0; i < w->smode_workspace_p->nfreq; ++i)
    {
      w->smode_idx[i] = w->nsmodes;
      w->nsmodes += w->smode_workspace_p->nmodes[i];
    }

  /* compute total number of model coefficients */
  w->p_complex = 0;
  for (i = 0; i < w->nfreq; ++i)
    {
      w->mode_idx[i] = w->p_complex;
      w->p_complex += w->tmode_workspace_p->nmodes[i] * w->smode_workspace_p->nmodes[i];
    }

  /*w->p_complex = w->ntmodes * w->nsmodes;*/
  w->p = 2 * w->p_complex;

  if (w->p == 0)
    {
      invert_free(w);
      GSL_ERROR_NULL("no model parameters to fit", GSL_EINVAL);
    }

  w->c = gsl_vector_calloc(w->p);
  w->c_copy = gsl_vector_alloc(w->p);

  /* covariance matrix of model parameters */
  w->covar = gsl_matrix_alloc(w->p, w->p);

  w->nobs_cnt = 0;

  /* these are computed later in invert_init() */
  w->t0_data = 0.0;

  w->niter = 0;

  w->JTJ_vec = gsl_matrix_alloc(w->p, w->p);
  w->choleskyL = gsl_matrix_alloc(w->p, w->p);

  w->L = gsl_spmatrix_alloc(w->p, w->p);
  w->Lambda = gsl_spmatrix_alloc(w->p, w->p);

  {
    /*
     * Number of components added to omp_J matrix:
     * X, Y, Z, F, DX_NS, DY_NS, DZ_NS, DF_NS, DX_EW, DY_EW, DZ_EW, DF_EW, dX/dt, dY/dt, dZ/dt
     */
    /*const size_t ncomp = 15;*/
    const size_t ncomp = 4;

    /*
     * maximum observations to accumulate at once in LS system, calculated to make
     * each omp_J matrix approximately of size 'INVERT_MATRIX_SIZE'
     */
    w->data_block = INVERT_MATRIX_SIZE / (ncomp * w->p * sizeof(double));

    w->omp_J = malloc(w->max_threads * sizeof(gsl_matrix *));
    w->omp_f = malloc(w->max_threads * sizeof(gsl_vector *));
    w->omp_B = malloc(w->max_threads * sizeof(gsl_matrix *));
    w->omp_rowidx = malloc(w->max_threads * sizeof(size_t));

    for (i = 0; i < w->max_threads; ++i)
      {
        w->omp_J[i] = gsl_matrix_alloc(ncomp * w->data_block, w->p);
        w->omp_f[i] = gsl_vector_alloc(ncomp * w->data_block);
        w->omp_B[i] = gsl_matrix_alloc(3, w->p);
      }

    fprintf(stderr, "invert_alloc: data_block = %zu\n", w->data_block);
    fprintf(stderr, "invert_alloc: nrows      = %zu\n", ncomp * w->data_block);
  }

  return w;
}

void
invert_free(invert_workspace *w)
{
  size_t i;

  if (w->c)
    gsl_vector_free(w->c);

  if (w->c_copy)
    gsl_vector_free(w->c_copy);

  if (w->covar)
    gsl_matrix_free(w->covar);

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

  if (w->nlinear_workspace_p)
    gsl_multilarge_nlinear_free(w->nlinear_workspace_p);

  if (w->multilarge_linear_p)
    gsl_multilarge_linear_free(w->multilarge_linear_p);

  for (i = 0; i < w->max_threads; ++i)
    {
      gsl_matrix_free(w->omp_J[i]);
      gsl_vector_free(w->omp_f[i]);
      gsl_matrix_free(w->omp_B[i]);
    }

  free(w->omp_J);
  free(w->omp_f);
  free(w->omp_B);
  free(w->omp_rowidx);

  if (w->tmode_workspace_p)
    invert_tmode_free(w->tmode_workspace_p);

  if (w->smode_workspace_p)
    invert_smode_free(w->smode_workspace_p);

  if (w->tmode_idx)
    free(w->tmode_idx);

  if (w->smode_idx)
    free(w->smode_idx);

  if (w->mode_idx)
    free(w->mode_idx);

  free(w);
}

/* initialize parameter structure */
int
invert_init_params(invert_parameters * params)
{
  params->R = -1.0;
  params->nsat = 0;
  params->max_iter = 0;
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
          const double u = satdata_epoch2year(mptr->t[j]);
          const double v = satdata_epoch2year(mptr->t_ns[j]);

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

  fwrite(&(w->c->size), sizeof(size_t), 1, fp);
  gsl_vector_fwrite(fp, w->c);

  fclose(fp);

  return s;
}

int
invert_read(const char *filename, gsl_vector *c)
{
  FILE *fp;
  size_t n;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "invert_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  fread(&n, sizeof(size_t), 1, fp);
  if (c->size != n)
    {
      fprintf(stderr, "mfield_read: coefficient vector has wrong size\n");
      return -1;
    }

  gsl_vector_fread(fp, c);

  fclose(fp);

  return 0;
}

int
invert_write_matrix(const char *filename, invert_workspace *w)
{
  int s = 0;
  const gsl_matrix * A = gsl_multilarge_linear_matrix_ptr(w->multilarge_linear_p);
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "invert_write_matrix: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fwrite(&(A->size1), sizeof(size_t), 1, fp);
  fwrite(&(A->size2), sizeof(size_t), 1, fp);
  gsl_matrix_fwrite(fp, A);

  fclose(fp);

  return s;
}

int
invert_write_rhs(const char *filename, invert_workspace *w)
{
  int s = 0;
  const gsl_vector * rhs = gsl_multilarge_linear_rhs_ptr(w->multilarge_linear_p);
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "invert_write_matrix: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fwrite(&(rhs->size), sizeof(size_t), 1, fp);
  gsl_vector_fwrite(fp, rhs);

  fclose(fp);

  return s;
}

/*
invert_write_ascii()
  Write ascii coefficient file

Inputs: filename - output file
        c        - coefficient vector, length at least p
        w        - workspace
*/

int
invert_write_ascii(const char *filename, const gsl_vector * c, invert_workspace *w)
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
  fprintf(fp, "%% radius: %.1f\n", w->R);

  gsl_vector_fprintf(fp, w->c, "%.12e");

  fclose(fp);

  return s;
}

int
invert_debug(const char *format, ...)
{
  int s = 0;
#if INVERT_DEBUG
  va_list args;

  va_start(args, format);

  vfprintf(stderr, format, args);

  va_end(args);

  fflush(stderr);
#endif

  return s;
}

/*******************************************************
 *      INTERNAL ROUTINES                              *
 *******************************************************/

/*
invert_coeff_idx()
  This function returns a unique index in [0,w->p_complex-1] corresponding
to a given set of (f,tmode,smode).

Inputs: f     - frequency band, in [0,nfreq-1]
        tmode - temporal mode number in band f, in [0,ntmodes[f]-1]
        smode - spatial mode number in band f, in [0,nsmodes[f]-1]
        w     - workspace

Return: index in [0,p_complex-1]
*/

size_t
invert_coeff_idx(const size_t f, const size_t tmode, const size_t smode,
                 const invert_workspace * w)
{
  invert_tmode_workspace * tmode_p = w->tmode_workspace_p;
  invert_smode_workspace * smode_p = w->smode_workspace_p;
  return w->mode_idx[f] + CIDX2(tmode, tmode_p->nmodes[f], smode, smode_p->nmodes[f]);
}

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
