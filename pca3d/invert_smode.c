/*
 * invert_smode.c
 *
 * Routines related to spatial modes in the least squares inversion
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_legendre.h>

#include <common/common.h>
#include <common/interp.h>
#include <satdata/satdata.h>
#include <magfield/magfield.h>
#include <magfield/magfield_eval.h>

#include "invert_smode.h"
#include "pca3d.h"
#include "io.h"

static size_t smode_idx(const size_t f, const size_t mode, const invert_smode_workspace * w);
static int smode_fill_arrays(const magfield_params * params, const gsl_vector_complex * u,
                             gsl_vector_complex * qtcoeff, gsl_vector_complex * qcoeff,
                             gsl_vector_complex * pcoeff);

/*
invert_smode_alloc()
  Allocate invert_smode_workspace

Inputs: nfreq  - number of frequency bands
        nmodes - number of modes per frequency band, length nfreq

Return: pointer to workspace
*/

invert_smode_workspace *
invert_smode_alloc(const size_t nfreq, const size_t nmodes[])
{
  invert_smode_workspace *w;
  magfield_params params;
  gsl_vector_complex * qtcoeff;
  gsl_vector_complex * qcoeff;
  gsl_vector_complex * pcoeff;
  size_t N, i, j, k;

  w = calloc(1, sizeof(invert_smode_workspace));
  if (!w)
    return 0;

  w->nfreq = nfreq;
  w->max_threads = (size_t) omp_get_max_threads();

  w->nmodes = malloc(nfreq * sizeof(size_t));
  if (!w->nmodes)
    {
      invert_smode_free(w);
      return 0;
    }

  w->modes_U = calloc(nfreq, sizeof(gsl_matrix_complex *));
  if (!w->modes_U)
    {
      invert_smode_free(w);
      return 0;
    }

  w->mode_idx = malloc(nfreq * sizeof(size_t));

  /* read spatial mode matrices and keep count of total modes */
  w->modes_tot = 0;
  for (i = 0; i < nfreq; ++i)
    {
      char buf[1024];

      w->nmodes[i] = nmodes[i];

      w->mode_idx[i] = w->modes_tot;
      w->modes_tot += nmodes[i];

      /*XXX this is currently fixed for the 1cpd file */
      sprintf(buf, "%s_%d", PCA3D_STAGE3B_U, 2);
      w->modes_U[i] = pca3d_read_matrix_complex(buf);
      if (!w->modes_U[i])
        {
          invert_smode_free(w);
          fprintf(stderr, "invert_smode_alloc: failed to read matrix %s\n", buf);
          return 0;
        }
    }

  /* allocate magfield_eval workspaces */
  params = pca3d_read_params(PCA3D_STAGE1B_DATA);
  w->magfield_eval_p = malloc(w->modes_tot * w->max_threads * sizeof(magfield_eval_workspace *));

  N = w->modes_U[0]->size1 / 3;
  qtcoeff = gsl_vector_complex_alloc(N);
  qcoeff = gsl_vector_complex_alloc(N);
  pcoeff = gsl_vector_complex_alloc(N);

  for (i = 0; i < nfreq; ++i)
    {
      for (j = 0; j < w->nmodes[i]; ++j)
        {
          size_t sidx = smode_idx(i, j, w);
          gsl_vector_complex_const_view v = gsl_matrix_complex_const_column(w->modes_U[i], j);

          /* separate U(:,j) into q~, q, p coefficients */
          smode_fill_arrays(&params, &v.vector, qtcoeff, qcoeff, pcoeff);

          for (k = 0; k < w->max_threads; ++k)
            {
              size_t idx = CIDX2(sidx, w->modes_tot, k, w->max_threads);
              w->magfield_eval_p[idx] = magfield_eval_alloc(&params);
              magfield_eval_init(qtcoeff, qcoeff, pcoeff, w->magfield_eval_p[idx]);
            }
        }
    }

  gsl_vector_complex_free(qtcoeff);
  gsl_vector_complex_free(qcoeff);
  gsl_vector_complex_free(pcoeff);

  w->plm_size = gsl_sf_legendre_array_n(params.lmax);
  w->mmax = params.mmax;

  w->Plm = malloc(w->plm_size * w->max_threads * sizeof(double));
  w->dPlm = malloc(w->plm_size * w->max_threads * sizeof(double));
  w->qlmr = malloc(w->plm_size * w->modes_tot * w->max_threads * sizeof(complex double));
  w->plmr = malloc(w->plm_size * w->modes_tot * w->max_threads * sizeof(complex double));
  w->drplmr = malloc(w->plm_size * w->modes_tot * w->max_threads * sizeof(complex double));
  w->expmphi = malloc((w->mmax + 1) * w->max_threads * sizeof(complex double));

  w->acc = malloc(w->max_threads * sizeof(gsl_interp_accel *));

  for (i = 0; i < w->max_threads; ++i)
    w->acc[i] = gsl_interp_accel_alloc();

  return w;
}

void
invert_smode_free(invert_smode_workspace *w)
{
  if (w->nmodes)
    free(w->nmodes);

  if (w->mode_idx)
    free(w->mode_idx);

  if (w->modes_U)
    {
      size_t i;

      for (i = 0; i < w->nfreq; ++i)
        {
          if (w->modes_U[i])
            gsl_matrix_complex_free(w->modes_U[i]);
        }

      free(w->modes_U);
    }

  if (w->magfield_eval_p)
    {
      size_t i;

      for (i = 0; i < w->modes_tot * w->max_threads; ++i)
        {
          if (w->magfield_eval_p[i])
            magfield_eval_free(w->magfield_eval_p[i]);
        }

      free(w->magfield_eval_p);
    }

  if (w->Plm)
    free(w->Plm);

  if (w->dPlm)
    free(w->dPlm);

  if (w->qlmr)
    free(w->qlmr);

  if (w->plmr)
    free(w->plmr);

  if (w->drplmr)
    free(w->drplmr);

  if (w->expmphi)
    free(w->expmphi);

  if (w->acc)
    {
      size_t i;

      for (i = 0; i < w->max_threads; ++i)
        {
          if (w->acc[i])
            gsl_interp_accel_free(w->acc[i]);
        }

      free(w->acc);
    }

  free(w);
}

/*
invert_smode_precompute()
  Prepare to calculate spatial modes for a given (r,theta,phi) by
initializing magfield arrays
*/

int
invert_smode_precompute(const double r, const double theta, const double phi, invert_smode_workspace * w)
{
  int status = 0;
  int thread_id = omp_get_thread_num();
  double * Plm = w->Plm + thread_id * w->plm_size;
  double * dPlm = w->dPlm + thread_id * w->plm_size;
  complex double * expmphi = w->expmphi + thread_id * (w->mmax + 1);
  size_t i, j;

  status += magfield_eval_B_precompute_theta(theta, Plm, dPlm, w->magfield_eval_p[0]);
  status += magfield_eval_B_precompute_phi(phi, expmphi, w->magfield_eval_p[0]);

  for (i = 0; i < w->nfreq; ++i)
    {
      for (j = 0; j < w->nmodes[i]; ++j)
        {
          size_t sidx = smode_idx(i, j, w);
          size_t idx = CIDX2(sidx, w->modes_tot, thread_id, w->max_threads);
          complex double * qlmr = w->qlmr + idx * w->plm_size;
          complex double * plmr = w->plmr + idx * w->plm_size;
          complex double * drplmr = w->drplmr + idx * w->plm_size;

          status = magfield_eval_B_precompute_r(r, qlmr, plmr, drplmr, w->acc[thread_id], w->magfield_eval_p[idx]);
          if (status)
            return status;
        }
    }

  return GSL_SUCCESS;
}

/*
invert_smode_get()
  Return value of spatial mode for a specified frequency band and mode number

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        f     - frequency band, in [0,nfreq-1]
        mode  - mode number in frequency band f, in [0,nmodes[f]-1]
        Phi   - (output) Phi_{f,mode}(r,theta,phi)
                Phi[0] = NEC X component
                Phi[1] = NEC Y component
                Phi[2] = NEC Z component
        w     - workspace

Return: success/error

Notes:
1) invert_smode_precompute() must be called prior to this function to setup magfield correctly
*/

int
invert_smode_get(const double r, const double theta, const double phi, const size_t f, const size_t mode,
                 gsl_complex Phi[3], invert_smode_workspace * w)
{
  (void) phi;

  if (f >= w->nfreq)
    {
      GSL_ERROR ("invalid frequency band", GSL_EINVAL);
    }
  else if (mode >= w->nmodes[f])
    {
      GSL_ERROR ("invalid mode number", GSL_EINVAL);
    }
  else
    {
      int thread_id = omp_get_thread_num();
      size_t sidx = smode_idx(f, mode, w);
      size_t idx = CIDX2(sidx, w->modes_tot, thread_id, w->max_threads);
      magfield_eval_workspace * magfield_eval_p = w->magfield_eval_p[idx];
      double * Plm = w->Plm + thread_id * w->plm_size;
      double * dPlm = w->dPlm + thread_id * w->plm_size;
      complex double * expmphi = w->expmphi + thread_id * (w->mmax + 1);
      complex double * qlmr = w->qlmr + idx * w->plm_size;
      complex double * plmr = w->plmr + idx * w->plm_size;
      complex double * drplmr = w->drplmr + idx * w->plm_size;
      double B[4];

#if 1
      magfield_eval_B_compute(r, theta, Plm, dPlm, expmphi, qlmr, plmr, drplmr, B, magfield_eval_p);
#else
      magfield_eval_B(r, theta, phi, B, magfield_eval_p);
#endif

      Phi[0] = gsl_complex_rect(-B[1], 0.0);
      Phi[1] = gsl_complex_rect(B[2], 0.0);
      Phi[2] = gsl_complex_rect(-B[0], 0.0);

      return GSL_SUCCESS;
    }
}

/*
invert_smode_get_J()
  Return value of current density spatial mode for a specified frequency band and mode number

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        f     - frequency band, in [0,nfreq-1]
        mode  - mode number in frequency band f, in [0,nmodes[f]-1]
        Phi   - (output) Phi_{f,mode}(r,theta,phi)
                Phi[0] = NEC X component
                Phi[1] = NEC Y component
                Phi[2] = NEC Z component
        w     - workspace

Return: success/error

Notes:
1) invert_smode_precompute() must be called prior to this function to setup magfield correctly
*/

int
invert_smode_get_J(const double r, const double theta, const double phi, const size_t f, const size_t mode,
                   gsl_complex Phi[3], invert_smode_workspace * w)
{
  (void) phi;

  if (f >= w->nfreq)
    {
      GSL_ERROR ("invalid frequency band", GSL_EINVAL);
    }
  else if (mode >= w->nmodes[f])
    {
      GSL_ERROR ("invalid mode number", GSL_EINVAL);
    }
  else
    {
      int thread_id = omp_get_thread_num();
      size_t sidx = smode_idx(f, mode, w);
      size_t idx = CIDX2(sidx, w->modes_tot, thread_id, w->max_threads);
      magfield_eval_workspace * magfield_eval_p = w->magfield_eval_p[idx];
      double * Plm = w->Plm + thread_id * w->plm_size;
      double * dPlm = w->dPlm + thread_id * w->plm_size;
      complex double * expmphi = w->expmphi + thread_id * (w->mmax + 1);
      complex double * qlmr = w->qlmr + idx * w->plm_size;
      complex double * plmr = w->plmr + idx * w->plm_size;
      complex double * drplmr = w->drplmr + idx * w->plm_size;
      double J[3];

      magfield_eval_J(r, theta, phi, J, magfield_eval_p);

      Phi[0] = gsl_complex_rect(-J[1], 0.0);
      Phi[1] = gsl_complex_rect(J[2], 0.0);
      Phi[2] = gsl_complex_rect(-J[0], 0.0);

      return GSL_SUCCESS;
    }
}

/*
smode_idx()
  Return index in [0,w->modes_tot-1] corresponding to frequency band 'f' and
mode number 'mode'
*/

static size_t
smode_idx(const size_t f, const size_t mode, const invert_smode_workspace * w)
{
  return w->mode_idx[f] + mode;
}

static int
smode_fill_arrays(const magfield_params * params, const gsl_vector_complex * u,
                  gsl_vector_complex * qtcoeff, gsl_vector_complex * qcoeff,
                  gsl_vector_complex * pcoeff)
{
  const size_t N = qtcoeff->size;
  const size_t nlm = mie_lmidx(params->lmax, params->mmax, params->mmax) + 1;
  size_t ir;

  for (ir = 0; ir < params->nr; ++ir)
    {
      size_t l;

      for (l = params->lmin; l <= params->lmax; ++l)
        {
          int M = (int) GSL_MIN(l, params->mmax);
          int m;

          for (m = 0; m <= M; ++m)
            {
              size_t lmidx = magfield_lmidx(l, m, params->mmax);
              size_t uidx = CIDX2(ir, params->nr, lmidx, nlm);
              /*size_t midx = MAG_COEFIDX(ir, lmidx, w);*/
              size_t midx = CIDX2(ir, params->nr, lmidx, nlm); /* XXX MAG_COEFIDX() */

              gsl_complex qt = gsl_vector_complex_get(u, uidx);
              gsl_complex p = gsl_vector_complex_get(u, uidx + N);
              gsl_complex q = gsl_vector_complex_get(u, uidx + 2*N);

              gsl_vector_complex_set(qtcoeff, midx, qt);
              gsl_vector_complex_set(pcoeff, midx, p);
              gsl_vector_complex_set(qcoeff, midx, q);
            }
        }
    }

  return GSL_SUCCESS;
}
