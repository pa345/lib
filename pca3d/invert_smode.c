/*
 * invert_smode.c
 *
 * Routines related to spatial modes in the least squares inversion
 *
 * Memory management
 * -----------------
 *
 * The loops and indexing are organized so the following code will execute
 * without errors:
 *
 * idx = 0;
 * for (i = 0; i < w->nfreq; ++i)
 *   {
 *     for (j = 0; j < w->nmodes[i]; ++j)
 *       {
 *         size_t sidx = smode_idx(i, j, w);
 *         assert(idx == sidx);
 *         ++idx;
 *       }
 *   }
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_legendre.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_interp.h>
#include <mainlib/ml_satdata.h>

#include <magfield/magfield.h>
#include <magfield/magfield_eval_complex.h>

#include "invert_smode.h"
#include "pca3d.h"
#include "io.h"

static size_t smode_idx(const size_t f, const size_t mode, const invert_smode_workspace * w);
static int smode_fill_arrays(const gsl_vector_complex * u,
                             gsl_vector_complex * qtcoeff, gsl_vector_complex * qcoeff,
                             gsl_vector_complex * pcoeff, const magfield_eval_complex_workspace * w);
static int smode_precompute_r(const int thread_id, const double r, invert_smode_workspace * w);
static int smode_precompute_theta(const int thread_id, const double theta, invert_smode_workspace * w);
static int smode_precompute_phi(const int thread_id, const double phi, invert_smode_workspace * w);
static int smode_precompute_J_r(const int thread_id, const double r, invert_smode_workspace * w);
static int smode_precompute_J_theta(const int thread_id, const double theta, invert_smode_workspace * w);
static int smode_precompute_J_phi(const int thread_id, const double phi, invert_smode_workspace * w);

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
  gsl_vector *band_freqs;
  size_t N, i, j;

  w = calloc(1, sizeof(invert_smode_workspace));
  if (!w)
    return 0;

  w->nfreq = nfreq;
  w->max_threads = (size_t) omp_get_max_threads();

  w->freqs = malloc(nfreq * sizeof(double));
  w->nmodes = malloc(nfreq * sizeof(size_t));
  if (!w->freqs || !w->nmodes)
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

  /* read band frequencies */
  band_freqs = pca3d_read_vector(PCA3D_STAGE3B_BANDS);
  if (!band_freqs)
    {
      invert_smode_free(w);
      return 0;
    }

  if (band_freqs->size != nfreq)
    {
      fprintf(stderr, "invert_smode_alloc: mismatch in number of frequency bands\n");
      invert_smode_free(w);
      return 0;
    }

  /* read spatial mode matrices and keep count of total modes */
  w->modes_tot = 0;
  for (i = 0; i < nfreq; ++i)
    {
      char buf[1024];

      w->freqs[i] = gsl_vector_get(band_freqs, i);
      w->nmodes[i] = nmodes[i];

      w->mode_idx[i] = w->modes_tot;
      w->modes_tot += nmodes[i];

      sprintf(buf, "%s_%zu", PCA3D_STAGE3B_U, i + 1);
      w->modes_U[i] = pca3d_read_matrix_complex(buf);
      if (!w->modes_U[i])
        {
          invert_smode_free(w);
          fprintf(stderr, "invert_smode_alloc: failed to read matrix %s\n", buf);
          return 0;
        }
    }

  /* allocate magfield_eval_complex workspaces */
  params = pca3d_read_params(PCA3D_STAGE1B_DATA);
  w->magfield_eval_p = malloc(w->modes_tot * sizeof(magfield_eval_complex_workspace *));

  N = w->modes_U[0]->size1 / 3;
  qtcoeff = gsl_vector_complex_alloc(N);
  qcoeff = gsl_vector_complex_alloc(N);
  pcoeff = gsl_vector_complex_alloc(N);

  for (i = 0; i < nfreq; ++i)
    {
      fprintf(stderr, "invert_smode_alloc: freq %zu/%zu\n", i + 1, nfreq);
      for (j = 0; j < w->nmodes[i]; ++j)
        {
          size_t sidx = smode_idx(i, j, w);
          gsl_vector_complex_const_view v = gsl_matrix_complex_const_column(w->modes_U[i], j);

          /* allocate complex evaluation workspace */
          w->magfield_eval_p[sidx] = magfield_eval_complex_alloc(&params);

          /* separate U(:,j) into q~, q, p coefficients */
          smode_fill_arrays(&v.vector, qtcoeff, qcoeff, pcoeff, w->magfield_eval_p[sidx]);

          /* initialize magfield_eval_complex workspace for this mode U(:,j) */
          magfield_eval_complex_init(qtcoeff, qcoeff, pcoeff, w->magfield_eval_p[sidx]);
        }
    }

  gsl_vector_complex_free(qtcoeff);
  gsl_vector_complex_free(qcoeff);
  gsl_vector_complex_free(pcoeff);

  w->plm_size = gsl_sf_legendre_array_n(params.lmax);
  w->nlm = w->magfield_eval_p[0]->nlm;
  w->mmax = params.mmax;

  w->Plm = malloc(w->plm_size * w->max_threads * sizeof(double));
  w->dPlm = malloc(w->plm_size * w->max_threads * sizeof(double));
  w->qtlmr = malloc(w->nlm * w->modes_tot * w->max_threads * sizeof(complex double));
  w->qlmr = malloc(w->nlm * w->modes_tot * w->max_threads * sizeof(complex double));
  w->plmr = malloc(w->nlm * w->modes_tot * w->max_threads * sizeof(complex double));
  w->drplmr = malloc(w->nlm * w->modes_tot * w->max_threads * sizeof(complex double));
  w->dqlmr = malloc(w->nlm * w->modes_tot * w->max_threads * sizeof(complex double));
  w->expmphi = malloc((w->mmax + 1) * w->max_threads * sizeof(complex double));

  w->acc = malloc(w->max_threads * sizeof(gsl_interp_accel *));

  for (i = 0; i < w->max_threads; ++i)
    w->acc[i] = gsl_interp_accel_alloc();

  gsl_vector_free(band_freqs);

  return w;
}

void
invert_smode_free(invert_smode_workspace *w)
{
  if (w->freqs)
    free(w->freqs);

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

      for (i = 0; i < w->modes_tot; ++i)
        {
          if (w->magfield_eval_p[i])
            magfield_eval_complex_free(w->magfield_eval_p[i]);
        }

      free(w->magfield_eval_p);
    }

  if (w->Plm)
    free(w->Plm);

  if (w->dPlm)
    free(w->dPlm);

  if (w->qtlmr)
    free(w->qtlmr);

  if (w->qlmr)
    free(w->qlmr);

  if (w->plmr)
    free(w->plmr);

  if (w->drplmr)
    free(w->drplmr);

  if (w->dqlmr)
    free(w->dqlmr);

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

Inputs: thread_id - thread id
        r         - radius (km)
        theta     - co-latitude (radians)
        phi       - longitude (radians)
        w         - workspace
*/

int
invert_smode_precompute(const int thread_id, const double r, const double theta, const double phi, invert_smode_workspace * w)
{
  int status;

  status = smode_precompute_r(thread_id, r, w);
  status += smode_precompute_theta(thread_id, theta, w);
  status += smode_precompute_phi(thread_id, phi, w);

  return status;
}

/*
invert_smode_precompute_J()
  Prepare to calculate J spatial modes for a given (r,theta,phi) by
initializing magfield arrays
*/

int
invert_smode_precompute_J(const int thread_id, const double r, const double theta, const double phi, invert_smode_workspace * w)
{
  int status;

  status = smode_precompute_J_r(thread_id, r, w);
  status += smode_precompute_J_theta(thread_id, theta, w);
  status += smode_precompute_J_phi(thread_id, phi, w);

  return status;
}

/*
invert_smode_get()
  Return value of spatial mode (B) for a specified frequency band and mode number

Inputs: thread_id - thread id
        r         - radius (km)
        theta     - colatitude (radians)
        phi       - longitude (radians)
        f         - frequency band, in [0,nfreq-1]
        mode      - mode number in frequency band f, in [0,nmodes[f]-1]
        Phi       - (output) Phi_{f,mode}(r,theta,phi)
                    Phi[0] = NEC X component
                    Phi[1] = NEC Y component
                    Phi[2] = NEC Z component
        w         - workspace

Return: success/error

Notes:
1) invert_smode_precompute() must be called prior to this function to setup magfield correctly
*/

int
invert_smode_get(const int thread_id, const double r, const double theta, const double phi, const size_t f, const size_t mode,
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
      const size_t thread_offset = thread_id * w->nlm * w->modes_tot; /* offset for this thread in qlmr,plmr,drplmr arrays */
      const size_t sidx = smode_idx(f, mode, w);
      const size_t idx = thread_offset + sidx * w->nlm;
      const magfield_eval_complex_workspace * magfield_eval_p = w->magfield_eval_p[sidx];
      const double * Plm = w->Plm + thread_id * w->plm_size;
      const double * dPlm = w->dPlm + thread_id * w->plm_size;
      const complex double * expmphi = w->expmphi + thread_id * (w->mmax + 1);
      const complex double * qlmr = w->qlmr + idx;
      const complex double * plmr = w->plmr + idx;
      const complex double * drplmr = w->drplmr + idx;
      complex double B[3];

      magfield_eval_complex_B_compute(r * 1.0e3, theta, Plm, dPlm, expmphi, qlmr, plmr, drplmr, B, magfield_eval_p);

      Phi[0] = gsl_complex_rect(-creal(B[1]), -cimag(B[1]));
      Phi[1] = gsl_complex_rect(creal(B[2]), cimag(B[2]));
      Phi[2] = gsl_complex_rect(-creal(B[0]), -cimag(B[0]));

      return GSL_SUCCESS;
    }
}

/*
invert_smode_get_J()
  Return value of current density spatial mode for a specified frequency band and mode number

Inputs: thread_id - thread id
        r         - radius (km)
        theta     - colatitude (radians)
        phi       - longitude (radians)
        f         - frequency band, in [0,nfreq-1]
        mode      - mode number in frequency band f, in [0,nmodes[f]-1]
        Phi       - (output) Phi_{f,mode}(r,theta,phi)
                    Phi[0] = NEC X component
                    Phi[1] = NEC Y component
                    Phi[2] = NEC Z component
        w         - workspace

Return: success/error

Notes:
1) invert_smode_precompute_J() must be called prior to this function to setup magfield correctly
*/

int
invert_smode_get_J(const int thread_id, const double r, const double theta, const double phi, const size_t f, const size_t mode,
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
      const size_t thread_offset = thread_id * w->nlm * w->modes_tot; /* offset for this thread in qlmr,plmr,drplmr arrays */
      const size_t sidx = smode_idx(f, mode, w);
      const size_t idx = thread_offset + sidx * w->nlm;
      magfield_eval_complex_workspace * magfield_eval_p = w->magfield_eval_p[sidx];
      double * Plm = w->Plm + thread_id * w->plm_size;
      double * dPlm = w->dPlm + thread_id * w->plm_size;
      complex double * expmphi = w->expmphi + thread_id * (w->mmax + 1);
      complex double * qtlmr = w->qtlmr + idx;
      complex double * qlmr = w->qlmr + idx;
      complex double * dqlmr = w->dqlmr + idx;
      complex double J[3];

      magfield_eval_complex_J_compute(r * 1.0e3, theta, Plm, dPlm, expmphi, qtlmr, qlmr, dqlmr, J, magfield_eval_p);

      Phi[0] = gsl_complex_rect(-creal(J[1]), -cimag(J[1]));
      Phi[1] = gsl_complex_rect(creal(J[2]), cimag(J[2]));
      Phi[2] = gsl_complex_rect(-creal(J[0]), -cimag(J[0]));

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

/*
smode_fill_arrays()
  For a given eigenmode, fill in the 'qtcoeff', 'pcoeff', and 'qcoeff' arrays

Inputs: u       - singular vector (3*N-by-1 where N = nr*nlm_complex)
        qtcoeff - (output) qt coefficients, length N
        qcoeff  - (output) q coefficients, length N
        pcoeff  - (output) p coefficients, length N
        w       - magfield_eval_complex workspace
*/

static int
smode_fill_arrays(const gsl_vector_complex * u,
                  gsl_vector_complex * qtcoeff, gsl_vector_complex * qcoeff,
                  gsl_vector_complex * pcoeff, const magfield_eval_complex_workspace * w)
{
  const size_t N = w->nr * w->nlm;
  size_t ir;

  for (ir = 0; ir < w->nr; ++ir)
    {
      size_t l;

      for (l = w->lmin; l <= w->lmax; ++l)
        {
          int M = (int) GSL_MIN(l, w->mmax);
          int m;

          for (m = -M; m <= M; ++m)
            {
              size_t lmidx = magfield_complex_lmidx(l, m, w->mmax);
              size_t uidx = CIDX2(ir, w->nr, lmidx, w->nlm);
              size_t midx = MAG_CMPLX_COEFIDX(ir, lmidx, w);

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

/*
invert_smode_print()
  Print magnetic field spatial modes to output directory
*/

int
invert_smode_print(const char * dir_prefix, invert_smode_workspace * w)
{
  int s = 0;
  int thread_id = omp_get_thread_num();
  const double r = R_EARTH_KM + 350.0;
  char buf[2048];
  size_t f;

  smode_precompute_r(thread_id, r, w);
  smode_precompute_J_r(thread_id, r, w);

  fprintf(stderr, "\n");

  for (f = 0; f < w->nfreq; ++f)
    {
      size_t mode;

      for (mode = 0; mode < w->nmodes[f]; ++mode)
        {
          FILE *fp;
          size_t i;
          double lat, lon;

          sprintf(buf, "%s/smode_%02zu_%02zu.txt", dir_prefix, f + 1, mode + 1);
          fp = fopen(buf, "w");
          if (!fp)
            continue;

          i = 1;
          fprintf(fp, "# Frequency:   %f [cpd]\n", w->freqs[f]);
          fprintf(fp, "# Mode number: %zu\n", mode);
          fprintf(fp, "# Radius:      %.2f [km]\n", r);
          fprintf(fp, "# Field %zu: Longitude (degrees)\n", i++);
          fprintf(fp, "# Field %zu: Latitude (degrees)\n", i++);
          fprintf(fp, "# Field %zu: Re B_r component of mode\n", i++);
          fprintf(fp, "# Field %zu: Re B_t component of mode\n", i++);
          fprintf(fp, "# Field %zu: Re B_p component of mode\n", i++);
          fprintf(fp, "# Field %zu: Im B_r component of mode\n", i++);
          fprintf(fp, "# Field %zu: Im B_t component of mode\n", i++);
          fprintf(fp, "# Field %zu: Im B_p component of mode\n", i++);
          fprintf(fp, "# Field %zu: Re J_r component of mode\n", i++);
          fprintf(fp, "# Field %zu: Re J_t component of mode\n", i++);
          fprintf(fp, "# Field %zu: Re J_p component of mode\n", i++);
          fprintf(fp, "# Field %zu: Im J_r component of mode\n", i++);
          fprintf(fp, "# Field %zu: Im J_t component of mode\n", i++);
          fprintf(fp, "# Field %zu: Im J_p component of mode\n", i++);

          for (lon = -180.0; lon <= 180.0; lon += 3.0)
            {
              double phi = lon * M_PI / 180.0;

              smode_precompute_phi(thread_id, phi, w);
              smode_precompute_J_phi(thread_id, phi, w);

              for (lat = -89.5; lat <= 89.5; lat += 3.0)
                {
                  double theta = M_PI / 2.0 - lat * M_PI / 180.0;
                  gsl_complex B[3], J[3];

                  smode_precompute_theta(thread_id, theta, w);
                  smode_precompute_J_theta(thread_id, theta, w);

                  invert_smode_get(thread_id, r, theta, phi, f, mode, B, w);
                  invert_smode_get_J(thread_id, r, theta, phi, f, mode, J, w);

                  fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                          lon,
                          lat,
                          -GSL_REAL(B[2]),
                          -GSL_REAL(B[0]),
                          GSL_REAL(B[1]),
                          -GSL_IMAG(B[2]),
                          -GSL_IMAG(B[0]),
                          GSL_IMAG(B[1]),
                          -GSL_REAL(J[2]),
                          -GSL_REAL(J[0]),
                          GSL_REAL(J[1]),
                          -GSL_IMAG(J[2]),
                          -GSL_IMAG(J[0]),
                          GSL_IMAG(J[1]));
                }

              fprintf(fp, "\n");
            }

          fclose(fp);
        }

      progress_bar(stderr, (double) (f + 1.0)  / (double) w->nfreq, 70);
    }

  progress_bar(stderr, 1.0, 70);

  return s;
}

/*
smode_precompute_r()
  Prepare to calculate spatial modes for a given r by
initializing magfield arrays

Inputs: thread_id - thread id
        r         - radius (km)
        w         - workspace
*/

static int
smode_precompute_r(const int thread_id, const double r, invert_smode_workspace * w)
{
  int status = 0;
  const double r_m = r * 1.0e3;
  const size_t thread_offset = thread_id * w->nlm * w->modes_tot; /* offset for this thread in qlmr,plmr,drplmr arrays */
  size_t j;

  /* for each mode, initialize qlmr, plmr, drplmr arrays */
  for (j = 0; j < w->modes_tot; ++j)
    {
      size_t idx = thread_offset + j * w->nlm;
      complex double * qlmr = w->qlmr + idx;
      complex double * plmr = w->plmr + idx;
      complex double * drplmr = w->drplmr + idx;

      status = magfield_eval_complex_B_precompute_r(r_m, qlmr, plmr, drplmr, w->acc[thread_id], w->magfield_eval_p[j]);
      if (status)
        return status;
    }

  return GSL_SUCCESS;
}

/*
smode_precompute_theta()
  Prepare to calculate spatial modes for a given theta by
initializing magfield arrays

Inputs: thread_id - thread id
        theta     - co-latitude (degrees)
        w         - workspace
*/

static int
smode_precompute_theta(const int thread_id, const double theta, invert_smode_workspace * w)
{
  int status = 0;
  double * Plm = w->Plm + thread_id * w->plm_size;
  double * dPlm = w->dPlm + thread_id * w->plm_size;

  status = magfield_eval_complex_precompute_theta(theta, Plm, dPlm, w->magfield_eval_p[0]);

  return status;
}

/*
smode_precompute_phi()
  Prepare to calculate spatial modes for a given phi by
initializing magfield arrays

Inputs: thread_id - thread id
        r         - radius (km)
        w         - workspace
*/

static int
smode_precompute_phi(const int thread_id, const double phi, invert_smode_workspace * w)
{
  int status;
  complex double * expmphi = w->expmphi + thread_id * (w->mmax + 1);

  status = magfield_eval_complex_precompute_phi(phi, expmphi, w->magfield_eval_p[0]);

  return status;
}

/*
smode_precompute_J_r()
  Prepare to calculate J spatial modes for a given r by
initializing magfield arrays
*/

static int
smode_precompute_J_r(const int thread_id, const double r, invert_smode_workspace * w)
{
  int status;
  const double r_m = r * 1.0e3;
  const size_t thread_offset = thread_id * w->nlm * w->modes_tot; /* offset for this thread in qtlmr,qlmr,dqlmr arrays */
  size_t j;

  /* for each mode, initialize qtlmr, qlmr, dqlmr arrays */
  for (j = 0; j < w->modes_tot; ++j)
    {
      size_t idx = thread_offset + j * w->nlm;
      complex double * qtlmr = w->qtlmr + idx;
      complex double * qlmr = w->qlmr + idx;
      complex double * dqlmr = w->dqlmr + idx;

      status = magfield_eval_complex_J_precompute_r(r_m, qtlmr, qlmr, dqlmr, w->acc[thread_id], w->magfield_eval_p[j]);
      if (status)
        return status;
    }

  return GSL_SUCCESS;
}

/*
smode_precompute_J_theta()
  Prepare to calculate J spatial modes for a given theta by
initializing magfield arrays
*/

static int
smode_precompute_J_theta(const int thread_id, const double theta, invert_smode_workspace * w)
{
  int status;
  double * Plm = w->Plm + thread_id * w->plm_size;
  double * dPlm = w->dPlm + thread_id * w->plm_size;

  status = magfield_eval_complex_precompute_theta(theta, Plm, dPlm, w->magfield_eval_p[0]);

  return status;
}

/*
smode_precompute_J_phi()
  Prepare to calculate J spatial modes for a given phi by
initializing magfield arrays
*/

static int
smode_precompute_J_phi(const int thread_id, const double phi, invert_smode_workspace * w)
{
  int status;
  complex double * expmphi = w->expmphi + thread_id * (w->mmax + 1);

  status = magfield_eval_complex_precompute_phi(phi, expmphi, w->magfield_eval_p[0]);

  return status;
}
