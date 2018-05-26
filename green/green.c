/*
 * green.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_legendre.h>

#include "green.h"

/*
green_alloc()
  Allocate Green's function workspace

Inputs: nmax - maximum spherical harmonic degree
        mmax - maximum spherical harmonic order
        R    - reference radius (km)
*/

green_workspace *
green_alloc(const size_t nmax, const size_t mmax, const double R)
{
  green_workspace *w;
  size_t plm_array_size = gsl_sf_legendre_array_n(nmax);

  if (mmax > nmax)
    {
      fprintf(stderr, "green_alloc: error: mmax > nmax\n");
      return 0;
    }

  w = calloc(1, sizeof(green_workspace));
  if (!w)
    return 0;

  /* total number of SH coefficients */
  w->nnm = green_calc_nnm(nmax, mmax);

  w->nmax = nmax;
  w->mmax = mmax;
  w->R = R;

  w->cosmphi = malloc((mmax + 1) * sizeof(double));
  w->sinmphi = malloc((mmax + 1) * sizeof(double));

  w->Plm = malloc(plm_array_size * sizeof(double));
  w->dPlm = malloc(plm_array_size * sizeof(double));
  if (!w->Plm || !w->dPlm)
    {
      green_free(w);
      return 0;
    }

  w->work = malloc(w->nnm * sizeof(double));
  if (!w->work)
    {
      green_free(w);
      return 0;
    }

  return w;
}

void
green_free(green_workspace *w)
{
  if (w->cosmphi)
    free(w->cosmphi);

  if (w->sinmphi)
    free(w->sinmphi);

  if (w->Plm)
    free(w->Plm);

  if (w->dPlm)
    free(w->dPlm);

  if (w->work)
    free(w->work);

  free(w);
}

/*
green_calc_int()
  Compute Green's functions for X,Y,Z spherical harmonic expansion. These
are simply the basis functions multiplying the g_{nm} and h_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        X     - (output) array of X Green's functions, size nnm
        Y     - (output) array of Y Green's functions, size nnm
        Z     - (output) array of Z Green's functions, size nnm
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi

2) X,Y,Z can be set to NULL if not desired
*/

int
green_calc_int(const double r, const double theta, const double phi,
               double *X, double *Y, double *Z, green_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = w->R / r;
  double term = ratio * ratio;     /* (a/r)^{n+2} */

  /* precompute cos(m phi) and sin(m phi) */
  for (m = 0; m <= (int) mmax; ++m)
    {
      w->cosmphi[m] = cos(m * phi);
      w->sinmphi[m] = sin(m * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_array(GSL_SF_LEGENDRE_SCHMIDT,
                              nmax, cost, w->Plm, w->dPlm);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(mmax, n);

      /* (a/r)^{n+2} */
      term *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m, w);
          size_t aidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* h_{nm} */
              if (X)
                X[cidx] = term * w->sinmphi[mabs] * w->dPlm[aidx] * (-sint);

              if (Y)
                Y[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[aidx];

              if (Z)
                Z[cidx] = -(n + 1.0) * term * w->sinmphi[mabs] * w->Plm[aidx];
            }
          else
            {
              /* g_{nm} */
              if (X)
                X[cidx] = term * w->cosmphi[mabs] * w->dPlm[aidx] * (-sint);

              if (Y)
                Y[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[aidx];

              if (Z)
                Z[cidx] = -(n + 1.0) * term * w->cosmphi[mabs] * w->Plm[aidx];
            }
        }
    }

  return s;
} /* green_calc_int() */

/*
green_calc_ext()
  Compute Green's functions for X,Y,Z spherical harmonic expansion due to
external current source. These are simply the basis functions multiplying
the k_{nm} and q_{nm} coefficients

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        X     - (output) array of X Green's functions, size nnm
        Y     - (output) array of Y Green's functions, size nnm
        Z     - (output) array of Z Green's functions, size nnm
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->dPlm
w->sinmphi
w->cosmphi

2) X,Y,Z can be set to NULL if not desired
*/

int
green_calc_ext(const double r, const double theta, const double phi,
               double *X, double *Y, double *Z, green_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  size_t n;
  int m;
  const double sint = sin(theta);
  const double cost = cos(theta);
  double ratio = r / w->R;
  double term = 1.0;     /* (r/a)^{n-1} */

  /* precompute cos(m phi) and sin(m phi) */
  for (m = 0; m <= (int) mmax; ++m)
    {
      w->cosmphi[m] = cos(m * phi);
      w->sinmphi[m] = sin(m * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cost,
                                  w->Plm, w->dPlm);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(mmax, n);

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* k_{nm} */
              if (X)
                X[cidx] = term * w->sinmphi[mabs] * w->dPlm[pidx];

              if (Y)
                Y[cidx] = -term / sint * mabs * w->cosmphi[mabs] * w->Plm[pidx];

              if (Z)
                Z[cidx] = (double) n * term * w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* q_{nm} */
              if (X)
                X[cidx] = term * w->cosmphi[mabs] * w->dPlm[pidx];

              if (Y)
                Y[cidx] = term / sint * mabs * w->sinmphi[mabs] * w->Plm[pidx];

              if (Z)
                Z[cidx] = (double) n * term * w->cosmphi[mabs] * w->Plm[pidx];
            }
        }

      /* (r/a)^{n-1} */
      term *= ratio;
    }

  return s;
} /* green_calc_ext() */

/*
green_potential_calc_ext()
  Compute Green's functions for potential spherical harmonic expansion due to
external current source:

V(r,theta,phi) = R sum_{nm} (r/R)^n [ q_{nm} cos(m phi) + k_{nm} sin(m phi) ] P_{nm}

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        V     - (output) array of Green's functions (coefficients of qnm and knm), size nnm
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->sinmphi
w->cosmphi
*/

int
green_potential_calc_ext(const double r, const double theta, const double phi,
                         double *V, green_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  size_t n;
  int m;
  const double cost = cos(theta);
  double ratio = r / w->R;
  double term = w->R * ratio; /* R (r/R)^n */

  /* precompute cos(m phi) and sin(m phi) */
  for (m = 0; m <= (int) mmax; ++m)
    {
      w->cosmphi[m] = cos(m * phi);
      w->sinmphi[m] = sin(m * phi);
    }

  /* compute associated legendres (XXX: don't need to compute dPlm) */
  gsl_sf_legendre_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cost, w->Plm);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(mmax, n);

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* k_{nm} */
              V[cidx] = term * w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* q_{nm} */
              V[cidx] = term * w->cosmphi[mabs] * w->Plm[pidx];
            }
        }

      /* R (r/R)^n */
      term *= ratio;
    }

  return s;
}

/*
green_Y_calc()
  Compute Green's functions for spherical harmonic expansion using
Schmidt-normalized harmonics

f(theta,phi) = sum_{nm} [ g_{nm} cos(m phi) + h_{nm} sin(m phi) ] P_{nm}

Inputs: theta - colatitude (radians)
        phi   - longitude (radians)
        Y     - (output) array of Green's functions (coefficients of gnm and hnm), size nnm
        w     - workspace

Notes:
1) On output, the following arrays are initialized
w->Plm
w->sinmphi
w->cosmphi
*/

int
green_Y_calc(const double theta, const double phi,
             double *Y, green_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  size_t n;
  int m;
  const double cost = cos(theta);

  /* precompute cos(m phi) and sin(m phi) */
  for (m = 0; m <= (int) mmax; ++m)
    {
      w->cosmphi[m] = cos(m * phi);
      w->sinmphi[m] = sin(m * phi);
    }

  /* compute associated legendres */
  gsl_sf_legendre_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cost, w->Plm);

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(mmax, n);

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t cidx = green_nmidx(n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);

          if (m < 0)
            {
              /* h_{nm} */
              Y[cidx] = w->sinmphi[mabs] * w->Plm[pidx];
            }
          else
            {
              /* g_{nm} */
              Y[cidx] = w->cosmphi[mabs] * w->Plm[pidx];
            }
        }
    }

  return s;
}

/*
green_nmidx()
  This function returns a unique index in [0,p-1] corresponding
to a given (n,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Things get a little more tricky when mmax != nmax, so the
base indices of each 'n' block are precomputed in _alloc and
stored for easy reference here. The offset of (n,m) in a given
'n' block is:

offset = m + min(n,mmax)

which defaults to the standard m + n for the case mmax = nmax

Inputs: n - SH degree in [1,nmax]
        m - SH order (-n <= m <= n)
        w - workspace

Return: index in [0,nnm-1]
*/

inline size_t
green_nmidx(const size_t n, const int m, const green_workspace *w)
{
  return green_idx(n, m, w->mmax);
}

/*
green_idx()
  This function returns a unique index in [0,p-1] corresponding
to a given (n,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Things get a little more tricky when mmax != nmax, so the
base indices of each 'n' block are precomputed in _alloc and
stored for easy reference here. The offset of (n,m) in a given
'n' block is:

offset = m + min(n,mmax)

which defaults to the standard m + n for the case mmax = nmax

Inputs: n    - SH degree in [1,nmax]
        m    - SH order (-mmax <= m <= mmax)
        mmax - maximum SH order

Return: index in [0,nnm-1]
*/

inline size_t
green_idx(const size_t n, const int m, const size_t mmax)
{
  size_t nmidx;

  if (n == 0)
    {
      fprintf(stderr, "green_idx: error: n = 0\n");
      return 0;
    }
  else if (abs(m) > (int) mmax)
    {
      fprintf(stderr, "green_idx: error: m = %d [mmax = %zu]\n", m, mmax);
      return 0;
    }

  if (n <= mmax)
    {
      size_t base = n * n; /* index of block for this n */
      int offset = m + n;  /* offset within block for this m */

      /* subtract 1 to exclude (0,0) coefficient */
      nmidx = base + offset - 1;
    }
  else
    {
      size_t base1 = (mmax + 1) * (mmax + 1);
      size_t base2 = (n - mmax - 1) * (2 * mmax + 1);
      int offset = m + (int)mmax;

      /* subtract 1 to exclude (0,0) coefficient */
      nmidx = base1 + base2 + offset - 1;
    }

  return nmidx;
}

size_t
green_nnm(const green_workspace *w)
{
  return w->nnm;
}

/*
green_calc_nnm()
  Compute total number of spherical harmonic coefficients
for a given nmax and mmax.

nnm = mmax * (mmax + 2) + (nmax - mmax) * (2*mmax + 1)

this formula starts at n = 1 (does not count n = 0 term)

Inputs: nmax - maximum SH degree
        mmax - maximum SH order

Return: number of total spherical harmonic coefficients
*/

size_t
green_calc_nnm(const size_t nmax, const size_t mmax)
{
  size_t nnm = mmax * (mmax + 2) + (nmax - mmax) * (2*mmax + 1);
  return nnm;
}

int
green_print_spectrum(const char *filename, const gsl_vector *c,
                     const green_workspace *w)
{
  const size_t nmax = w->nmax;
  size_t n;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_spectrum: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
  fprintf(fp, "# Field %zu: power (nT^2)\n", n++);

  for (n = 1; n <= nmax; ++n)
    {
      double sum = 0.0;
      int M = (int) GSL_MIN(w->mmax, n);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m, w);
          double gnm = gsl_vector_get(c, cidx);

          sum += gnm * gnm;
        }

      sum *= (n + 1.0);

      fprintf(fp, "%zu %.12e\n", n, sum);
    }

  fclose(fp);

  return 0;
}

int
green_print_spectrum_azim(const char *filename, const gsl_vector * c,
                          const green_workspace *w)
{
  size_t n;
  FILE *fp;
  
  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "green_print_spectrum_azim: cannot open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic degree n\n", n++);
  fprintf(fp, "# Field %zu: spherical harmonic order m\n", n++);
  fprintf(fp, "# Field %zu: m/n\n", n++);
  fprintf(fp, "# Field %zu: n + 0.25*(m/n)\n", n++);
  fprintf(fp, "# Field %zu: gnm^2 (nT^2)\n", n++);

  for (n = 1; n <= w->nmax; ++n)
    {
      int M = (int) GSL_MIN(w->mmax, n);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m, w);
          double gnm = gsl_vector_get(c, cidx);
          double ratio = (double)m / (double)n;

          fprintf(fp, "%4zu %4d %8.4f %8.4f %.12e\n",
                  n,
                  m,
                  ratio,
                  n + 0.25 * ratio,
                  gnm * gnm);
        }
    }

  fclose(fp);

  return 0;
}

/*
green_k2g()
  Convert external Gauss coefficients knm to internal coefficients
gnm using the relation

gnm = -(n / (n+1)) (b/R)^{2n + 1} knm

Inputs: b - radius of current shell (km)
        k - knm, size nnm
        g - (output) gnm, size nnm
        w - workspace

Return: success/error

Notes:
1) It is allowed for k = g for an in-place transform
*/

int
green_k2g(const double b, const gsl_vector *k, gsl_vector *g,
          const green_workspace *w)
{
  const size_t nmax = w->nmax;
  const size_t mmax = w->mmax;
  const double ratio = b / w->R;
  const double ratio_sq = ratio * ratio;
  double rfac = ratio * ratio_sq; /* (b/R)^3 */
  size_t n;

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;
      double nfac = -(double)n / (n + 1.0);

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = green_nmidx(n, m, w);
          double knm = gsl_vector_get(k, cidx);

          gsl_vector_set(g, cidx, nfac * rfac * knm);
        }

      /* (b/R)^{2n + 1} */
      rfac *= ratio_sq;
    }

  return 0;
}
