/*
 * invert_smode.h
 */

#ifndef INCLUDED_invert_smode_h
#define INCLUDED_invert_smode_h

#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_complex.h>

#include <magfield/magfield_eval.h>

/* temporal modes data structure */
typedef struct
{
  size_t nfreq;          /* number of frequency bands */
  size_t *nmodes;        /* number of spatial modes in each frequency band, length nfreq */
  size_t modes_tot;      /* total number of modes in all frequency bands */
  size_t *mode_idx;      /* array of length nfreq, mode_idx[i] = index of start of frequency band i */

  size_t max_threads;    /* total threads available */
  size_t plm_size;       /* size of (l,m) arrays */
  size_t mmax;           /* maximum spherical harmonic order */

  double *Plm;             /* Legendre functions, plm_size * max_threads */
  double *dPlm;            /* Legendre function derivatives, plm_size * max_threads */
  complex double *expmphi; /* exp(i m phi) values, (mmax+1) * max_threads */
  complex double *qlmr;    /* q_l^m(r) values, plm_size * max_threads */
  complex double *plmr;    /* p_l^m(r) values, plm_size * max_threads */
  complex double *drplmr;  /* d/dr [r p_l^m(r)] values, plm_size * max_threads */
  gsl_interp_accel **acc;  /* accelerator objects, max_threads */

  /*
   * modes[i] is a matrix of size 3*N-by-T containing all spatial
   * modes for frequency band i, one mode per column. T is the total number of
   * modes, but we will only use nmodes[i] of them. N is nr*nlm
   */
  gsl_matrix_complex ** modes_U;

  magfield_eval_workspace ** magfield_eval_p;
} invert_smode_workspace;

/*
 * Prototypes
 */

invert_smode_workspace *invert_smode_alloc(const size_t nfreq, const size_t nmodes[]);
void invert_smode_free(invert_smode_workspace *w);
int invert_smode_precompute(const double r, const double theta, const double phi, invert_smode_workspace * w);
int invert_smode_get(const double r, const double theta, const double phi, const size_t f, const size_t mode,
                     gsl_complex Phi[3], invert_smode_workspace * w);
int invert_smode_get_J(const double r, const double theta, const double phi, const size_t f, const size_t mode,
                       gsl_complex Phi[3], invert_smode_workspace * w);

#endif /* INCLUDED_invert_smode_h */
