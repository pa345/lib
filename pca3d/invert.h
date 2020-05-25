/*
 * invert.h
 */

#ifndef INCLUDED_invert_h
#define INCLUDED_invert_h

#include <stdarg.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_complex.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_spatwt.h>
#include <mainlib/ml_track_weight.h>
#include <bspline2/gsl_bspline2.h>

#include "invert_data.h"
#include "invert_tmode.h"
#include "invert_smode.h"

/* cast gsl_complex to complex double */
#define CastComplex(x) ((complex double) CMPLX(x.dat[0], x.dat[1]))

#define INVERT_MAX_BUFFER          2048

/*
 * approximate matrix size in bytes for precomputing J^T J; each
 * thread gets its own matrix (200 MB)
 */
#define INVERT_MATRIX_SIZE    (2e8)

#define INVERT_DEBUG              1

typedef struct
{
  double R;                             /* reference radius (km) */
  size_t nmax;                          /* total nmax including core and crustal field */
  size_t nsat;                          /* number of satellites */
  size_t nfreq;                         /* number of frequency bands */
  char tmode_file[INVERT_MAX_BUFFER];   /* temporal mode file */

  size_t max_iter;                      /* number of robust iterations */

  int use_weights;                      /* use weights in the fitting */

  int regularize;                       /* regularize the solution vector */

  double weight_X;                      /* relative weighting for X component */
  double weight_Y;                      /* relative weighting for Y component */
  double weight_Z;                      /* relative weighting for Z component */
  double weight_F;                      /* relative weighting for F component */
  double weight_DXDT;                   /* relative weighting for dX/dt component */
  double weight_DYDT;                   /* relative weighting for dY/dt component */
  double weight_DZDT;                   /* relative weighting for dZ/dt component */
  double weight_DX;                     /* relative weighting for DX component */
  double weight_DY;                     /* relative weighting for DY component */
  double weight_DZ;                     /* relative weighting for DZ component */

  double qdlat_fit_cutoff;              /* QD latitude separating high and mid-latitudes for fitting components (degrees) */

  /* synthetic data parameters */
  int synth_data;                       /* replace real data with synthetic for testing */
  int synth_noise;                      /* add gaussian noise to synthetic model */
  size_t synth_nmin;                    /* minimum spherical harmonic degree for synthetic model */

  invert_data_workspace *invert_data_p; /* satellite data */
} invert_parameters;

typedef struct
{
  size_t nsat;        /* number of different satellites */
  size_t nfreq;       /* number of frequency bins */

  invert_parameters params;

  size_t *smode_idx;  /* length nfreq; smode_idx[i] = index of start of frequency band i */
  size_t nsmodes;     /* total number of spatial modes in all frequency bands */

  size_t *tmode_idx;  /* length nfreq; tmode_idx[i] = index of start of frequency band i */
  size_t ntmodes;     /* number of temporal modes for all frequency bands */

  size_t *mode_idx;  /* length nfreq; mode_idx[i] = index of start of frequency band i */

  size_t p;          /* number of real model coefficients */
  size_t p_complex;  /* number of complex model coefficients (p/2) */

  size_t nobs_cnt;

  double *t;        /* data timestamps in units of years */

  double t0_data;   /* time of first data input (CDF_EPOCH) */

  double R;         /* reference radius (km) */

  /*
   * The model coefficients are partitioned as follows:
   *
   * c = [ MF | SV | SA | Euler | External | fluxcal ]
   */
  gsl_vector *c;       /* model coefficients */
  gsl_vector *c_copy;  /* model coefficients in physical units */

  gsl_matrix *covar;   /* coefficient covariance matrix */

  size_t niter;        /* number of robust LS iterations */

  /* nonlinear least squares parameters */
  gsl_vector *wts_spatial; /* spatial weights, nres-by-1 */
  gsl_vector *wts_robust;  /* robust weights, nres-by-1 */
  gsl_vector *wts_final;   /* final weights (robust x spatial), nres-by-1 */
  gsl_vector *sqrt_wts_final; /* sqrt(wts_final), length nres */
  gsl_multifit_nlinear_workspace *multifit_nlinear_p;
  gsl_multilarge_nlinear_workspace *nlinear_workspace_p;
  size_t ndata;            /* number of unique data points in LS system */
  size_t nres_tot;         /* total number of residuals to minimize, including regularization terms */
  size_t nres;             /* number of residuals to minimize (data only) */
  size_t nres_vec;         /* number of vector residuals to minimize */
  size_t nres_vec_SV;      /* number of secular variation vector residuals to minimize */
  size_t nres_vec_grad;    /* number of vector gradient residuals to minimize */
  size_t data_block;       /* maximum observations to accumulate at once in LS system */
  size_t data_block_tot;   /* maximum observations to accumulate at once in J_int matrix */

  /* regularization parameters */
  gsl_spmatrix *L;         /* regularization matrix Cholesky factor, p-by-p */
  gsl_spmatrix *Lambda;    /* regularization matrix (L * L'), p-by-p */
  int old_fdf;             /* use multifit instead of multilarge */

  /*
   * The Jacobian is organized as follows:
   *
   * J = [ J_mf    | J_sv    | J_sa    | J_euler(x) | J_ext(x) ] vector
   *     [ J_mf(x) | J_sv(x) | J_sa(x) |     0      | J_ext(x) ] scalar
   *
   * J_mf, J_sv, and J_sa are constant for vector
   * residuals, and depend on the model parameters x
   * for scalar residuals. J_euler is 0 for scalar
   * residuals and depends on x for vector.
   * J_ext depends on x for both vector and scalar residuals.
   * J_euler and J_ext have significant sparse structure.
   *
   * For each iteration, we need to compute J^T J. This is
   * organized as:
   *
   * J^T J = [ JTJ_11 |    x   |    x   ]
   *         [ JTJ_21 | JTJ_22 |    x   ]
   *         [ JTJ_31 | JTJ_32 | JTJ_33 ]
   *
   * where we only need to compute the lower triangle since
   * the matrix is symmetric
   */
  gsl_matrix *JTJ_vec;     /* J_mf^T J_mf for vector measurements, p-by-p */
  gsl_matrix *choleskyL;   /* Cholesky factor for JTJ_vec if using linear system, p-by-p */

  size_t max_threads;      /* maximum number of threads/processors available */
  gsl_matrix **omp_J;      /* max_threads matrices, each 4*data_block-by-p */
  gsl_vector **omp_f;      /* max_threads vectors, each 4*data_block-by-1 */
  gsl_matrix **omp_B;      /* max_threads matrices, each 3-by-p */
  size_t *omp_rowidx;      /* row indices for omp_J */

  int lls_solution;        /* 1 if inverse problem is linear (no scalar residuals or Euler angles) */

  gsl_vector *fvec;        /* residual vector for robust weights */
  gsl_vector *wfvec;       /* weighted residual vector */
  gsl_multifit_robust_workspace *robust_workspace_p;
  gsl_multilarge_linear_workspace *multilarge_linear_p;

  invert_tmode_workspace *tmode_workspace_p;
  invert_smode_workspace *smode_workspace_p;
  invert_data_workspace *data_workspace_p;
  track_weight_workspace *weight_workspace_p;
  spatwt_workspace *spatwtMF_workspace_p; /* spatial weights for observatory MF measurements */
  spatwt_workspace *spatwtSV_workspace_p; /* spatial weights for observatory SV measurements */
} invert_workspace;

/*
 * Prototypes
 */

invert_workspace *invert_alloc(const invert_parameters *params);
void invert_free(invert_workspace *w);
int invert_init_params(invert_parameters * params);
invert_workspace *invert_copy(const invert_workspace *w);
int invert_init(invert_workspace *w);
int invert_calc_linear(gsl_vector *c, invert_workspace *w);
int invert_calc_nonlinear(gsl_vector *c, invert_workspace *w);
int invert_calc_JTJ_multilarge(const gsl_vector *c, gsl_matrix * JTJ, invert_workspace *w);
gsl_vector *invert_residual(const gsl_vector *c, invert_workspace *w);
int invert_reset(invert_workspace *w);
int invert_eval(const double t, const double r, const double theta, const double phi,
                double B[4], invert_workspace *w);
int invert_write(const char *filename, invert_workspace *w);
int invert_read(const char *filename, gsl_vector *c);
int invert_write_matrix(const char *filename, invert_workspace *w);
int invert_write_rhs(const char *filename, invert_workspace *w);
int invert_write_ascii(const char *filename, const gsl_vector * c, invert_workspace *w);
size_t invert_coeff_idx(const size_t f, const size_t tmode, const size_t smode, const invert_workspace * w);
int invert_debug(const char *format, ...);

#endif /* INCLUDED_invert_h */
