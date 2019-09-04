/*
 * pde.h
 * Patrick Alken
 */

#ifndef INCLUDED_pde_h
#define INCLUDED_pde_h

#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spmatrix.h>

#include <mainlib/ml_msynth.h>
#include <mainlib/ml_mageq.h>

#include <pde/gsl_pde2d.h>

#include "hwm.h"

#include "sigma.h"
#include "superlu.h"

/*****************************************
 * PDE solution parameters               *
 *****************************************/

/* initial eastward electric field in V/m */
#define EEF_PHI_0            (1.0e-3)

/* define to force zonal winds to be symmetric about dip equator */
#define PDE_SYMMETRIC_WINDS

/* minimum conductivity allowed in SI units */
#define PDE_MIN_CONDUCTIVITY (1.0e-10)

/*
 * define to taper sigma tensor to zero outside a window specified by
 * PDE_SIGMA_TAPER_RANGE
 */
#define PDE_SIGMA_TAPER
#define PDE_SIGMA_TAPER_RANGE (10.0 * M_PI / 180.0)

/* theta range to use for height integrated current calculation */
#define PDE_HI_THETA_MIN     (75.0 * M_PI / 180.0)
#define PDE_HI_THETA_MAX     (105.0 * M_PI / 180.0)

/*****************************************
 * Inversion parameters                  *
 *****************************************/

/* define to invert for electric field */
#define PDE_INVERT_E_FIELD

/* define to allow a DC shift in inversion */
#define PDE_ALLOW_DC_SHIFT

/*
 * define to constrain inversion to have solution match satellite
 * profile at 90 degrees
 */
#define PDE_USE_LSE

/* indices of basis variables for matrices */
#define IDX_R                0
#define IDX_THETA            1
#define IDX_PHI              2

/* PDE flags */
#define PDE_FLG_RHS_G1       (1 << 0) /* compute g1 part of rhs */
#define PDE_FLG_RHS_G2       (1 << 1) /* compute g2 part of rhs */

typedef struct
{
  size_t nr;         /* number of radial grid points */
  size_t ntheta;     /* number of theta grid points */
  double rmin;       /* minimum radius in m */
  double rmax;       /* maximum radius in m */
  double theta_min;  /* minimum theta in radians */
  double theta_max;  /* maximum theta in radians */
  double R;          /* reference radius in m */
  char *f107_file;   /* f10.7 data file */
} pde_parameters;

typedef struct
{
  size_t nr;         /* number of radial grid points */
  size_t ntheta;     /* number of theta grid points */

  double rmin;       /* minimum radius in m */
  double rmax;       /* maximum radius in m */
  double dr;         /* radius step size in m */
  double R;          /* reference radius in m */

  double theta_min;  /* minimum theta in radians */
  double theta_max;  /* maximum theta in radians */
  double dtheta;     /* theta step size in radians */

  double longitude;  /* current longitude */
  double lat_eq;     /* magnetic equator latitude */
  time_t t;          /* current time */

  double *r_grid;     /* array of r grid values (m), size nr */
  double *theta_grid; /* array of theta grid values (radians), size ntheta */

  double eej_angle;  /* angle EEJ makes with geographic eastward */

  gsl_matrix **sigma;

  double *zwind;    /* zonal wind (u_phi) */
  double *mwind;    /* meridional wind (u_theta) */
  double *vwind;    /* vertical wind (u_r) */

  double *merid;    /* meridional wind array for HWM */
  double *zonal;    /* zonal wind array for HWM */

  double *Br_main;  /* main field r component in T */
  double *Bt_main;  /* main field theta component in T */
  double *Bp_main;  /* main field phi component in T */
  double *Bf_main;  /* total main field intensity in T */

  gsl_spmatrix *S; /* sparse pde matrix */
  gsl_matrix *G;   /* rhs vectors, nr*ntheta-by-nrhs; G(:,1) = g1, G(:,2) = g2 */
  gsl_matrix *B;   /* rhs vectors with boundary conditions applied, nr*ntheta-by-nrhs */
  gsl_matrix *PSI; /* PDE solution, nr*ntheta-by-nrhs */
  double residual; /* residual ||A*psi - b|| */
  double rrnorm;   /* relative residual ||b - A*psi|| / ||b|| */
  double rcond;    /* reciprical condition number */
  size_t nrhs;     /* number of right hand sides */

  gsl_matrix *WR;      /* [sigma*(UxB)]_r matrix, nr*ntheta-by-nrhs */
  gsl_matrix *WTHETA;  /* [sigma*(UxB)]_theta matrix, nr*ntheta-by-nrhs */
  gsl_matrix *WPHI;    /* [sigma*(UxB)]_phi matrix, nr*ntheta-by-nrhs */
  gsl_matrix *JR;      /* J_r(psi_i) matrix, nr*ntheta-by-nrhs */
  gsl_matrix *JTHETA;  /* J_theta(psi_i) matrix, nr*ntheta-by-nrhs */
  gsl_matrix *JPHI;    /* J_phi(psi_i) matrix, nr*ntheta-by-nrhs */
  gsl_matrix *JPHI_HI; /* height-integrated J_phi(psi_i) matrix, ntheta-by-nrhs */
  gsl_matrix *ER;      /* E_r(psi_i) matrix, nr*ntheta-by-nrhs */
  gsl_matrix *ETHETA;  /* E_theta(psi_i) matrix, nr*ntheta-by-nrhs */
  gsl_matrix *EPHI;    /* E_phi(psi_i) matrix, nr*ntheta-by-nrhs */

  gsl_matrix *J_r;     /* J_r current solution (A/m^2) */
  gsl_matrix *J_theta; /* J_theta current solution (A/m^2) */
  gsl_matrix *J_phi;   /* J_phi current solution (A/m^2) */
  gsl_vector *J_lat;   /* height integrated current density profile (A/m) */

  gsl_vector *J_lat_E; /* current profile for u = 0 (A/m) */
  gsl_vector *J_lat_u; /* current profile for E = 0 (A/m) */

  gsl_matrix *E_r;     /* E_r electric field in V/m */
  gsl_matrix *E_theta; /* E_theta electric field in V/m */
  gsl_matrix *E_phi;   /* E_phi electric field in V/m */

  /* pde parameters */
  double *alpha;
  double *beta;
  double *gamma;

  /* pde coefficients, size nr*ntheta */
  double *f1;
  double *f2;
  double *f3;
  double *f4;
  double *f5;

  double E_phi0;     /* eastward electric field in V/m */
  int compute_winds; /* use winds in PDE solution? */

  /* sparse blas parameters */
  int nprocs;    /* number of processors to use */

  /* scaling factors for non-dimensionalization */
  double r_s;
  double sigma_s;
  double B_s;
  double U_s;
  double E_s;
  double psi_s;
  double J_s;

  size_t flags; /* PDE_FLG_xxx */

  int myid;   /* MPI id */

  hwm_workspace *hwm_workspace_p;
  mageq_workspace *mageq_workspace_p;
  msynth_workspace *msynth_workspace_p;
  sigma_workspace *sigma_workspace_p;
  slu_workspace *superlu_workspace_p;
  gsl_pde2d_workspace *gsl_pde2d_workspace_p;
} pde_workspace;

#define PDE_IDX(i, j, w)     ((i) * (w)->ntheta + (j))

/* index into 'psi' vector */
#define PSI_GET(p, i, j, w)  (gsl_vector_get(p, PDE_IDX((i), (j), (w))))

/*
 * Prototypes
 */

pde_workspace *pde_alloc(pde_parameters *params);
void pde_free(pde_workspace *w);
int pde_proc(const time_t t, const double longitude,
             pde_workspace *w);

#endif /* INCLUDED_pde_h */
