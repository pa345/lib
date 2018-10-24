/*
 * euler.h
 */

#ifndef INCLUDED_euler_h
#define INCLUDED_euler_h

#include <satdata/satdata.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "magdata.h"

#define EULER_FLG_ZYX               (1 << 0) /* ZYX Euler convention (Swarm) */
#define EULER_FLG_ZYZ               (1 << 1) /* ZYZ Euler convention (Swarm ASM-V) */
#define EULER_FLG_DERIV_ALPHA       (1 << 2) /* first derivative wrt alpha */
#define EULER_FLG_DERIV_BETA        (1 << 3) /* first derivative wrt beta */
#define EULER_FLG_DERIV_GAMMA       (1 << 4) /* first derivative wrt gamma */
#define EULER_FLG_DERIV2_ALPHA      (1 << 5) /* second derivative wrt alpha */
#define EULER_FLG_DERIV2_BETA       (1 << 6) /* second derivative wrt beta */
#define EULER_FLG_DERIV2_GAMMA      (1 << 7) /* second derivative wrt gamma */
#define EULER_FLG_ROTSC             (1 << 8) /* rotate to NEC with star camera quaternions */
#define EULER_FLG_RINV              (1 << 9) /* apply R_inv matrix prior to Euler rotation (Swarm ASM-V) */

#define EULER_IDX_ALPHA             0
#define EULER_IDX_BETA              1
#define EULER_IDX_GAMMA             2

/* maximum number of bins for Euler angles */
#define EULER_MAX_BINS             1000

#define EULER_MAX_BUFFER           2048

typedef struct
{
  double t[EULER_MAX_BINS];     /* timestamp (CDF_EPOCH) sorted low to high */
  double alpha[EULER_MAX_BINS]; /* alpha (radians) */
  double beta[EULER_MAX_BINS];  /* beta (radians) */
  double gamma[EULER_MAX_BINS]; /* gamma (radians) */
  size_t n;                     /* number of data */
  size_t flags;                 /* EULER_FLG_xxx for Euler convention */
} euler_workspace;

/*
 * Prototypes
 */

euler_workspace *euler_alloc(const size_t flags);
void euler_free(euler_workspace *w);
euler_workspace *euler_read(const char *filename);
int euler_write(const char *filename, const euler_workspace *w);
int euler_write_swarm(const double fday_start, const double fday_end,
                      const double fday_step, const char *filename,
                      const euler_workspace *w);
int euler_add(const double t, const gsl_vector *x, euler_workspace *w);
int euler_apply(satdata_mag *data, const euler_workspace *w);
int euler_magdata_apply(magdata *data, const euler_workspace *w);
int euler_nec2vfm_t(const double t, const double q[],
                    const double B_in[3], double B_out[3],
                    const euler_workspace *w);
int euler_vfm2nec_t(const double t, const double q[],
                    const double B_in[3], double B_out[3],
                    const euler_workspace *w);
int euler_vfm2nec(const size_t flags, const double alpha, const double beta,
                  const double gamma, const double q[],
                  const double B_in[3], double B_out[3]);
int euler_matrix_vfm2nec(const size_t flags, const double alpha, const double beta,
                         const double gamma, const double q[],
                         const gsl_matrix *B_in, gsl_matrix *B_out);
int euler_nec2vfm(const size_t flags, const double alpha, const double beta,
                  const double gamma, const double q[],
                  const double B_in[3], double B_out[3]);
int euler_apply_R3(const size_t flags, const double alpha, const double beta,
                   const double gamma, const double B_in[3], double B_out[3]);
int euler_matrix_apply_R3(const size_t flags, const double alpha, const double beta,
                          const double gamma, const gsl_matrix *B_in, gsl_matrix *B_out);
int euler_Rq(const double *q, gsl_matrix *Rq);

#endif /* INCLUDED_euler_h */
