/*
 * mfield_align.h
 */

#ifndef INCLUDED_mfield_align_h
#define INCLUDED_mfield_align_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <mainlib/ml_att.h>

#include "mfield.h"

/* number of fit parameters */
#define ALIGN_P              3

/*
 * Prototypes
 */

int mfield_align_print(const char *filename, const size_t sat_idx, const mfield_workspace *w);
int mfield_align_nec2vfm(const gsl_vector * align_params, const double * q, const double B_nec[3], double B_vfm[3],
                         const att_workspace * att_p, const mfield_workspace * w);
int mfield_align_vfm2nec(const gsl_vector * align_params, const double * q, const double B_vfm[3], double B_nec[3],
                         const att_workspace * att_p, const mfield_workspace * w);
int mfield_align_matrix_vfm2nec(const gsl_vector * align_params, const double * q, const gsl_matrix * in, gsl_matrix * out, const att_workspace * att_p,
                                const mfield_workspace * w);
int mfield_align_deriv_vfm2nec(const size_t idx, const gsl_vector * align_params, const double * q, const double B_vfm[3], double B_nec[3],
                               const att_workspace * att_p, const mfield_workspace * w);

#endif /* INCLUDED_mfield_align_h */
