/*
 * mfield_fluxcal.h
 */

#ifndef INCLUDED_mfield_fluxcal_h
#define INCLUDED_mfield_fluxcal_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <satdata/satdata.h>

#include "mfield.h"

/* number of fit parameters */
#define FLUXCAL_P              9

#define FLUXCAL_IDX_SX         0
#define FLUXCAL_IDX_SY         1
#define FLUXCAL_IDX_SZ         2
#define FLUXCAL_IDX_OX         3
#define FLUXCAL_IDX_OY         4
#define FLUXCAL_IDX_OZ         5
#define FLUXCAL_IDX_U1         6
#define FLUXCAL_IDX_U2         7
#define FLUXCAL_IDX_U3         8

/*
 * Prototypes
 */

int mfield_fluxcal_print(const char *filename, const size_t sat_idx, const mfield_workspace *w);
int mfield_fluxcal_apply_datum(const gsl_vector *m, const double E[3], double B[4]);
int mfield_fluxcal_invapply_datum(const gsl_vector *m, const double B[3], double E[3]);
int mfield_fluxcal_jac(const gsl_vector *m, const double E[3], gsl_matrix *jac);

#endif /* INCLUDED_mfield_fluxcal_h */
