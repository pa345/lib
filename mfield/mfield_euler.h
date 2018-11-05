/*
 * mfield_euler.h
 */

#ifndef INCLUDED_mfield_euler_h
#define INCLUDED_mfield_euler_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <satdata/satdata.h>

#include "mfield.h"

/* number of fit parameters */
#define EULER_P              3

/*
 * Prototypes
 */

int mfield_euler_print(const char *filename, const size_t sat_idx, const mfield_workspace *w);

#endif /* INCLUDED_mfield_euler_h */
