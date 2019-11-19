/*
 * mfield_error.h
 */

#ifndef INCLUDED_mfield_error_h
#define INCLUDED_mfield_error_h

/*
 * Prototypes
 */

int mfield_print_uncertainties(const char * filename, const double t, const gsl_matrix * covar, mfield_workspace * w);
int mfield_covariance(gsl_matrix * covar, mfield_workspace *w);
int mfield_correlation(gsl_matrix * C);

#endif /* INCLUDED_mfield_error_h */
