/*
 * window.h
 */

#ifndef INCLUDED_window_h
#define INCLUDED_window_h

/*
 * Prototypes
 */

int apply_ps1(const gsl_vector * in, gsl_vector * out);
int apply_hamming(const gsl_vector * in, gsl_vector * out);
int apply_kaiser(const gsl_vector * in, gsl_vector * out);
int apply_modsinsq(const gsl_vector * in, gsl_vector * out);

#endif /* INCLUDED_window_h */
