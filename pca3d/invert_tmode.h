/*
 * invert_tmode.h
 */

#ifndef INCLUDED_invert_tmode_h
#define INCLUDED_invert_tmode_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

/* temporal modes data structure */
typedef struct
{
  size_t nfreq;          /* number of frequency bands */
  size_t *nmodes;        /* number of temporal modes in each frequency band, length nfreq */
  size_t N;              /* length of time series for each mode */

  /*
   * modes[i] is a matrix of size N-by-nmodes[i] containing all temporal
   * modes for frequency band i, one mode per column. N is the length of the
   * time series
   */
  gsl_matrix_complex ** modes;
} invert_tmode_workspace;

/*
 * Prototypes
 */

invert_tmode_workspace *invert_tmode_alloc(const size_t nfreq, const size_t nmodes[], const size_t N);
void invert_tmode_free(invert_tmode_workspace *w);
int invert_tmode_read_ascii(const char * filename, const size_t freq, const size_t mode, invert_tmode_workspace * w);
int invert_tmode_write_binary(const char * filename, invert_tmode_workspace * w);
invert_tmode_workspace * invert_tmode_read_binary(const char * filename);

#endif /* INCLUDED_invert_tmode_h */
