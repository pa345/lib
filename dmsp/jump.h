/*
 * jump.h
 */

#ifndef INCLUDED_jump_h
#define INCLUDED_jump_h

#include <gsl/gsl_filter.h>
#include <track/track.h>

#include "peak.h"

typedef struct
{
  double cumjump[3];             /* total sum of all jumps detected so far */
  gsl_vector *input[3];          /* input to first filter (residuals in VFM frame) */
  gsl_vector *output_impulse[3]; /* impulse-filtered first differences */
  gsl_vector_int *ioutlier[3];   /* outliers detected in Hampel filter */
  gsl_vector_int *ijump[3];      /* jump detected */
  size_t n;                      /* maximum length of time series */
  peak_workspace *peak_workspace_p;
  gsl_filter_impulse_workspace *impulse_workspace_p;
  gsl_filter_gaussian_workspace *gaussian_workspace_p;
} jump_workspace;

/*
 * Prototypes
 */

jump_workspace *jump_alloc(const size_t n, const size_t K);
void jump_free(jump_workspace * w);
int jump_proc(const int header, FILE *fp, const size_t start_idx, const size_t end_idx, size_t njump[3],
              satdata_mag * data, jump_workspace * w);
int jump_proc2(const int header, FILE *fp, size_t njump[3], track_data *tptr,
               satdata_mag * data, jump_workspace * w);

#endif /* INCLUDED_jump_h */
