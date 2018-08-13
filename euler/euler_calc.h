/*
 * euler_calc.h
 */

#ifndef INCLUDED_euler_calc_h
#define INCLUDED_euler_calc_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

typedef struct
{
  size_t nmax;     /* maximum number of observations */
  size_t n;        /* number of observations added */
  size_t p;        /* number of Euler angles to compute */

  double *B_VFM;   /* vector measurements in VFM frame, size 3*nmax */
  double *B_model; /* vector model in NEC frame, size 3*nmax */
  double *q;       /* rotation quaternions from spacecraft-fixed to NEC, size 4*nmax */
  double *t;       /* timestamps (CDF_EPOCH), size nmax */
  double *qdlat;   /* QD latitudes (degrees), size nmax */
  size_t flags;    /* EULER_xxx flags */

  gsl_vector *c;   /* Euler angles, size p */
} euler_calc_workspace;

/*
 * Prototypes
 */

euler_calc_workspace * euler_calc_alloc(const size_t n);
void euler_calc_free(euler_calc_workspace * w);
int euler_calc_add(const double t, const double qdlat, const double B_VFM[3],
                   const double B_model[3], const double q[4], euler_calc_workspace * w);
int euler_calc_proc(euler_calc_workspace * w);
int euler_calc_print_residuals(const char * filename, euler_calc_workspace * w);

#endif /* INCLUDED_euler_calc_h */
