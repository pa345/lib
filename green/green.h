/*
 * green.h
 */

#ifndef INCLUDED_green_h
#define INCLUDED_green_h

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

/* mu_0 in units of: nT / (kA km^{-1}) */
#define GREEN_MU_0                  (400.0 * M_PI)

typedef struct
{
  size_t nmax;     /* maximum spherical harmonic degree */
  size_t mmax;     /* maximum spherical harmonic order */
  size_t nnm;      /* number of total Green's functions for nmax */
  double R;        /* reference radius (km) */
  double *cosmphi; /* array of cos(m phi) values, size mmax + 1 */
  double *sinmphi; /* array of sin(m phi) values, size mmax + 1 */
  double *Plm;     /* associated Legendre functions */
  double *dPlm;    /* derivatives of associated Legendre functions */
  double *work;    /* workspace, size nnm */
} green_workspace;

/*
 * Prototypes
 */

green_workspace *green_alloc(const size_t nmax, const size_t mmax, const double R);
void green_free(green_workspace *w);
int green_calc_int(const double r, const double theta, const double phi,
                   double *X, double *Y, double *Z, green_workspace *w);
int green_calc_int2(const double r, const double theta, const double phi,
                    gsl_vector *X, gsl_vector *Y, gsl_vector *Z, green_workspace *w);
int green_calc_ext(const double r, const double theta, const double phi,
                   double *X, double *Y, double *Z, green_workspace *w);
int green_potential_calc_ext(const double r, const double theta, const double phi,
                             double *V, green_workspace *w);
int green_Y_calc(const double theta, const double phi,
                 double *Y, green_workspace *w);
extern inline size_t green_nmidx(const size_t n, const int m, const green_workspace *w);
extern inline size_t green_idx(const size_t n, const int m, const size_t mmax);
size_t green_nnm(const green_workspace *w);
size_t green_calc_nnm(const size_t nmax, const size_t mmax);
int green_print_spectrum(const char *filename, const gsl_vector *c,
                         const green_workspace *w);
int green_print_spectrum_azim(const char *filename, const gsl_vector * c,
                              const green_workspace *w);
int green_k2g(const double b, const gsl_vector *k, gsl_vector *g,
              const green_workspace *w);
int green_g2k(const double b, const gsl_vector *g, gsl_vector *k,
              const green_workspace *w);

/* current.c */
int green_eval_sheet_int(const double b, const double theta, const double phi,
                         const gsl_vector *g, double K[3], green_workspace *w);
double green_eval_chi_int(const double b, const double theta, const double phi,
                          const gsl_vector *g, green_workspace *w);
double green_eval_chi_ext(const double b, const double theta, const double phi,
                          const gsl_vector *k, green_workspace *w);

#endif /* INCLUDED_green_h */
