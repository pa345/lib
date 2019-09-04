/*
 * mfield_sht_test.c
 *
 * This program tests the NFFT method for computing spherical
 * harmonic expansions of scattered points. It generates
 * random (theta,phi) points and computes the vector field components
 * using the NFFT and direct methods for comparison.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include <mainlib/ml_common.h>

#include "mfield_eval.h"
#include "mfield_sht.h"

int
create_random_vector(gsl_vector *v, const double lower,
                     const double upper, gsl_rng *r)
{
  const size_t n = v->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double x = gsl_rng_uniform(r); /* in [0,1] */
      gsl_vector_set(v, i, (upper - lower) * x + lower);
    }

  return 0;
} /* create_random_vector() */

int
compute_nfft(const gsl_vector *g, const gsl_vector *r,
             const gsl_vector *theta, const gsl_vector *phi,
             gsl_vector *Y, gsl_vector *Z, mfield_sht_workspace *w)
{
  int s = 0;
  double r0 = gsl_vector_get(r, 0);

  s += mfield_sht_calc(g, r0, theta, phi, Y, Z, w);

  return s;
}

int
compute_direct(const gsl_vector *r, const gsl_vector *theta,
               const gsl_vector *phi, gsl_vector *Y, gsl_vector *Z,
               mfield_eval_workspace *w)
{
  const size_t M = theta->size;
  size_t i;
  const double t = w->epoch;

  for (i = 0; i < M; ++i)
    {
      double ri = gsl_vector_get(r, i) * w->R;
      double ti = gsl_vector_get(theta, i);
      double pi = gsl_vector_get(phi, i);
      double B[4];

      mfield_eval(t, ri, ti, pi, B, w);

      /* store output vector components */
      gsl_vector_set(Y, i, B[1]);
      gsl_vector_set(Z, i, B[2]);
    }

  return 0;
} /* compute_direct() */

int
main(int argc, char *argv[])
{
  const size_t M = 100000; /* number of point nodes */
  mfield_sht_workspace *sht_p;
  mfield_eval_workspace *eval_p;

#if 0
  eval_p = mfield_eval_read("dmsp_mag_1.txt");
#else
  eval_p = mfield_eval_read("dmsp_mag_MF7.txt");
#endif
  if (!eval_p)
    exit(1);

  sht_p = mfield_sht_alloc(eval_p->nmax, M);
  if (!sht_p)
    exit(1);

  fprintf(stderr, "main: number of scattered node points: %zu\n", M);
  fprintf(stderr, "main: using model to degree %zu\n", eval_p->nmax);

  {
    const size_t nnm = eval_p->nnm;
    size_t i;
    gsl_vector *r = gsl_vector_alloc(M);
    gsl_vector *theta = gsl_vector_alloc(M);
    gsl_vector *phi = gsl_vector_alloc(M);
    gsl_vector *Y = gsl_vector_alloc(M);
    gsl_vector *Y_direct = gsl_vector_alloc(M);
    gsl_vector *Z = gsl_vector_alloc(M);
    gsl_vector *Z_direct = gsl_vector_alloc(M);
    gsl_vector_view gv = gsl_vector_view_array(eval_p->c, nnm);
    gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);
    double max_yerr = 0.0, max_zerr = 0.0;
    struct timeval tv0, tv1;

#if 0
    gsl_vector_set_zero(&gv.vector);
    gsl_vector_set(&gv.vector, mfield_eval_nmidx(1,-1), 100.43);
    gsl_vector_set(&gv.vector, mfield_eval_nmidx(1,1), -388.56);
    gsl_vector_set(&gv.vector, mfield_eval_nmidx(1,0), 2000.38);
#endif

    create_random_vector(theta, 0.0, M_PI, rng_p);
    create_random_vector(phi, 0.0, 2.0 * M_PI, rng_p);

    /*
     * currently only fixed radii points are allowed, due to
     * limitations in NFFT; this r vector is scaled so that
     * r~ = r / a
     */
    gsl_vector_set_all(r, 0.5);

    /* compute Y,Z components with NFFT */
    fprintf(stderr, "Computing via FFT method...");
    gettimeofday(&tv0, NULL);
    compute_nfft(&gv.vector, r, theta, phi, Y, Z, sht_p);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

    fprintf(stderr, "Computing via direct method...");
    gettimeofday(&tv0, NULL);
    compute_direct(r, theta, phi, Y_direct, Z_direct, eval_p);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

    for (i = 0; i < M; ++i)
      {
        double zdi = gsl_vector_get(Z_direct, i);
        double zi = gsl_vector_get(Z, i);
        double zerr = (zi - zdi) / zdi;
        double ydi = gsl_vector_get(Y_direct, i);
        double yi = gsl_vector_get(Y, i);
        double yerr = (yi - ydi) / ydi;

        /*printf("%.12e %.12e %.12e\n", zdi, zi, zerr);*/

        if (fabs(yerr) > max_yerr)
          max_yerr = fabs(yerr);
        if (fabs(zerr) > max_zerr)
          max_zerr = fabs(zerr);
      }

    fprintf(stderr, "Y max relative error = %.12e\n", max_yerr);
    fprintf(stderr, "Z max relative error = %.12e\n", max_zerr);

    gsl_vector_free(r);
    gsl_vector_free(theta);
    gsl_vector_free(phi);
    gsl_vector_free(Y);
    gsl_vector_free(Y_direct);
    gsl_vector_free(Z);
    gsl_vector_free(Z_direct);
    gsl_rng_free(rng_p);
  }

  mfield_eval_free(eval_p);
  mfield_sht_free(sht_p);

  return 0;
}
