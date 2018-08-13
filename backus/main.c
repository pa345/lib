/*
 * main.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include <common/common.h>
#include <common/oct.h>

#include "green.h"
#include "lapack_wrapper.h"

/* generate random vector with elements in [-1,1] */
static void
random_vector(gsl_vector * v, gsl_rng * r)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      double vi = 2.0 * gsl_rng_uniform(r) - 1.0; /* in [-1,1] */
      gsl_vector_set(v, i, vi);
    }
}

int
random_orthogonal_matrix(gsl_matrix * A, gsl_rng * r)
{
  int s = 0;
  const size_t N = A->size1;
  gsl_vector *v = gsl_vector_alloc(N);
  size_t i;

  gsl_matrix_set_identity(A);

  for (i = 0; i < N; ++i)
    {
      double tau;

      random_vector(v, r);
      tau = gsl_linalg_householder_transform(v);
      gsl_linalg_householder_hm(tau, v, A);
    }

  gsl_vector_free(v);

  return s;
}

int
build_matrix(const size_t npts, gsl_matrix * A, green_workspace *green_p, gsl_rng *rng_p)
{
  const double r = R_EARTH_KM;
  size_t i;

  for (i = 0; i < npts; ++i)
    {
      gsl_vector_view X = gsl_matrix_row(A, 3*i);
      gsl_vector_view Y = gsl_matrix_row(A, 3*i + 1);
      gsl_vector_view Z = gsl_matrix_row(A, 3*i + 2);
      double theta = gsl_rng_uniform(rng_p) * M_PI;
      double phi = gsl_rng_uniform(rng_p) * 2.0 * M_PI;

      green_calc_int(r, theta, phi, X.vector.data, Y.vector.data, Z.vector.data, green_p);
    }

  return 0;
}

int
build_X(const gsl_vector * eval, const gsl_matrix * evec, gsl_matrix * X, gsl_rng * r)
{
  int s = 0;
  const size_t N = eval->size;
  gsl_matrix *Q = gsl_matrix_alloc(N, N);
  size_t i;

  random_orthogonal_matrix(Q, r);
  print_octave(Q, "Q");
  print_octave(evec, "P");
  printv_octave(eval, "eval");

  /* compute Q = Q D^{1/2} */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_column(Q, i);
      double di = gsl_vector_get(eval, i);

      gsl_vector_scale(&v.vector, sqrt(di));
    }

  /* compute Q = D^{-1/2} Q D^{1/2} */
  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(Q, i);
      double di = gsl_vector_get(eval, i);

      gsl_vector_scale(&v.vector, 1.0 / sqrt(di));
    }

  /* compute X = D^{-1/2} Q D^{1/2} P^T */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Q, evec, 0.0, X);

  gsl_matrix_memcpy(Q, X);

  /* compute X = P D^{-1/2} Q D^{1/2} P^T */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, evec, Q, 0.0, X);

  print_octave(X, "X");

  gsl_matrix_free(Q);

  return s;
}

int
main(int argc, char *argv[])
{
  const size_t nmax = 8;
  const size_t npts = 10;
  green_workspace *green_p = green_alloc(nmax, nmax, R_EARTH_KM);
  const size_t nnm = green_p->nnm;
  gsl_matrix *A = gsl_matrix_alloc(3 * npts, nnm);
  gsl_matrix *ATA = gsl_matrix_alloc(nnm, nnm);
  gsl_matrix *X = gsl_matrix_alloc(nnm, nnm);
  struct timeval tv0, tv1;
  gsl_vector *eval = gsl_vector_alloc(nnm);
  gsl_matrix *evec = gsl_matrix_alloc(nnm, nnm);
  gsl_vector *g = gsl_vector_alloc(nnm);       /* gauss coefficients */
  gsl_vector *Xg = gsl_vector_alloc(nnm);      /* X*g */
  gsl_vector *B1 = gsl_vector_alloc(3 * npts); /* A*g */
  gsl_vector *B2 = gsl_vector_alloc(3 * npts); /* A*X*g */
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  int eval_found;
  size_t i;

  fprintf(stderr, "main: building matrix (%zu-by-%zu)...", A->size1, A->size2);
  gettimeofday(&tv0, NULL);
  build_matrix(npts, A, green_p, r);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: computing A^T A...");
  gettimeofday(&tv0, NULL);
  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, A, 0.0, ATA);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: computing eigendecomposition of A^T A...");
  gettimeofday(&tv0, NULL);
  lapack_eigen_symmv(ATA, eval, evec, &eval_found);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: sorting eigenvalues...");
  gettimeofday(&tv0, NULL);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /*gsl_vector_fprintf(stdout, eval, "%.12e");*/

  fprintf(stderr, "main: building X matrix...");
  gettimeofday(&tv0, NULL);
  build_X(eval, evec, X, r);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  gsl_vector_set_zero(g);
  gsl_vector_set(g, green_nmidx(1, 0, green_p), 1.0);

  /* B1 = A * g */
  gsl_blas_dgemv(CblasNoTrans, 1.0, A, g, 0.0, B1);

  /* B2 = A * X * g */
  gsl_blas_dgemv(CblasNoTrans, 1.0, X, g, 0.0, Xg);
  gsl_blas_dgemv(CblasNoTrans, 1.0, A, Xg, 0.0, B2);

  fprintf(stderr, "main: || A g || = %.12e\n", gsl_blas_dnrm2(B1));
  fprintf(stderr, "main: || A X g || = %.12e\n", gsl_blas_dnrm2(B2));
  fprintf(stderr, "main: || A g || - || A X g || = %.12e\n", gsl_blas_dnrm2(B1) - gsl_blas_dnrm2(B2));

  printv_octave(B1, "B1");
  printv_octave(B2, "B2");

  print_octave(A, "A");
  printv_octave(g, "g");
  printv_octave(Xg, "Xg");

  for (i = 0; i < npts; ++i)
    {
      gsl_vector_view x1 = gsl_vector_subvector(B1, 3*i, 3);
      gsl_vector_view x2 = gsl_vector_subvector(B2, 3*i, 3);
      double F1 = gsl_blas_dnrm2(&x1.vector);
      double F2 = gsl_blas_dnrm2(&x2.vector);

      fprintf(stdout, "%.12e %.12e\n", F1, F2);
    }

  gsl_matrix_free(A);
  gsl_matrix_free(ATA);
  gsl_matrix_free(X);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
  gsl_vector_free(g);
  gsl_vector_free(Xg);
  gsl_vector_free(B1);
  gsl_vector_free(B2);
  green_free(green_p);
  gsl_rng_free(r);

  return 0;
}
