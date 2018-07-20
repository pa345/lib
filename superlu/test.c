/*
 * test.c
 * Patrick Alken
 *
 * Test superlu module by creating random sparse matrices and rhs
 * vectors, solving the systems and comparing with GSL output
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spmatrix.h>

#include "superlu.h"

void
test_vectors(const double tol, const gsl_vector *x, const gsl_vector *x_exact)
{
  size_t i;

  for (i = 0; i < x->size; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double yi = gsl_vector_get(x_exact, i);
      gsl_test_rel(xi, yi, tol, "n = %zu, i = %d", x->size, i);
    }
}

/*
create_random_sparse()
  Create a random sparse matrix with approximately
M*N*density non-zero entries

Inputs: M       - number of rows
        N       - number of columns
        density - sparse density \in [0,1]
                  0 = no non-zero entries
                  1 = all m*n entries are filled
        r       - random number generator

Return: pointer to sparse matrix in triplet format (must be freed by caller)

Notes:
1) non-zero matrix entries are uniformly distributed in [0,1]
*/

static gsl_spmatrix *
create_random_sparse(const size_t M, const size_t N, const double density,
                     const gsl_rng *r)
{
  size_t nnzwanted = (size_t) floor(M * N * GSL_MIN(density, 1.0));
  gsl_spmatrix *m = gsl_spmatrix_alloc_nzmax(M, N,
                                             nnzwanted,
                                             GSL_SPMATRIX_TRIPLET);
  size_t i, j;

  /* set all diagonal elements */
  for (i = 0; i < M; ++i)
    {
      double x = gsl_rng_uniform(r);
      gsl_spmatrix_set(m, i, i, x);
    }

  while (gsl_spmatrix_nnz(m) < nnzwanted)
    {
      double x = gsl_rng_uniform(r);

      /* generate a random row and column */
      i = gsl_rng_uniform(r) * M;
      j = gsl_rng_uniform(r) * N;

      gsl_spmatrix_set(m, i, j, x);
    }

  return m;
}

void
create_random_sparse_matrix(gsl_matrix *m, gsl_rng *r, double lower,
                            double upper)
{
  size_t N = m->size1;
  size_t M = m->size2;
  size_t i, j;
  int numel; /* number of non-zero elements in a row */
  double x;

  gsl_matrix_set_zero(m);

  if (N == 1)
    {
      x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_matrix_set(m, 0, 0, x);
      return;
    }

  for (i = 0; i < N; ++i)
    {
      /* pick a random number between 1 and M/2 - this is how many
       * nonzero elements are in this row
       */
      numel = (int) (gsl_rng_uniform(r) * (M / 2 - 1) + 1);
      for (j = 0; j < (size_t) numel; ++j)
        {
          int k = (int) (gsl_rng_uniform(r) * (M - 2));
          x = gsl_rng_uniform(r) * (upper - lower) + lower;
          gsl_matrix_set(m, i, k, x);
        }

      /* always set the diagonal element */
      x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_matrix_set(m, i, i, x);
    }
} /* create_random_sparse_matrix() */

void
create_random_vector(gsl_vector *v, gsl_rng *r, double lower, double upper)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      double x = gsl_rng_uniform(r) * (upper - lower) + lower;
      gsl_vector_set(v, i, x);
    }
}

void
test_superlu(const size_t nrhs)
{
  const int N_max = 50;
  const int nprocs = 1;
  int n, i;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  for (n = nrhs; n <= N_max; ++n)
    {
      gsl_matrix *A = gsl_matrix_alloc(n, n);
      gsl_matrix *B = gsl_matrix_alloc(n, nrhs);
      gsl_matrix *X = gsl_matrix_alloc(n, nrhs);
      gsl_matrix *T = gsl_matrix_alloc(n, nrhs);
      slu_workspace *w = slu_alloc(n, n, nprocs, nrhs);

      for (i = 0; i < 20; ++i)
        {
          gsl_spmatrix *S = create_random_sparse(n, n, 0.1, r);
          gsl_spmatrix *C = gsl_spmatrix_ccs(S);
          size_t j;

          for (j = 0; j < nrhs; ++j)
            {
              gsl_vector_view bj = gsl_matrix_column(B, j);
              create_random_vector(&bj.vector, r, -10.0, 10.0);
            }

          slu_proc(C, B, X, w);

          /* form T = A X */
          gsl_spmatrix_sp2d(A, S);
          gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, X, 0.0, T);

          for (j = 0; j < nrhs; ++j)
            {
              gsl_vector_view bj = gsl_matrix_column(B, j);
              gsl_vector_view tj = gsl_matrix_column(T, j);
              test_vectors(1.0e-7, &tj.vector, &bj.vector);
            }

          /* T = A X - B */
          gsl_matrix_sub(T, B);

          for (j = 0; j < nrhs; ++j)
            {
              gsl_vector_view tj = gsl_matrix_column(T, j);
              double rnorm = gsl_blas_dnrm2(&tj.vector);

              gsl_test_abs(w->rnorm[j], rnorm, 1.0e-8, "residual n = %zu, i = %zu cond = %.12e",
                           n, i, 1.0 / w->rcond);

              gsl_test_abs(w->rnorm[j], 0.0, 1.0e-8, "residual n = %zu, i = %zu cond = %.12e",
                           n, i, 1.0 / w->rcond);
            }

          gsl_spmatrix_free(S);
          gsl_spmatrix_free(C);
        }

      gsl_matrix_free(A);
      gsl_matrix_free(B);
      gsl_matrix_free(X);
      gsl_matrix_free(T);
      slu_free(w);
    }

  gsl_rng_free(r);
}

int
main()
{
  slu_workspace *sw;
  gsl_spmatrix *A, *C;
  size_t i;
  const int m = 5;
  const int n = 5;
  double rhs[m];
  double sol[n];
  gsl_matrix_view B = gsl_matrix_view_array(rhs, m, 1);
  gsl_matrix_view X = gsl_matrix_view_array(sol, n, 1);
  int nprocs;
  double s, u, p, e, r, l;
  double min, max;
  double x[] = { -0.0312500000000000, 0.0654761904761905, 0.0133928571428571,
                  0.0625000000000000, 0.0327380952380952 };

  nprocs = 1;
  s = 19.0;
  u = 21.0;
  p = 16.0;
  e = 5.0;
  r = 18.0;
  l = 12.0;

  for (i = 0; i < (size_t) n; ++i)
    rhs[i] = 1.0;

  sw = slu_alloc(m, n, nprocs, 1);
  A = gsl_spmatrix_alloc(m, n);

  gsl_spmatrix_set(A, 0, 0, s);
  gsl_spmatrix_set(A, 1, 0, l);
  gsl_spmatrix_set(A, 4, 0, l);
  gsl_spmatrix_set(A, 1, 1, u);
  gsl_spmatrix_set(A, 2, 1, l);
  gsl_spmatrix_set(A, 4, 1, l);
  gsl_spmatrix_set(A, 0, 2, u);
  gsl_spmatrix_set(A, 2, 2, p);
  gsl_spmatrix_set(A, 0, 3, u);
  gsl_spmatrix_set(A, 3, 3, e);
  gsl_spmatrix_set(A, 3, 4, u);
  gsl_spmatrix_set(A, 4, 4, r);

  gsl_spmatrix_minmax(A, &min, &max);
  fprintf(stderr, "min = %f, max = %f\n", min, max);

  C = gsl_spmatrix_ccs(A);

  slu_proc(C, &B.matrix, &X.matrix, sw);

  {
    gsl_vector_view xv = gsl_vector_view_array(x, n);
    gsl_vector_view solv = gsl_matrix_column(&X.matrix, 0);
    test_vectors(1.0e-8, &solv.vector, &xv.vector);
  }

  printf("sol = [\n");
  for (i = 0; i < (size_t) n; ++i)
    printf("%.12e\n", sol[i]);
  printf("]\n");

  printf("residual = %.12e\n", slu_residual(sw));
  printf("cond(A)  = %.12e\n", 1.0 / sw->rcond);

  test_superlu(1);
  test_superlu(2);
  test_superlu(3);
  test_superlu(10);

  slu_free(sw);
  gsl_spmatrix_free(A);
  gsl_spmatrix_free(C);

  exit (gsl_test_summary());

  return 0;
}
