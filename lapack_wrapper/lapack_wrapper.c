#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <lapacke/lapacke.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include "lapack_wrapper.h"

/* minimize || A X - B ||
 *
 * Inputs: A    - least squares matrix
 *         B    - right hand sides
 *         X    - solution matrix
 *         rank - (output) effective rank of A
 */
int
lapack_lls(const gsl_matrix * A, const gsl_matrix * B, gsl_matrix * X, int *rank)
{
  int s;
  lapack_int m = A->size1;
  lapack_int n = A->size2;
  lapack_int nrhs = B->size2;
  lapack_int lda = A->size1;
  lapack_int ldb = B->size1;
  lapack_int lrank;
  lapack_int *jpvt = malloc(n * sizeof(lapack_int));
  gsl_matrix *work_A = gsl_matrix_alloc(A->size2, A->size1);
  gsl_matrix *work_B = gsl_matrix_alloc(B->size2, B->size1);
  double rcond = 1.0e-6;

  gsl_matrix_transpose_memcpy(work_A, A);
  gsl_matrix_transpose_memcpy(work_B, B);

  s = LAPACKE_dgelsy(LAPACK_COL_MAJOR,
                     m,
                     n,
                     nrhs,
                     work_A->data,
                     lda,
                     work_B->data,
                     ldb,
                     jpvt,
                     rcond,
                     &lrank);

  /* store solution in X */
  {
    gsl_matrix_view m = gsl_matrix_submatrix(work_B, 0, 0, X->size2, X->size1);
    gsl_matrix_transpose_memcpy(X, &m.matrix);
  }

#if 0
  /* 21 Jan 2017: computing all residual norms is way too slow */

  /* compute residual norms */
  {
    int i;
    gsl_vector_view r = gsl_matrix_row(work_B, 0);

    for (i = 0; i < nrhs; ++i)
      {
        gsl_vector_view xi = gsl_matrix_column(X, i);
        gsl_vector_const_view bi = gsl_matrix_const_column(B, i);
        double norm;

        gsl_vector_memcpy(&r.vector, &bi.vector);
        gsl_blas_dgemv(CblasNoTrans, -1.0, A, &xi.vector, 1.0, &r.vector);
        norm = gsl_blas_dnrm2(&r.vector);

        gsl_vector_set(rnorm, i, norm);
      }
  }
#endif

  *rank = lrank;

  gsl_matrix_free(work_A);
  gsl_matrix_free(work_B);
  free(jpvt);

  return s;
}

/* minimize || A x - b || */
int
lapack_lls2(const gsl_matrix * A, const gsl_vector * b, gsl_vector * x,
            int *rank)
{
  int s;
  lapack_int m = A->size1;
  lapack_int n = A->size2;
  lapack_int nrhs = 1;
  lapack_int lda = A->size1;
  lapack_int ldb = b->size;
  lapack_int lrank;
  lapack_int *jpvt = malloc(n * sizeof(lapack_int));
  gsl_matrix *work_A = gsl_matrix_alloc(A->size2, A->size1);
  gsl_vector *work_b = gsl_vector_alloc(b->size);
  double rcond = 1.0e-6;

  gsl_matrix_transpose_memcpy(work_A, A);
  gsl_vector_memcpy(work_b, b);

  s = LAPACKE_dgelsy(LAPACK_COL_MAJOR,
                     m,
                     n,
                     nrhs,
                     work_A->data,
                     lda,
                     work_b->data,
                     ldb,
                     jpvt,
                     rcond,
                     &lrank);

  /* store solution in X */
  {
    gsl_vector_view v = gsl_vector_subvector(work_b, 0, x->size);
    gsl_vector_memcpy(x, &v.vector);
  }

  *rank = lrank;

  gsl_matrix_free(work_A);
  gsl_vector_free(work_b);
  free(jpvt);

  return s;
}

int
lapack_complex_lls(const gsl_matrix_complex * A, const gsl_matrix_complex * B,
                   gsl_matrix_complex * X, int *rank)
{
  int s;
  lapack_int m = A->size1;
  lapack_int n = A->size2;
  lapack_int nrhs = B->size2;
  lapack_int lda = A->size1;
  lapack_int ldb = B->size1;
  lapack_int lrank;
  lapack_int *jpvt = malloc(n * sizeof(lapack_int));
  gsl_matrix_complex *work_A = gsl_matrix_complex_alloc(A->size2, A->size1);
  gsl_matrix_complex *work_B = gsl_matrix_complex_alloc(B->size2, B->size1);
  double rcond = 1.0e-6;

  gsl_matrix_complex_transpose_memcpy(work_A, A);
  gsl_matrix_complex_transpose_memcpy(work_B, B);

  s = LAPACKE_zgelsy(LAPACK_COL_MAJOR,
                     m,
                     n,
                     nrhs,
                     (lapack_complex_double *) work_A->data,
                     lda,
                     (lapack_complex_double *) work_B->data,
                     ldb,
                     jpvt,
                     rcond,
                     &lrank);

  /* store solution in X */
  {
    gsl_matrix_complex_view m = gsl_matrix_complex_submatrix(work_B, 0, 0, X->size2, X->size1);
    gsl_matrix_complex_transpose_memcpy(X, &m.matrix);
  }

  *rank = lrank;

  gsl_matrix_complex_free(work_A);
  gsl_matrix_complex_free(work_B);
  free(jpvt);

  return s;
}

/* symmetric eigenvalues using lower triangle */
int
lapack_eigen_symm(const gsl_matrix * m, gsl_vector *eval, int *eval_found)
{
  int s;
  const lapack_int N = m->size1;
  gsl_matrix *A = gsl_matrix_alloc(N, N);
  double vl = 0.0, vu = 0.0;
  lapack_int il = 0, iu = 0;
  lapack_int lda = A->size1;
  lapack_int ldz = A->size1;
  double abstol = 0.0;
  lapack_int M = 0;
  lapack_int *isuppz = malloc(2*N*sizeof(lapack_int));

  gsl_matrix_transpose_memcpy(A, m);

  s = LAPACKE_dsyevr(LAPACK_COL_MAJOR,
                     'N',
                     'A',
                     'L',
                     N,
                     A->data,
                     lda,
                     vl,
                     vu,
                     il,
                     iu,
                     abstol,
                     &M,
                     eval->data,
                     NULL,
                     ldz,
                     isuppz);

  /* sort into descending order */
  gsl_vector_reverse(eval);

  *eval_found = M;

  free(isuppz);
  gsl_matrix_free(A);

  return s;
}

/* symmetric eigenvalues and eigenvectors using lower triangle */
int
lapack_eigen_symmv(const gsl_matrix * m, gsl_vector *eval, gsl_matrix *evec,
                   int *eval_found)
{
  int s;
  const lapack_int N = m->size1;
  gsl_matrix *A = gsl_matrix_alloc(N, N);
  double vl = 0.0, vu = 0.0;
  lapack_int il = 0, iu = 0;
  lapack_int lda = A->size1;
  lapack_int ldz = evec->size1;
  double abstol = 0.0;
  lapack_int M = 0;
  lapack_int *isuppz = malloc(2*N*sizeof(lapack_int));

  gsl_matrix_transpose_memcpy(A, m);

  s = LAPACKE_dsyevr(LAPACK_COL_MAJOR,
                     'V',
                     'A',
                     'L',
                     N,
                     A->data,
                     lda,
                     vl,
                     vu,
                     il,
                     iu,
                     abstol,
                     &M,
                     eval->data,
                     evec->data,
                     ldz,
                     isuppz);

  gsl_matrix_transpose(evec);

  /* sort into descending order */
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  *eval_found = M;

  free(isuppz);
  gsl_matrix_free(A);

  return s;
}

int
lapack_eigen_herm(const gsl_matrix_complex * m, gsl_vector *eval, gsl_matrix_complex *evec,
                  int *eval_found)
{
  int s;
  const lapack_int N = m->size1;
  gsl_matrix_complex *A = gsl_matrix_complex_alloc(N, N);
  double vl = 0.0, vu = 0.0;
  lapack_int il = 0, iu = 0;
  lapack_int lda = A->size1;
  lapack_int ldz = evec->size1;
  double abstol = 0.0;
  lapack_int M = 0;
  lapack_int *isuppz = malloc(2*N*sizeof(lapack_int));

  gsl_matrix_complex_transpose_memcpy(A, m);

  s = LAPACKE_zheevr(LAPACK_COL_MAJOR,
                     'V',
                     'A',
                     'L',
                     N,
                     (lapack_complex_double *) A->data,
                     lda,
                     vl,
                     vu,
                     il,
                     iu,
                     abstol,
                     &M,
                     eval->data,
                     (lapack_complex_double *) evec->data,
                     ldz,
                     isuppz);

  gsl_matrix_complex_transpose(evec);

  /* sort into descending order */
  gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  *eval_found = M;

  free(isuppz);
  gsl_matrix_complex_free(A);

  return s;
}

/*
lapack_svd()
  Compute SVD of M-by-N matrix A

Inputs: A - M-by-N matrix
        S - (output) min(M,N) vector of singular values
        U - (output) M-by-M matrix of left singular vectors
        V - (output) N-by-N matrix of right singular vectors

Return: success/error
*/

int
lapack_svd(const gsl_matrix * A, gsl_vector * S, gsl_matrix * U, gsl_matrix * V)
{
  int s;
  lapack_int M = A->size1;
  lapack_int N = A->size2;
  lapack_int lda = A->size1;
  lapack_int ldu = U->size1;
  lapack_int ldvt = V->size1;
  gsl_matrix *work_A = gsl_matrix_alloc(A->size2, A->size1);

  gsl_matrix_transpose_memcpy(work_A, A);

  s = LAPACKE_dgesdd(LAPACK_COL_MAJOR,
                     'A',
                     M,
                     N,
                     work_A->data,
                     lda,
                     S->data,
                     U->data,
                     ldu,
                     V->data,
                     ldvt);

  /* no need to transpose V since V^T is computed already by
   * lapack routine */
  gsl_matrix_transpose(U);

  gsl_matrix_free(work_A);

  return s;
}

/*
lapack_complex_svd()
  Compute SVD of complex M-by-N matrix A

Inputs: A - M-by-N complex matrix
        S - (output) min(M,N) real vector of singular values
        U - (output) M-by-M complex matrix of left singular vectors
        V - (output) N-by-N complex matrix of right singular vectors

Return: success/error
*/

int
lapack_complex_svd(const gsl_matrix_complex * A, gsl_vector * S,
                   gsl_matrix_complex * U, gsl_matrix_complex * V)
{
  size_t M = A->size1;
  size_t N = A->size2;

  if (U->size1 != U->size2)
    {
      GSL_ERROR ("U matrix must be square", GSL_ENOTSQR);
    }
  else if (U->size1 != M)
    {
      GSL_ERROR ("U matrix must be M-by-M", GSL_EBADLEN);
    }
  else if (S->size != GSL_MIN(M, N))
    {
      GSL_ERROR ("S must have length MIN(M,N)", GSL_EBADLEN);
    }
  else if (V->size1 != V->size2)
    {
      GSL_ERROR ("V matrix must be square", GSL_ENOTSQR);
    }
  else if (V->size1 != N)
    {
      GSL_ERROR ("V matrix must be N-by-N", GSL_EBADLEN);
    }
  else
    {
      int s;
      lapack_int lda = A->size1;
      lapack_int ldu = U->size1;
      lapack_int ldvt = V->size1;
      gsl_matrix_complex *work_A = gsl_matrix_complex_alloc(N, M);

      gsl_matrix_complex_transpose_memcpy(work_A, A);

      s = LAPACKE_zgesdd(LAPACK_COL_MAJOR,
                         'A',
                         (lapack_int) M,
                         (lapack_int) N,
                         (lapack_complex_double *) work_A->data,
                         lda,
                         S->data,
                         (lapack_complex_double *) U->data,
                         ldu,
                         (lapack_complex_double *) V->data,
                         ldvt);

      /* no need to transpose V since V^H is computed already by
       * lapack routine */
      gsl_matrix_complex_transpose(U);

      gsl_matrix_complex_free(work_A);

      return s;
    }
}

/*
lapack_complex_svd_thin()
  Compute thin SVD of complex M-by-N matrix A

Inputs: A - M-by-N complex matrix
        S - (output) min(M,N) real vector of singular values
        U - (output) M-by-min(M,N) complex matrix of left singular vectors
        V - (output) min(M,N)-by-N complex matrix of right singular vectors

Return: success/error
*/

int
lapack_complex_svd_thin(const gsl_matrix_complex * A, gsl_vector * S,
                        gsl_matrix_complex * U, gsl_matrix_complex * V)
{
  const size_t M = A->size1;
  const size_t N = A->size2;
  const size_t minMN = GSL_MIN(M, N);

  if ((U->size1 != M) || (U->size2 != minMN))
    {
      GSL_ERROR ("U matrix must be M-by-min(M,N)", GSL_EBADLEN);
    }
  else if (S->size != minMN)
    {
      GSL_ERROR ("S must have length min(M,N)", GSL_EBADLEN);
    }
  else if ((V->size1 != minMN) || (V->size2 != N))
    {
      GSL_ERROR ("V matrix must be min(M,N)-by-N", GSL_EBADLEN);
    }
  else
    {
      int s;
      lapack_int lda = A->size1;
      lapack_int ldu = U->size1;
      lapack_int ldvt = V->size1;
      gsl_matrix_complex *work_A = gsl_matrix_complex_alloc(N, M);
      gsl_matrix_complex *work_U = gsl_matrix_complex_alloc(minMN, M);

      gsl_matrix_complex_transpose_memcpy(work_A, A);

      s = LAPACKE_zgesdd(LAPACK_COL_MAJOR,
                         'S',
                         (lapack_int) M,
                         (lapack_int) N,
                         (lapack_complex_double *) work_A->data,
                         lda,
                         S->data,
                         (lapack_complex_double *) work_U->data,
                         ldu,
                         (lapack_complex_double *) V->data,
                         ldvt);

      /* no need to transpose V since V^H is computed already by
       * lapack routine */
      gsl_matrix_complex_transpose_memcpy(U, work_U);

      gsl_matrix_complex_free(work_A);
      gsl_matrix_complex_free(work_U);

      return s;
    }
}

/*
lapack_cholesky_solve()
  solve A x = b using Cholesky factorization of A

Inputs: A     - matrix A
        b     - rhs vector b
        x     - (output) solution vector
        rcond - (output) reciprocal condition number
        L     - (output) Cholesky factor of A in lower triangle;
                can be set to NULL
*/

int
lapack_cholesky_solve(const gsl_matrix * A, const gsl_vector * b, gsl_vector * x,
                      double * rcond, gsl_matrix * L)
{
  int s;
  lapack_int N = A->size1;
  lapack_int nrhs = 1;
  lapack_int lda = A->size1;
  lapack_int ldb = b->size;
  lapack_int ldx = x->size;
  gsl_vector *S = gsl_vector_alloc(N);
  gsl_vector *work_b = gsl_vector_alloc(N);
  gsl_matrix *work_A = gsl_matrix_alloc(N, N);
  gsl_matrix *AF = gsl_matrix_alloc(N, N);
  lapack_int ldaf = AF->size1;
  char fact = (L == NULL) ? 'E' : 'N';
  char equed;
  double ferr, berr;

  gsl_vector_memcpy(work_b, b);
  gsl_matrix_transpose_memcpy(work_A, A);

  s = LAPACKE_dposvx(LAPACK_COL_MAJOR,
                     fact,
                     'L',
                     N,
                     nrhs,
                     work_A->data,
                     lda,
                     AF->data,
                     ldaf,
                     &equed,
                     S->data,
                     work_b->data,
                     ldb,
                     x->data,
                     ldx,
                     rcond,
                     &ferr,
                     &berr);

  if (L != NULL)
    {
      /* store Cholesky factor in L */
      gsl_matrix_transpose_memcpy(L, AF);
    }

  gsl_vector_free(S);
  gsl_vector_free(work_b);
  gsl_matrix_free(AF);
  gsl_matrix_free(work_A);

  return s;
}

/*
lapack_cholesky_decomp()
  Decompose A using Cholesky factorization

Inputs: A - on input, the SPD matrix in lower triangle
            on output, the Cholesky factor in lower triangle
*/

int
lapack_cholesky_decomp(gsl_matrix * A)
{
  int s;
  lapack_int N = A->size1;
  lapack_int lda = A->size1;
  gsl_matrix *work_A = gsl_matrix_alloc(N, N);

  gsl_matrix_transpose_tricpy('L', 1, work_A, A);

  s = LAPACKE_dpotrf(LAPACK_COL_MAJOR,
                     'L',
                     N,
                     work_A->data,
                     lda);

  gsl_matrix_transpose_tricpy('U', 1, A, work_A);

  gsl_matrix_free(work_A);

  return s;
}

/*
lapack_cholesky_invert()
  Invert A using Cholesky factorization, which was previously
computed with lapack_cholesky_solve() (the L factor output)

Inputs: A - on input, the Cholesky factor L in lower triangle
            on output, (L L^T)^{-1} in lower triangle
*/

int
lapack_cholesky_invert(gsl_matrix * A)
{
  int s;
  lapack_int N = A->size1;
  lapack_int lda = A->size1;
  gsl_matrix *work_A = gsl_matrix_alloc(N, N);

  gsl_matrix_transpose_tricpy('L', 1, work_A, A);

  s = LAPACKE_dpotri(LAPACK_COL_MAJOR,
                     'L',
                     N,
                     work_A->data,
                     lda);

  gsl_matrix_transpose_tricpy('U', 1, A, work_A);

  gsl_matrix_free(work_A);

  return s;
}

/* solve A x = b for square Hermitian system and compute condition number */
int
lapack_zposv(const gsl_matrix_complex * A, const gsl_vector_complex * b, gsl_vector_complex *x,
             double * rcond)
{
  int s = 0;
  lapack_int n = A->size1;
  lapack_int nrhs = 1;
  lapack_int lda = A->size1;
  lapack_int ldb = b->size;
  lapack_int ldx = x->size;
  gsl_matrix_complex *work_A = gsl_matrix_complex_alloc(n, n);
  gsl_vector_complex *work_b = gsl_vector_complex_alloc(n);
  gsl_matrix_complex *AF = gsl_matrix_complex_alloc(n, n);
  gsl_vector *S = gsl_vector_alloc(n);
  double ferr, berr;
  lapack_int ldaf = n;
  char equed = 'N';

  gsl_matrix_complex_transpose_memcpy(work_A, A);
  gsl_vector_complex_memcpy(work_b, b);

  /* use expert driver to get condition number estimate */
  s = LAPACKE_zposvx(LAPACK_COL_MAJOR,
                     'E',
                     'U',
                     n,
                     nrhs,
                     (lapack_complex_double *) work_A->data,
                     lda,
                     (lapack_complex_double *) AF->data,
                     ldaf,
                     &equed,
                     S->data,
                     (lapack_complex_double *) work_b->data,
                     ldb,
                     (lapack_complex_double *) x->data,
                     ldx,
                     rcond,
                     &ferr,
                     &berr);

  gsl_matrix_complex_free(work_A);
  gsl_matrix_complex_free(AF);
  gsl_vector_complex_free(work_b);
  gsl_vector_free(S);

  return s;
}
