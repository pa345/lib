/*
solve:

min_x || [ QTb ] - [    R     ] x ||^2
      || [  0  ]   [ lambda L ]   ||
*/

static int
tsqr_reg_solution(const double lambda, const gsl_vector * L, const gsl_matrix * R,
                  const gsl_vector * QTb, gsl_vector * x, double * rnorm, double * snorm)
{
  const size_t p = R->size2;
  gsl_matrix * U = gsl_matrix_alloc(p, p);
  gsl_matrix * Y = gsl_matrix_alloc(p, p);
  gsl_vector * D = gsl_vector_alloc(p);
  gsl_vector * b = gsl_vector_calloc(2 * p);
  gsl_matrix * T = gsl_matrix_alloc(p, p);
  gsl_vector * sol = gsl_vector_alloc(2 * p);
  gsl_vector * work = gsl_vector_alloc(p);
  gsl_vector_view v;

  gsl_matrix_tricpy(CblasUpper, CblasNonUnit, U, R);

  /* D = lambda*L */
  gsl_vector_axpby(lambda, L, 0.0, D);

  /* b = [ QTb ; 0 ] */
  v = gsl_vector_subvector(b, 0, p);
  gsl_vector_memcpy(&v.vector, QTb);

  gsl_linalg_QR_TD_decomp(U, D, Y, T);
  gsl_linalg_QR_TD_lssolve(U, Y, T, b, sol, work);

  v = gsl_vector_subvector(sol, 0, p);
  gsl_vector_memcpy(x, &v.vector);

  /* compute ||L x|| */
  gsl_vector_memcpy(work, x);
  gsl_vector_mul(work, L);
  *snorm = gsl_blas_dnrm2(work);

  /* compute rnorm = || QTb - R x || */
#if 1
  gsl_vector_memcpy(work, x);
  gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, R, work);
  gsl_vector_sub(work, QTb);
  *rnorm = gsl_blas_dnrm2(work);
#else
  v = gsl_vector_subvector(sol, p, p);
  *rnorm = gsl_blas_dnrm2(&v.vector);
#endif

  gsl_matrix_free(U);
  gsl_matrix_free(Y);
  gsl_matrix_free(T);
  gsl_vector_free(b);
  gsl_vector_free(sol);
  gsl_vector_free(D);
  gsl_vector_free(work);

  return 0;
}
