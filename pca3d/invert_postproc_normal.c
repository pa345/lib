/*
solve:

min_x || b - A x ||^2 + lambda^2 || L x ||^2
*/

static int
normal_reg_solution(const double lambda, const gsl_vector * L, const gsl_matrix * ATA,
                    const gsl_vector * ATb, const double bnorm, gsl_vector * x,
                    double * rnorm, double * snorm)
{
  const size_t p = ATA->size2;
  const double lambda_sq = lambda * lambda;
  const double bnorm_sq = bnorm * bnorm; /* ||b||^2 */
  gsl_matrix * work_ATA = gsl_matrix_alloc(p, p);
  gsl_vector * S = gsl_vector_alloc(p);
  gsl_vector * v = gsl_vector_alloc(p);
  double term1, term2;
  size_t i;

  /* work_ATA := ATA */
  gsl_matrix_tricpy(CblasLower, CblasNonUnit, work_ATA, ATA);

  /* work_ATA := ATA + lambda^2 L^T L */
  for (i = 0; i < p; ++i)
    {
      double Li = gsl_vector_get(L, i);
      double * ptr = gsl_matrix_ptr(work_ATA, i, i);

      *ptr += lambda_sq * Li * Li;
    }

  /* compute Cholesky decomposition */
  gsl_linalg_cholesky_decomp2(work_ATA, S);

  /* solve system */
  gsl_linalg_cholesky_solve2(work_ATA, S, ATb, x);

  /* compute solution norm ||L x|| */
  gsl_vector_memcpy(v, x);
  gsl_vector_mul(v, L);
  *snorm = gsl_blas_dnrm2(v);

  /* term1 := x^T ATA x */
  gsl_blas_dsymv(CblasLower, 1.0, ATA, x, 0.0, v);
  gsl_blas_ddot(x, v, &term1);

  /* term2 := x^T ATb */
  gsl_blas_ddot(ATb, x, &term2);

  /* compute residual norm ||b - Ax|| */
  *rnorm = GSL_MAX(bnorm_sq + term1 - 2.0 * term2, GSL_DBL_EPSILON);
  *rnorm = sqrt(*rnorm);

  fprintf(stderr, "lambda = %g, term1 = %g, term2 = %g, rsq = %g\n", lambda, term1, term2, term1 - 2.0*term2);

  gsl_matrix_free(work_ATA);
  gsl_vector_free(S);
  gsl_vector_free(v);

  return 0;
}
