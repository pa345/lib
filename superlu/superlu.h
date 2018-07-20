/*
 * superlu.h
 * Patrick Alken
 */

#ifndef INCLUDED_superlu_h
#define INCLUDED_superlu_h

#include <gsl/gsl_spmatrix.h>

#include <superlu_mt/slu_mt_ddefs.h>

typedef struct
{
  size_t size1;  /* number of rows in A */
  size_t size2;  /* number of columns in A */

  SuperMatrix A;
  SuperMatrix B;
  SuperMatrix X;
  SuperMatrix L;
  SuperMatrix U;
  size_t nrhs;
  int *perm_r;
  int *perm_c;
  int *rind;
  int *cptr;
  double *R;
  double *C;
  double *ferr;
  double *berr;
  int nprocs; /* number of processors */

  gsl_vector *rhs_copy;
  gsl_matrix *B_copy; /* copy of right hand side matrix */
  gsl_matrix *X_copy; /* copy of solution matrix */

  double *rnorm; /* residual norm ||A*x - b||, size nrhs */
  double rcond;

  superlumt_options_t options;
} slu_workspace;

/* Prototypes */

slu_workspace *slu_alloc(const size_t size1, const size_t size2, const int nprocs, const size_t nrhs);
void slu_free(slu_workspace *w);
double slu_residual(slu_workspace *w);
int slu_proc(const gsl_spmatrix *A, const gsl_matrix *B, gsl_matrix *X, slu_workspace *w);

#endif /* INCLUDED_superlu_h */
