/*
 * superlu_mt.c
 * Patrick Alken
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <superlu_mt/slu_mt_ddefs.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

#include "superlu.h"

/*
slu_alloc()
  Allocate a SuperLU workspace

Inputs: size1  - number of rows in sparse matrix
        size2  - number of columns in sparse matrix
        nprocs - number of processors to use
        nrhs   - number of right hand side vectors

Return: pointer to workspace
*/

slu_workspace *
slu_alloc(const size_t size1, const size_t size2, const int nprocs, const size_t nrhs)
{
  slu_workspace *w;

  w = calloc(1, sizeof(slu_workspace));
  if (!w)
    {
      fprintf(stderr, "slu_alloc failed\n");
      return 0;
    }

  w->perm_r = malloc(size1 * sizeof(int));
  w->perm_c = malloc(size2 * sizeof(int));

  if (w->perm_r == 0 || w->perm_c == 0)
    {
      slu_free(w);
      return 0;
    }

  w->rhs_copy = gsl_vector_alloc(size1);
  w->B_copy = gsl_matrix_alloc(nrhs, size1);
  w->X_copy = gsl_matrix_alloc(nrhs, size2);
  w->cptr = malloc((size2 + 1) * sizeof(int));

  w->R = malloc(size1 * sizeof(double));
  w->C = malloc(size2 * sizeof(double));

  w->ferr = malloc(nrhs * sizeof(double));
  w->berr = malloc(nrhs * sizeof(double));
  w->rnorm = malloc(nrhs * sizeof(double));

  w->size1 = size1;
  w->size2 = size2;

  w->nprocs = nprocs;
  w->nrhs = nrhs;
  w->rcond = 0.0;

  w->options.nprocs = nprocs;
  w->options.fact = EQUILIBRATE;
  w->options.trans = NOTRANS;
  w->options.refact = NO;
  w->options.panel_size = sp_ienv(1);
  w->options.relax = sp_ienv(2);
  w->options.usepr = NO;
  w->options.drop_tol = 0.0;
  w->options.diag_pivot_thresh = 1.0;
  w->options.SymmetricMode = NO;
  w->options.PrintStat = NO;
  w->options.perm_c = w->perm_c;
  w->options.perm_r = w->perm_r;
  w->options.work = NULL;
  w->options.lwork = 0;

  w->options.etree = malloc(size2 * sizeof(int));
  w->options.colcnt_h = malloc(size2 * sizeof(int));
  w->options.part_super_h = malloc(size2 * sizeof(int));

  return (w);
} /* slu_alloc() */

void
slu_free(slu_workspace *w)
{
  if (w->perm_r)
    free(w->perm_r);

  if (w->perm_c)
    free(w->perm_c);

  if (w->rhs_copy)
    gsl_vector_free(w->rhs_copy);

  if (w->B_copy)
    gsl_matrix_free(w->B_copy);

  if (w->X_copy)
    gsl_matrix_free(w->X_copy);

  if (w->cptr)
    free(w->cptr);

  if (w->R)
    free(w->R);

  if (w->C)
    free(w->C);

  if (w->ferr)
    free(w->ferr);

  if (w->berr)
    free(w->berr);

  if (w->rnorm)
    free(w->rnorm);

  if (w->options.etree)
    free(w->options.etree);

  if (w->options.colcnt_h)
    free(w->options.colcnt_h);

  if (w->options.part_super_h)
    free(w->options.part_super_h);

  free(w);
} /* slu_free() */

double
slu_residual(slu_workspace *w)
{
  return w->rnorm[0];
}

/*
slu_proc()
  Solve a sparse linear system A X = B

Inputs: A  - sparse matrix in CCS format
        B  - rhs matrix (length w->size1-by-nrhs)
        X  - (output) solution matrix (length w->size2-by-nrhs)
        w  - workspace

Return: 0 on success, non-zero on error
*/

int
slu_proc(const gsl_spmatrix *A, const gsl_matrix *B, gsl_matrix *X, slu_workspace *w)
{
  if (A->size1 != w->size1 || A->size2 != w->size2)
    {
      GSL_ERROR("sparse matrix does not match workspace", GSL_EBADLEN);
    }
  else if (B->size1 != w->size1 || B->size2 != w->nrhs)
    {
      GSL_ERROR("rhs matrix does not match workspace", GSL_EBADLEN);
    }
  else if (X->size1 != w->size2 || X->size2 != w->nrhs)
    {
      GSL_ERROR("solution matrix does not match workspace", GSL_EBADLEN);
    }
  else if (!GSL_SPMATRIX_ISCCS(A))
    {
      GSL_ERROR("sparse matrix must be in CCS format", GSL_EINVAL);
    }
  else
    {
      const size_t nnz = gsl_spmatrix_nnz(A);
      int info = 0;
      int *rind = malloc(nnz * sizeof(int));
      double *data = malloc(nnz * sizeof(double));
      size_t i;

      /* make copy of rhs matrix */
      gsl_matrix_transpose_memcpy(w->B_copy, B);

      /* have to copy arrays since sizeof(int) != sizeof(size_t) */
      for (i = 0; i < nnz; ++i)
        {
          rind[i] = (int) A->i[i];
          data[i] = A->data[i];
        }

      for (i = 0; i < w->size2 + 1; ++i)
        w->cptr[i] = (int) A->p[i];

      dCreate_CompCol_Matrix(&(w->A),
                             (int) A->size1,
                             (int) A->size2,
                             nnz,
                             data,
                             rind,
                             w->cptr,
                             SLU_NC,
                             SLU_D,
                             SLU_GE);

      /* rhs matrix */
      dCreate_Dense_Matrix(&(w->B),
                           (int) w->size1,
                           (int) w->nrhs,
                           w->B_copy->data,
                           (int) w->size1,
                           SLU_DN,
                           SLU_D,
                           SLU_GE);

      /* solution matrix */
      dCreate_Dense_Matrix(&(w->X),
                           (int) w->size2,
                           (int) w->nrhs,
                           w->X_copy->data,
                           (int) w->size2,
                           SLU_DN,
                           SLU_D,
                           SLU_GE);

      get_perm_c(1, &(w->A), w->perm_c);

      {
        equed_t equed = NOEQUIL;
        double rpg;
        superlu_memusage_t memusage;

        pdgssvx(w->nprocs, &(w->options), &(w->A), w->perm_c, w->perm_r,
                &equed, w->R, w->C, &(w->L), &(w->U), &(w->B), &(w->X), &rpg, &(w->rcond),
	              w->ferr, w->berr, &memusage, &info);
      }

      if (info != 0)
        {
          fprintf(stderr, "slu_proc: error in pdgssv: info = %d\n", info);
          return info; /* error */
        }

      /* store solution matrix */
      gsl_matrix_transpose_memcpy(X, w->X_copy);

      /* compute residual */
      gsl_matrix_transpose_memcpy(w->B_copy, B);

      for (i = 0; i < w->nrhs; ++i)
        {
          gsl_vector_view xi = gsl_matrix_column(X, i);
          gsl_vector_view bi = gsl_matrix_row(w->B_copy, i);

          gsl_spblas_dgemv(CblasNoTrans, -1.0, A, &xi.vector, 1.0, &bi.vector);
          w->rnorm[i] = gsl_blas_dnrm2(&bi.vector);
        }

      Destroy_SuperMatrix_Store(&(w->A));
      Destroy_SuperMatrix_Store(&(w->B));
      Destroy_SuperMatrix_Store(&(w->X));
      Destroy_SuperNode_SCP(&(w->L));
      Destroy_CompCol_NCP(&(w->U));

      free(rind);
      free(data);

      return GSL_SUCCESS;
    }
}
