/* lse.c
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * This module is a wrapper for the LAPACK linear equality-constrained
 * least squares problem (LSE).
 */

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void
dgglse_(int *m, int *n, int *p, double *a, int *lda, double *b,
        int *ldb, double *c, double *d, double *x, double *work,
        int *lwork, int *info);

/*
lse()
  This routine solves

min_x || A x - c ||_2 subject to: B x = d
*/

int
lse(const gsl_matrix *A, const gsl_matrix *B, const gsl_vector *c,
    const gsl_vector *d, gsl_vector *x)
{
  int m = (int) A->size1;
  int n = (int) A->size2;
  int p = (int) B->size1;
  int info;
  int lwork;
  gsl_matrix *aa, *bb;
  gsl_vector *cc, *dd;
  gsl_vector *work;
  double q[1];

  aa = gsl_matrix_alloc(A->size2, A->size1);
  bb = gsl_matrix_alloc(B->size2, B->size1);
  cc = gsl_vector_alloc(c->size);
  dd = gsl_vector_alloc(d->size);

  gsl_matrix_transpose_memcpy(aa, A);
  gsl_matrix_transpose_memcpy(bb, B);
  gsl_vector_memcpy(cc, c);
  gsl_vector_memcpy(dd, d);

  lwork = -1;
  dgglse_(&m,
          &n,
          &p,
          aa->data,
          (int *) &(aa->tda),
          bb->data,
          (int *) &(bb->tda),
          cc->data,
          dd->data,
          x->data,
          q,
          &lwork,
          &info);

  lwork = (int) q[0];
  work = gsl_vector_alloc((size_t) lwork);

  dgglse_(&m,
          &n,
          &p,
          aa->data,
          (int *) &(aa->tda),
          bb->data,
          (int *) &(bb->tda),
          cc->data,
          dd->data,
          x->data,
          work->data,
          &lwork,
          &info);

  gsl_matrix_free(aa);
  gsl_matrix_free(bb);
  gsl_vector_free(cc);
  gsl_vector_free(dd);
  gsl_vector_free(work);

  if (info)
    fprintf(stderr, "ERROR: lse: info = %d\n", info);

  return (info);
} /* lse() */

/*
wlse()
  This routine solves

min_x || A x - c ||_2 subject to: B x = d

where we are performing a weighted least squares fit with weights w
*/

int
wlse(const gsl_matrix *A, const gsl_matrix *B, const gsl_vector *c,
     const gsl_vector *d, const gsl_vector *w, gsl_vector *x)
{
  size_t i;
  int m = (int) A->size1;
  int n = (int) A->size2;
  int p = (int) B->size1;
  int info;
  int lwork;
  gsl_matrix *aa, *bb;
  gsl_vector *cc, *dd;
  gsl_vector *work;
  double q[1];

  aa = gsl_matrix_alloc(A->size2, A->size1);
  bb = gsl_matrix_alloc(B->size2, B->size1);
  cc = gsl_vector_alloc(c->size);
  dd = gsl_vector_alloc(d->size);

  gsl_matrix_transpose_memcpy(aa, A);
  gsl_matrix_transpose_memcpy(bb, B);
  gsl_vector_memcpy(dd, d);

  for (i = 0; i < A->size1; ++i)
    {
      double wi = gsl_vector_get(w, i);
      gsl_vector_view v = gsl_matrix_column(aa, i);

      gsl_vector_scale(&v.vector, sqrt(wi));
      gsl_vector_set(cc, i, gsl_vector_get(c, i) * sqrt(wi));
    }

  lwork = -1;
  dgglse_(&m,
          &n,
          &p,
          aa->data,
          (int *) &(aa->tda),
          bb->data,
          (int *) &(bb->tda),
          cc->data,
          dd->data,
          x->data,
          q,
          &lwork,
          &info);

  lwork = (int) q[0];
  work = gsl_vector_alloc((size_t) lwork);

  dgglse_(&m,
          &n,
          &p,
          aa->data,
          (int *) &(aa->tda),
          bb->data,
          (int *) &(bb->tda),
          cc->data,
          dd->data,
          x->data,
          work->data,
          &lwork,
          &info);

  gsl_matrix_free(aa);
  gsl_matrix_free(bb);
  gsl_vector_free(cc);
  gsl_vector_free(dd);
  gsl_vector_free(work);

  if (info)
    fprintf(stderr, "ERROR: wlse: info = %d\n", info);

  return (info);
} /* wlse() */
