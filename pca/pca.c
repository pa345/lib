/*
 * pca.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_green.h>

#include "io.h"
#include "pca.h"

static int pca_calc_G(const double b, pca_workspace *w);
static double pca_pc_chi(const size_t pcidx, const double theta, const double phi,
                         pca_workspace *w);
static double pca_pc_calc_chi(const double b, const double theta, const double phi,
                              const gsl_vector *k, pca_workspace *w);

pca_workspace *
pca_alloc(const pca_source_t source_type)
{
  pca_workspace *w;
  char buf[2048];
  size_t i;

  w = calloc(1, sizeof(pca_workspace));
  if (!w)
    return 0;

  w->R = R_EARTH_KM;
  w->ut = 0;
  w->source_type = source_type;

  if (source_type == PCA_SRC_EXTERNAL)
    {
      /* read knm coefficients */

      w->b = w->R + 110.0;

      fprintf(stderr, "pca_alloc: reading %s...", PCA_STAGE1_DATA);
      pca_read_data(PCA_STAGE1_DATA, &(w->nmax), &(w->mmax), NULL, NULL);
      fprintf(stderr, "done (nmax = %zu mmax = %zu)\n", w->nmax, w->mmax);

      fprintf(stderr, "pca_alloc: reading %s...", PCA_STAGE1_KNM);
      w->K = pca_read_matrix(PCA_STAGE1_KNM);
      fprintf(stderr, "done (%zu-by-%zu matrix)\n", w->K->size1, w->K->size2);
    }
  else if (source_type == PCA_SRC_INDUCED)
    {
      /* read gnm coefficients */

      /*XXXw->b = w->R - 10.0;*/
      w->b = w->R - 110.0;

      fprintf(stderr, "pca_alloc: reading %s...", PCA_STAGE1_DATA_IND);
      pca_read_data(PCA_STAGE1_DATA_IND, &(w->nmax), &(w->mmax), NULL, NULL);
      fprintf(stderr, "done (nmax = %zu mmax = %zu)\n", w->nmax, w->mmax);

      fprintf(stderr, "pca_alloc: reading %s...", PCA_STAGE1_GNM);
      w->K = pca_read_matrix(PCA_STAGE1_GNM);
      fprintf(stderr, "done (%zu-by-%zu matrix)\n", w->K->size1, w->K->size2);
    }

  w->green_workspace_p = green_alloc(w->nmax, w->mmax, w->R);
  w->nnm = green_nnm(w->green_workspace_p);
  w->X = malloc(w->nnm * sizeof(double));
  w->Y = malloc(w->nnm * sizeof(double));
  w->Z = malloc(w->nnm * sizeof(double));
  w->work = gsl_vector_alloc(w->nnm);

  /* read SVD components for each UT hour */
  fprintf(stderr, "pca_alloc: reading singular values and left singular vectors...");

  if (source_type == PCA_SRC_EXTERNAL)
    {
      for (i = 0; i < 24; ++i)
        {
          sprintf(buf, "%s_%02zuUT.dat", PCA_STAGE2B_SVAL, i);
          w->S[i] = pca_read_vector(buf);

          sprintf(buf, "%s_%02zuUT.dat", PCA_STAGE2B_U, i);
          w->U[i] = pca_read_matrix(buf);

          w->G[i] = gsl_matrix_alloc(w->nnm, w->nnm);

          sprintf(buf, "%s_%02zuUT.dat", PCA_STAGE2B_MU, i);
          w->mu[i] = pca_read_vector(buf);

          /* convert mu in knm to gnm */
          w->mu_gnm[i] = gsl_vector_alloc(w->nnm);
          green_k2g(w->b, w->mu[i], w->mu_gnm[i], w->green_workspace_p);
        }

      /* calculate gnm from knm coefficients, so we can compute field above current shell b */
      pca_calc_G(w->b, w);
    }
  else if (source_type == PCA_SRC_INDUCED)
    {
      for (i = 0; i < 24; ++i)
        {
          sprintf(buf, "%s_%02zuUT.dat", PCA_STAGE2B_SVAL_IND, i);
          w->S[i] = pca_read_vector(buf);

          sprintf(buf, "%s_%02zuUT.dat", PCA_STAGE2B_U_IND, i);
          w->U[i] = pca_read_matrix(buf);

          sprintf(buf, "%s_%02zuUT.dat", PCA_STAGE2B_MU_IND, i);
          w->mu[i] = pca_read_vector(buf);
        }
    }

  fprintf(stderr, "done\n");

  return w;
}

void
pca_free(pca_workspace *w)
{
  size_t i;

  if (w->K)
    gsl_matrix_free(w->K);

  for (i = 0; i < 24; ++i)
    {
      if (w->S[i])
        gsl_vector_free(w->S[i]);

      if (w->U[i])
        gsl_matrix_free(w->U[i]);

      if (w->G[i])
        gsl_matrix_free(w->G[i]);

      if (w->mu[i])
        gsl_vector_free(w->mu[i]);

      if (w->mu_gnm[i])
        gsl_vector_free(w->mu_gnm[i]);
    }

  if (w->green_workspace_p)
    green_free(w->green_workspace_p);

  if (w->X)
    free(w->X);

  if (w->Y)
    free(w->Y);

  if (w->Z)
    free(w->Z);

  if (w->work)
    gsl_vector_free(w->work);

  free(w);
}

/*
pca_set_UT()
  Set UT hour for PCA modes

Inputs: ut - UT hour in [0,23]
        w  - workspace
*/

int
pca_set_UT(const size_t ut, pca_workspace *w)
{
  w->ut = ut;

  if (w->ut >= 24)
    w->ut -= 24;

  return 0;
}

/*
pca_sval()
  Return some singular values for current UT hour in a vector

Inputs: T - (output) vector of singular values, length of vector
            determines how many singular values are returned,
            sorted largest to smallest
        w - workspace
*/

int
pca_sval(gsl_vector * T, pca_workspace *w)
{
  const gsl_vector *S = w->S[w->ut];
  const size_t nsing = T->size; /* number of requested singular values */

  if (nsing > S->size)
    {
      GSL_ERROR ("T has wrong size", GSL_EBADLEN);
    }
  else
    {
      /* copy largest nsing singular values into T */
      gsl_vector_const_view v = gsl_vector_const_subvector(S, 0, nsing);
      gsl_vector_memcpy(T, &v.vector);
      return GSL_SUCCESS;
    }
}

/*
pca_variance()
  Print cumulative variance based on singular values to a file.
Also determine number of singular vectors required to explain
a desired threshold of variance in the data.

Inputs: filename - output file
        S        - singular values
        thresh   - desired variance threshold in [0,1]
        nsing    - (output) number of singular values required to
                   explain a variance of 'thresh'
*/

int
pca_variance(const char *filename, const double thresh, size_t * nsing,
             pca_workspace *w)
{
  int s = 0;
  const gsl_vector *S = w->S[w->ut];
  const size_t n = S->size;
  FILE *fp = fopen(filename, "w");
  size_t i;
  double sum_all = 0.0;
  double cumsum = 0.0;
  size_t nvec = 0;

  i = 1;
  fprintf(fp, "# %02zu UT\n", w->ut);
  fprintf(fp, "# Field %zu: eigenvalue number i\n", i++);
  fprintf(fp, "# Field %zu: ( sum_{j=1}^i lambda_j ) / ( sum_j lambda_j )\n", i++);

  /* compute sum_i lambda_i, with lambda_i = sigma_i^2 */
  gsl_blas_ddot(S, S, &sum_all);

  for (i = 0; i < n; ++i)
    {
      double sigma = gsl_vector_get(S, i);
      double cumvar;

      cumsum += sigma * sigma;
      cumvar = cumsum / sum_all;

      /* compute number of singular vectors / eigenvectors needed
       * to explain desired variance threshold */
      if (nvec == 0 && cumvar > thresh)
        nvec = i + 1;

      fprintf(fp, "%zu %.12e\n", i + 1, cumvar);
    }

  fclose(fp);

  *nsing = nvec;

  return s;
}

/*
pca_print_pc_map()
  Print lat/lon map of one PC

Inputs: filename - output file
        r        - radius of shell for magnetic field maps
        pcidx    - which principal component, in [0,nnm-1]
        w        - workspace
*/

int
pca_print_pc_map(const char *filename, const double r, const size_t pcidx,
                 pca_workspace *w)
{
  int s = 0;
  FILE *fp;
  double lat, lon;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_print_pc_map: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Principal component %zu/%zu\n", pcidx + 1, w->nnm);
  fprintf(fp, "# Current shell radius: %g km\n", w->b);
  fprintf(fp, "# Magnetic field shell radius: %g km\n", r);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: chi (kA / nT)\n", i++);
  fprintf(fp, "# Field %zu: B_x (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: B_y (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: B_z (dimensionless)\n", i++);

  for (lon = -180.0; lon <= 180.0; lon += 1.0)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -89.9; lat <= 89.9; lat += 1.0)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double B_pc[3];
          double chi;

          /* compute current stream function for this PC */
          chi = pca_pc_chi(pcidx, theta, phi, w);

          /* compute magnetic field vector for this PC */
          pca_pc_B(pcidx, r, theta, phi, B_pc, w);

          fprintf(fp, "%f %f %f %f %f %f\n",
                  lon,
                  lat,
                  chi,
                  B_pc[0],
                  B_pc[1],
                  B_pc[2]);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/*
pca_print_mean_map()
  Print lat/lon map of mean ionospheric signal

Inputs: filename - output file
        r        - radius of shell for magnetic field maps
        w        - workspace
*/

int
pca_print_mean_map(const char *filename, const double r, pca_workspace *w)
{
  int s = 0;
  FILE *fp;
  double lat, lon;
  size_t i;
  gsl_vector *mu = w->mu[w->ut];

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "pca_print_mean_map: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Mean component\n");
  fprintf(fp, "# Current shell radius: %g km\n", w->b);
  fprintf(fp, "# Magnetic field shell radius: %g km\n", r);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: chi (kA)\n", i++);
  fprintf(fp, "# Field %zu: B_x (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z (nT)\n", i++);

  for (lon = -180.0; lon <= 180.0; lon += 1.0)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -89.9; lat <= 89.9; lat += 1.0)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double B[3];
          double chi;

          /* compute current stream function */
          chi = pca_pc_calc_chi(w->b, theta, phi, mu, w);

          /* compute magnetic field vector */
          pca_mean_B(r, theta, phi, B, w);

          fprintf(fp, "%f %f %f %f %f %f\n",
                  lon,
                  lat,
                  chi,
                  B[0],
                  B[1],
                  B[2]);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/*
pca_print_map()
  Print lat/lon map of linear combination of PCs

Inputs: fp    - output file
        r     - radius of shell for magnetic field maps
        alpha - vector specifying linear combination of PCs
        w     - workspace
*/

int
pca_print_map(FILE *fp, const double r, const gsl_vector * alpha,
              pca_workspace *w)
{
  int s = 0;
  double lat, lon;
  size_t i;

  i = 1;
  fprintf(fp, "# Number of principal components: %zu/%zu\n", alpha->size, w->nnm);
  fprintf(fp, "# Current shell radius: %g km\n", w->b);
  fprintf(fp, "# Magnetic field shell radius: %g km\n", r);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: chi (kA / nT)\n", i++);
  fprintf(fp, "# Field %zu: B_x (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z (nT)\n", i++);

  for (lon = -180.0; lon <= 180.0; lon += 1.0)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -89.9; lat <= 89.9; lat += 1.0)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double B[3];
          double chi;

          /* compute current stream function */
          chi = pca_chi(alpha, theta, phi, w);

          /* compute magnetic field vector */
          pca_B(alpha, r, theta, phi, B, w);

          fprintf(fp, "%f %f %f %f %f %f\n",
                  lon,
                  lat,
                  chi,
                  B[0],
                  B[1],
                  B[2]);
        }

      fprintf(fp, "\n");
    }

  fprintf(fp, "\n\n");

  fflush(fp);

  return s;
}

/*
pca_pc_B()
  Compute the (unit) magnetic field vector at a given
point (r,theta,phi) due to a single PC

Inputs: pcidx - PC index in [0,nnm-1]
        r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        B     - (output) unit magnetic field vector B for pcidx
        w     - workspace
*/

int
pca_pc_B(const size_t pcidx, const double r, const double theta, const double phi,
         double B[3], pca_workspace *w)
{
  int s = 0;
  const size_t nnm = w->nnm;
  gsl_vector_view Xv = gsl_vector_view_array(w->X, nnm);
  gsl_vector_view Yv = gsl_vector_view_array(w->Y, nnm);
  gsl_vector_view Zv = gsl_vector_view_array(w->Z, nnm);

  if (w->source_type == PCA_SRC_EXTERNAL)
    {
      /*
       * If r < b, the current shell is an external source so
       * we can directly use the knm coefficients in the U matrix.
       *
       * If r > b, the current shell is an internal source, and
       * we must first compute the gnm coefficients from knm (done
       * in pca_alloc via pca_calc_G), and then use internal Green's functions
       * for the dot product
       */
      if (r <= w->b)
        {
          gsl_vector_const_view pc_knm = gsl_matrix_const_column(w->U[w->ut], pcidx);

          green_calc_ext(r, theta, phi, w->X, w->Y, w->Z, w->green_workspace_p);

          gsl_blas_ddot(&pc_knm.vector, &Xv.vector, &B[0]);
          gsl_blas_ddot(&pc_knm.vector, &Yv.vector, &B[1]);
          gsl_blas_ddot(&pc_knm.vector, &Zv.vector, &B[2]);
        }
      else
        {
          gsl_vector_const_view pc_gnm = gsl_matrix_const_column(w->G[w->ut], pcidx);

          green_calc_int(r, theta, phi, w->X, w->Y, w->Z, w->green_workspace_p);

          gsl_blas_ddot(&pc_gnm.vector, &Xv.vector, &B[0]);
          gsl_blas_ddot(&pc_gnm.vector, &Yv.vector, &B[1]);
          gsl_blas_ddot(&pc_gnm.vector, &Zv.vector, &B[2]);
        }
    }
  else if (w->source_type == PCA_SRC_INDUCED)
    {
      /*
       * For an induced field source, we should always be above the current shell (which
       * is somewhere inside the Earth). Calculate internal Green's functions to compute
       * the field
       */
      if (r >= w->b)
        {
          gsl_vector_const_view pc_gnm = gsl_matrix_const_column(w->U[w->ut], pcidx);

          green_calc_int(r, theta, phi, w->X, w->Y, w->Z, w->green_workspace_p);

          gsl_blas_ddot(&pc_gnm.vector, &Xv.vector, &B[0]);
          gsl_blas_ddot(&pc_gnm.vector, &Yv.vector, &B[1]);
          gsl_blas_ddot(&pc_gnm.vector, &Zv.vector, &B[2]);
        }
      else
        {
          fprintf(stderr, "pca_pc_B: error: r < b for induced field\n");
          return -1;
        }
    }

  return s;
}

/*
pca_B()
  Compute total magnetic field vector at a given point, using
a linear combination of principal components:

B(r) = sum_i alpha_i B_i(r)

where B_i(r) is the unit magnetic field of the ith principal component:

B_i(r) = [ dX^T(r) ; dY^T(r) ; dZ^T(r) ] U_i, 0 <= i < nnm

Inputs: alpha - coefficients of sum
        r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        B     - (output) magnetic field vector
        w     - workspace
*/

int
pca_B(const gsl_vector *alpha, const double r, const double theta, const double phi,
      double B[3], pca_workspace *w)
{
  const size_t p = alpha->size; /* number of principal components to use to compute B */
  const size_t nnm = w->nnm;
  gsl_vector_view Xv = gsl_vector_view_array(w->X, nnm);
  gsl_vector_view Yv = gsl_vector_view_array(w->Y, nnm);
  gsl_vector_view Zv = gsl_vector_view_array(w->Z, nnm);
  gsl_matrix *U;
  size_t i;

  gsl_vector_set_zero(w->work);

  if (w->source_type == PCA_SRC_EXTERNAL)
    {
      if (r <= w->b)
        {
          /* external current source, use U matrix and external Green's functions */
          U = w->U[w->ut];
          green_calc_ext(r, theta, phi, w->X, w->Y, w->Z, w->green_workspace_p);
        }
      else
        {
          /* internal current source, use G matrix and internal Green's functions */
          U = w->G[w->ut];
          green_calc_int(r, theta, phi, w->X, w->Y, w->Z, w->green_workspace_p);
        }
    }
  else if (w->source_type == PCA_SRC_INDUCED)
    {
      /* induced source, use U matrix and internal Green's functions */
      U = w->U[w->ut];
      green_calc_int(r, theta, phi, w->X, w->Y, w->Z, w->green_workspace_p);

      /* sanity check */
      if (r < w->b)
        {
          fprintf(stderr, "pca_B: error: r < b for induced source\n");
          return -1;
        }
    }

  /* compute: work = sum_i alpha_i U_i over largest p principal components */
  for (i = 0; i < p; ++i)
    {
      double ai = gsl_vector_get(alpha, i);
      gsl_vector_view Ui = gsl_matrix_column(U, i);

      gsl_blas_daxpy(ai, &Ui.vector, w->work);
    }

  /* compute B = [ dX^T ; dY^T ; dZ^T ] * sum_i alpha_i U_i */
  gsl_blas_ddot(w->work, &Xv.vector, &B[0]);
  gsl_blas_ddot(w->work, &Yv.vector, &B[1]);
  gsl_blas_ddot(w->work, &Zv.vector, &B[2]);

  return 0;
}

/*
pca_mean_B()
  Compute the magnetic field vector at a given
point (r,theta,phi) due to the mean

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        B     - (output) mean magnetic field vector B in NEC (nT)
        w     - workspace
*/

int
pca_mean_B(const double r, const double theta, const double phi,
           double B[3], pca_workspace *w)
{
  int s = 0;
  const size_t nnm = w->nnm;
  gsl_vector_view Xv = gsl_vector_view_array(w->X, nnm);
  gsl_vector_view Yv = gsl_vector_view_array(w->Y, nnm);
  gsl_vector_view Zv = gsl_vector_view_array(w->Z, nnm);
  gsl_vector *mu;

  if (w->source_type == PCA_SRC_EXTERNAL)
    {
      /*
       * If r < b, the current shell is an external source so
       * we can directly use the knm coefficients in the U matrix.
       *
       * If r > b, the current shell is an internal source, and
       * we must first compute the gnm coefficients from knm (done
       * in pca_alloc via pca_calc_G), and then use internal Green's functions
       * for the dot product
       */
      if (r <= w->b)
        {
          mu = w->mu[w->ut];
          green_calc_ext(r, theta, phi, w->X, w->Y, w->Z, w->green_workspace_p);
        }
      else
        {
          mu = w->mu_gnm[w->ut];
          green_calc_int(r, theta, phi, w->X, w->Y, w->Z, w->green_workspace_p);
        }
    }
  else if (w->source_type == PCA_SRC_INDUCED)
    {
      mu = w->mu[w->ut];
      green_calc_int(r, theta, phi, w->X, w->Y, w->Z, w->green_workspace_p);

      /* sanity check */
      if (r < w->b)
        {
          fprintf(stderr, "pca_mean_B: error: r < b for induced source\n");
          return -1;
        }
    }

  gsl_blas_ddot(mu, &Xv.vector, &B[0]);
  gsl_blas_ddot(mu, &Yv.vector, &B[1]);
  gsl_blas_ddot(mu, &Zv.vector, &B[2]);

  return s;
}

/*
pca_K()
  Compute sheet current density at a given point, using
a linear combination of principal components

Inputs: alpha - coefficients of sum over PCs
        theta - colatitude (radians)
        phi   - longitude (radians)
        K     - (output) sheet current density in A/km
        w     - workspace

Return: success/error
*/

int
pca_K(const gsl_vector *alpha, const double theta, const double phi, double K[3], pca_workspace *w)
{
  int status;
  const size_t p = alpha->size; /* number of principal components to use to compute K */
  size_t i;
  gsl_matrix *U;

  gsl_vector_set_zero(w->work);

  if (w->source_type == PCA_SRC_EXTERNAL)
    {
      U = w->G[w->ut];
    }
  else
    {
      U = w->U[w->ut];
    }

  /* compute: work = sum_i alpha_i G_i over largest p principal components in gnm basis */
  for (i = 0; i < p; ++i)
    {
      double ai = gsl_vector_get(alpha, i);
      gsl_vector_view Gi = gsl_matrix_column(U, i);

      gsl_blas_daxpy(ai, &Gi.vector, w->work);
    }

  status = green_eval_sheet_int(w->b, theta, phi, w->work, K, w->green_workspace_p);

  return status;
}

/*
pca_mean_K()
  Compute sheet current density at a given point, using mean current

Inputs: theta - colatitude (radians)
        phi   - longitude (radians)
        K     - (output) sheet current density in A/km
        w     - workspace

Return: success/error
*/

int
pca_mean_K(const double theta, const double phi, double K[3], pca_workspace *w)
{
  int status;
  gsl_vector *mu;

  if (w->source_type == PCA_SRC_EXTERNAL)
    {
      mu = w->mu_gnm[w->ut];
    }
  else if (w->source_type == PCA_SRC_INDUCED)
    {
      mu = w->mu[w->ut];
    }

  status = green_eval_sheet_int(w->b, theta, phi, mu, K, w->green_workspace_p);

  return status;
}

/*
pca_chi()
  Compute total current stream function at a given point, using
a linear combination of principal components:

chi(r) = sum_i alpha_i chi_i(r)

where chi_i(r) is the unit stream function of the ith principal component:

Inputs: alpha - coefficients of sum
        theta - colatitude (radians)
        phi   - longitude (radians)
        w     - workspace

Return: current stream function
*/

double
pca_chi(const gsl_vector *alpha, const double theta, const double phi, pca_workspace *w)
{
  const size_t p = alpha->size; /* number of principal components to use to compute B */
  size_t i;
  double chi;

  gsl_vector_set_zero(w->work);

  /* compute: work = sum_i alpha_i U_i over largest p principal components */
  for (i = 0; i < p; ++i)
    {
      double ai = gsl_vector_get(alpha, i);
      gsl_vector_view Ui = gsl_matrix_column(w->U[w->ut], i);

      gsl_blas_daxpy(ai, &Ui.vector, w->work);
    }

  chi = pca_pc_calc_chi(w->b, theta, phi, w->work, w);

  return chi;
}

/*
pca_mean_chi()
  Compute current stream function for mean

Inputs: theta - colatitude (radians)
        phi   - longitude (radians)
        w     - workspace
*/

double
pca_mean_chi(const double theta, const double phi, pca_workspace *w)
{
  gsl_vector *mu = w->mu[w->ut];
  return pca_pc_calc_chi(w->b, theta, phi, mu, w);
}

/*
pca_calc_G()
  Convert PCs in a basis of external Gauss coefficients representing
a current shell at radius b into internal Gauss coefficients, suitable
for computing the magnetic field above the current shell.

gnm = - (n / (n+1)) (b/R)^{2n + 1} knm

Inputs: b - radius of current shell (km)
        w - workspace
*/

static int
pca_calc_G(const double b, pca_workspace *w)
{
  size_t n, ut;
  green_workspace *green_p = w->green_workspace_p;

  for (ut = 0; ut < 24; ++ut)
    {
      for (n = 1; n <= w->nmax; ++n)
        {
          int M = (int) GSL_MIN(n, w->mmax);
          int m;
          double nfac = -(double)n / (n + 1.0);
          double rfac = pow(b / w->R, 2.0*n + 1.0);

          for (m = -M; m <= M; ++m)
            {
              size_t cidx = green_nmidx(n, m, green_p);
              gsl_vector_view k = gsl_matrix_row(w->U[ut], cidx);
              gsl_vector_view g = gsl_matrix_row(w->G[ut], cidx);

              gsl_vector_memcpy(&g.vector, &k.vector);
              gsl_vector_scale(&g.vector, nfac * rfac);
            }
        }
    }

  return 0;
}

/*
pca_pc_chi()
  Compute current stream function for a single PC

Inputs: pcidx - index of PC
        theta - colatitude (radians)
        phi   - longitude (radians)
        w     - workspace
*/

static double
pca_pc_chi(const size_t pcidx, const double theta, const double phi, pca_workspace *w)
{
  gsl_vector_const_view pc = gsl_matrix_const_column(w->U[w->ut], pcidx);
  return pca_pc_calc_chi(w->b, theta, phi, &pc.vector, w);
}

/*
pca_pc_calc_chi()
  Compute current stream function for a single PC

chi = -(b/mu0) (b/R) sum_{nm} qnm Ynm

where:

qnm = (2n + 1)/(n + 1) (b/R)^{n-2} knm

Inputs: b     - radius of current shell (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        k     - coefficient vector knm for this PC (or gnm for induced source)
        w     - workspace
*/

static double
pca_pc_calc_chi(const double b, const double theta, const double phi,
                const gsl_vector *k, pca_workspace *w)
{
  const double mu0 = 400.0 * M_PI; /* units of nT / (kA km^{-1}) */
  const size_t nmax = w->nmax;
  green_workspace *green_p = w->green_workspace_p;
  double *Plm = green_p->Plm;
  double *dPlm = green_p->dPlm;
  size_t n;
  double chi = 0.0;

  /* compute associated legendres */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax, cos(theta), Plm, dPlm);

  if (w->source_type == PCA_SRC_EXTERNAL)
    {
      /* external source, use:
       *
       * qnm = (2n+1)/(n+1) (b/R)^{n-2} knm
       */

      const double ratio = b / w->R;
      double rfac = 1.0;

      for (n = 1; n <= nmax; ++n)
        {
          int M = (int) GSL_MIN(n, w->mmax);
          int m;
          double nfac = (2.0 * n + 1.0) / (n + 1.0);

          for (m = -M; m <= M; ++m)
            {
              int mabs = abs(m);
              size_t cidx = green_nmidx(n, m, green_p);
              double knm = gsl_vector_get(k, cidx);
              double qnm = nfac * rfac * knm;
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              double Snm = Plm[pidx];

              if (m < 0)
                Snm *= sin(mabs * phi);
              else
                Snm *= cos(mabs * phi);

              chi += qnm * Snm;
            }

          /* (b/R)^{n-2} */
          rfac *= ratio;
        }
    }
  else if (w->source_type == PCA_SRC_INDUCED)
    {
      /* induced (internal) source, use:
       *
       * qnm = -(2n+1)/n (R/b)^{n+3} gnm
       */
      const double ratio = w->R / b;
      double rfac = ratio * ratio * ratio;

      for (n = 1; n <= nmax; ++n)
        {
          int M = (int) GSL_MIN(n, w->mmax);
          int m;
          double nfac = -(2.0 * n + 1.0) / (double)n;

          /* (R/b)^{n+3} */
          rfac *= ratio;

          for (m = -M; m <= M; ++m)
            {
              int mabs = abs(m);
              size_t cidx = green_nmidx(n, m, green_p);
              double gnm = gsl_vector_get(k, cidx);
              double qnm = nfac * rfac * gnm;
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              double Snm = Plm[pidx];

              if (m < 0)
                Snm *= sin(mabs * phi);
              else
                Snm *= cos(mabs * phi);

              chi += qnm * Snm;
            }
        }
    }

  chi *= (b / w->R);
  chi *= -(b / mu0); /* units of kA */

  return chi;
}
