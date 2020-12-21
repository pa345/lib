/*
 * stage3b.c
 *
 * Similar to stage3 but for time domain PCA
 *
 * 1. Print cumulative variance curve to determine how
 *    many eigenvectors are needed
 * 2. Solve U*alpha = knm where U is a matrix of principal
 *    eigenvectors
 * 3. Print lat/lon maps of magnetic field and current flow
 *    due to first 3 PCs
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_legendre.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_green.h>
#include <mainlib/ml_lapack_wrapper.h>

#include "io.h"
#include "pca.h"

static double
norm_fro(const gsl_matrix * m)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;
  double sum = 0.0;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double mij = gsl_matrix_get(m, i, j);
          sum += fabs(mij);
        }
    }

  sum = sqrt(sum);

  return sum;
}

/*
print_variance()
  Print cumulative variance based on singular values to a file.
Also determine number of singular vectors required to explain
a desired threshold of variance in the data.
*/

int
print_variance(const char *filename, const gsl_vector *S, const double thresh,
               size_t * P)
{
  int s = 0;
  const size_t n = S->size;
  FILE *fp = fopen(filename, "w");
  size_t i;
  double sum_all = 0.0;
  double cumsum = 0.0;
  size_t nvec = 0;

  i = 1;
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

  *P = nvec;

  return s;
}

/*
solve_PCA()
  Solve PCA problem by computing alpha, such that

U alpha = Q

in a least squares sense, where U are the principal eigenvectors
of the SDM

Inputs: P        - number of eigenvectors of SDM to use
        knm      - SH time series matrix, nnm-by-nt
        U        - left singular vector (eigenvector) matrix, nnm-by-nnm
        alpha    - (output) alpha matrix, P-by-nt
        knmt     - (output) k~ = U*alpha, nnm-by-nt
*/

int
solve_PCA(const size_t P, const gsl_matrix * knm,
          const gsl_matrix * U, gsl_matrix * alpha,
          gsl_matrix * knmt)
{
  int status = 0;
  const size_t nt = knm->size2; /* number of time stamps */
  const size_t nnm = U->size1;
  gsl_matrix *R;                /* R = knm - U*alpha */
  struct timeval tv0, tv1;
  double residual;           /* || knm - U*alpha || */
  int rank;

  /* select largest P eigenvectors of SDM */
  gsl_matrix_const_view Uv = gsl_matrix_const_submatrix(U, 0, 0, nnm, P);

  /* solve: U*alpha = Q */
  fprintf(stderr, "solve_PCA: solving PCA problem for alpha...");
  gettimeofday(&tv0, NULL);
  status = lapack_lls(&Uv.matrix, knm, alpha, &rank);
  gettimeofday(&tv1, NULL);

  /* compute: knm~ = U*alpha */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Uv.matrix, alpha, 0.0, knmt);

  /* compute: R = knm - knm~ */
  R = gsl_matrix_alloc(nnm, nt);
  gsl_matrix_memcpy(R, knm);
  gsl_matrix_sub(R, knmt);

  residual = norm_fro(R);

  fprintf(stderr, "done (%g seconds, status = %d, rank = %d, residual = %.12e)\n",
          time_diff(tv0, tv1), status, rank, residual);

  gsl_matrix_free(R);

  return status;
}

int
print_pc_maps(const char *filename, const gsl_matrix * U,
              green_workspace *green_p)
{
  int s = 0;
  FILE *fp;
  const double r = R_EARTH_KM + 0.0;   /* radius of magnetic field maps */
  const double b = R_EARTH_KM + 110.0; /* radius of current shell */
  double lat, lon;
  size_t nnm = green_nnm(green_p);
  double *X = malloc(nnm * sizeof(double));
  double *Y = malloc(nnm * sizeof(double));
  double *Z = malloc(nnm * sizeof(double));
  gsl_vector_view Xv = gsl_vector_view_array(X, nnm);
  gsl_vector_view Yv = gsl_vector_view_array(Y, nnm);
  gsl_vector_view Zv = gsl_vector_view_array(Z, nnm);
  gsl_vector_const_view pc1 = gsl_matrix_const_column(U, 0);
  gsl_vector_const_view pc2 = gsl_matrix_const_column(U, 1);
  gsl_vector_const_view pc3 = gsl_matrix_const_column(U, 2);
  gsl_vector *gnm1 = gsl_vector_alloc(nnm);
  gsl_vector *gnm2 = gsl_vector_alloc(nnm);
  gsl_vector *gnm3 = gsl_vector_alloc(nnm);
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_pc_maps: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  /* compute gnm coefficients in case r > b */
  green_k2g(b, &pc1.vector, gnm1, green_p);
  green_k2g(b, &pc2.vector, gnm2, green_p);
  green_k2g(b, &pc3.vector, gnm3, green_p);

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: PC1 chi (kA / nT)\n", i++);
  fprintf(fp, "# Field %zu: PC1 B_x (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: PC1 B_y (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: PC1 B_z (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: PC2 chi (kA / nT)\n", i++);
  fprintf(fp, "# Field %zu: PC2 B_x (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: PC2 B_y (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: PC2 B_z (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: PC3 chi (kA / nT)\n", i++);
  fprintf(fp, "# Field %zu: PC3 B_x (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: PC3 B_y (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: PC3 B_z (dimensionless)\n", i++);

  for (lon = -180.0; lon <= 180.0; lon += 1.0)
    {
      double phi = lon * M_PI / 180.0;

      for (lat = -89.9; lat <= 89.9; lat += 1.0)
        {
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double B_pc1[3], B_pc2[3], B_pc3[3];
          double chi1, chi2, chi3;

          green_ext_chi(b, theta, phi, &pc1.vector, &chi1, green_p);
          green_ext_chi(b, theta, phi, &pc2.vector, &chi2, green_p);
          green_ext_chi(b, theta, phi, &pc3.vector, &chi3, green_p);

          /*
           * If r < b, the current shell is an external source so
           * we can directly use the knm coefficients in the U matrix.
           *
           * If r > b, the current shell is an internal source, and
           * we must first compute the gnm coefficients from knm (done
           * above via green_k2g, and then use internal Green's functions
           * for the dot product
           */
          if (r <= b)
            {
              green_ext(r, theta, phi, &Xv.vector, &Yv.vector, &Zv.vector, green_p);

              gsl_blas_ddot(&pc1.vector, &Xv.vector, &B_pc1[0]);
              gsl_blas_ddot(&pc1.vector, &Yv.vector, &B_pc1[1]);
              gsl_blas_ddot(&pc1.vector, &Zv.vector, &B_pc1[2]);

              gsl_blas_ddot(&pc2.vector, &Xv.vector, &B_pc2[0]);
              gsl_blas_ddot(&pc2.vector, &Yv.vector, &B_pc2[1]);
              gsl_blas_ddot(&pc2.vector, &Zv.vector, &B_pc2[2]);

              gsl_blas_ddot(&pc3.vector, &Xv.vector, &B_pc3[0]);
              gsl_blas_ddot(&pc3.vector, &Yv.vector, &B_pc3[1]);
              gsl_blas_ddot(&pc3.vector, &Zv.vector, &B_pc3[2]);
            }
          else
            {
              green_int(r, theta, phi, &Xv.vector, &Yv.vector, &Zv.vector, green_p);

              gsl_blas_ddot(gnm1, &Xv.vector, &B_pc1[0]);
              gsl_blas_ddot(gnm1, &Yv.vector, &B_pc1[1]);
              gsl_blas_ddot(gnm1, &Zv.vector, &B_pc1[2]);

              gsl_blas_ddot(gnm2, &Xv.vector, &B_pc2[0]);
              gsl_blas_ddot(gnm2, &Yv.vector, &B_pc2[1]);
              gsl_blas_ddot(gnm2, &Zv.vector, &B_pc2[2]);

              gsl_blas_ddot(gnm3, &Xv.vector, &B_pc3[0]);
              gsl_blas_ddot(gnm3, &Yv.vector, &B_pc3[1]);
              gsl_blas_ddot(gnm3, &Zv.vector, &B_pc3[2]);
            }

          fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                  lon,
                  lat,
                  chi1,
                  B_pc1[0],
                  B_pc1[1],
                  B_pc1[2],
                  chi2,
                  B_pc2[0],
                  B_pc2[1],
                  B_pc2[2],
                  chi3,
                  B_pc3[0],
                  B_pc3[1],
                  B_pc3[2]);
        }

      fprintf(fp, "\n");
    }

  free(X);
  free(Y);
  free(Z);
  fclose(fp);
  gsl_vector_free(gnm1);
  gsl_vector_free(gnm2);
  gsl_vector_free(gnm3);

  return s;
}

int
main(int argc, char *argv[])
{
  const size_t nmax = 60;
  const size_t mmax = GSL_MIN(nmax, 30);
  const double R = R_EARTH_KM;
  green_workspace *green_p = green_alloc(nmax, mmax, R);
  pca_workspace *pca_p = pca_alloc(PCA_SRC_EXTERNAL);
  char *knm_file = "data/stage1_knm.dat";

  char *pc_file = "pc_time.txt";
  char *recon_file = "recon_time.txt";
  char buf[2048];

  const double var_thresh = 0.99;

  gsl_vector *S;             /* singular values of SDM */
  gsl_matrix *U;             /* left singular vectors of SDM */
  gsl_matrix *alpha;         /* alpha matrix, P-by-nt */
  gsl_matrix *knmt;          /* knm~ = U*alpha, nnm-by-nt */

  gsl_matrix *knm;           /* knm(t) matrix */

  size_t nnm;
  size_t nt;                 /* number of time stamps */
  size_t P;                  /* number of principal eigenvectors to use (<= T) */

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          default:
            fprintf(stderr, "Usage: %s <-i stage1_matrix_file>\n", argv[0]);
            break;
        }
    }

  fprintf(stderr, "main: reading knm matrix from %s...", knm_file);
  knm = pca_read_matrix(knm_file);
  fprintf(stderr, "done (%zu-by-%zu matrix read)\n", knm->size1, knm->size2);

  {
    size_t i;

    /* plot a variance curve for each UT to help decide how many eigenvectors to keep */
    for (i = 0; i < 24; ++i)
      {
        pca_set_UT(i, pca_p);

        sprintf(buf, "%s_time_%02zuUT.txt", PCA_STAGE3B_VAR, i);

        fprintf(stderr, "main: writing variance curve for %02zu UT to %s...", i, buf);
        pca_variance(buf, var_thresh, &P, pca_p);
        fprintf(stderr, "done (%zu singular vectors needed to explain %.1f%% of variance)\n",
                P, var_thresh * 100.0);
      }
  }

  pca_set_UT(12, pca_p);
  U = pca_p->U[pca_p->ut];
  S = pca_p->S[pca_p->ut];

  nnm = U->size1;
  nt = knm->size2;

  fprintf(stderr, "main: using %zu largest eigenvectors\n", P);

  alpha = gsl_matrix_alloc(P, nt);
  knmt = gsl_matrix_alloc(nnm, nt);

  /* find alpha such that || knm - U*alpha || is
   * minimized in a least squares sense */
  solve_PCA(P, knm, U, alpha, knmt);

  /* plot reconstructed time series using dominant PCs */
  {
    const size_t n = 3;
    const int m = 1;
    const size_t cidx = green_nmidx(n, m, green_p);
    FILE *fp = fopen(recon_file, "w");
    size_t i;

    fprintf(stderr, "main: writing reconstructed (%zu,%d) time series to %s...",
            n, m, recon_file);

    for (i = 0; i < nt; ++i)
      {
        double t = (double) i;

        fprintf(fp, "%f %f %f\n",
                t / 24.0,
                gsl_matrix_get(knm, cidx, i),
                gsl_matrix_get(knmt, cidx, i));
      }

    fprintf(stderr, "done\n");

    fclose(fp);
  }

  fprintf(stderr, "main: printing principle component maps to %s...",
          pc_file);
  print_pc_maps(pc_file, U, green_p);
  fprintf(stderr, "done\n");

  gsl_matrix_free(alpha);
  gsl_matrix_free(knmt);
  gsl_matrix_free(knm);
  green_free(green_p);
  pca_free(pca_p);

  return 0;
}
