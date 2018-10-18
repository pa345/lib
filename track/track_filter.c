/*
 * track_filter.c
 *
 * Filter/smooth tracks
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include <common/common.h>
#include <common/interp.h>

#include "track.h"

#define TRACK_FILTER_EXCLUDE_LAT        (20.0)

typedef struct
{
  gsl_vector *c;   /* filter coefficients */
  gsl_vector *rhs;
  gsl_matrix *cov;
  gsl_matrix *X;
} track_filter_params;

static int track_filter_proc(const size_t track_idx, track_filter_params *params,
                             satdata_mag *data, track_workspace *w);
static size_t track_filter_matrix_row_int(const double t, const double r, const double theta,
                                          const double phi, gsl_vector *v, track_workspace *w);
static size_t track_filter_matrix_row_ext(const double t, const double r, const double theta,
                                          const double phi, gsl_vector *v, track_workspace *w);
static size_t track_filter_matrix_row(const double t, const double r, const double theta,
                                      const double phi, gsl_vector *v, track_workspace *w);

int
track_filter(const char *filename, track_workspace *w)
{
  int s = 0;
  const size_t ntot = 10000; /* max data in a track */
  const size_t p = 20;       /* model coefficients */
  size_t i, j;
  satdata_mag *data = w->data;
  track_filter_params params;
  gsl_vector_view v;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "track_filter: error opening %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  params.X = gsl_matrix_alloc(ntot, p);
  params.c = gsl_vector_alloc(p);
  params.rhs = gsl_vector_alloc(ntot);
  params.cov = gsl_matrix_alloc(p, p);

  v = gsl_matrix_row(params.X, 0);

  i = 1;
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: X data (nT)\n", i++);
  fprintf(fp, "# Field %zu: X internal model (nT)\n", i++);
  fprintf(fp, "# Field %zu: X external model (nT)\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);

      if (tptr->flags)
        continue;

      s = track_filter_proc(i, &params, data, w);
      if (s)
        continue;

      for (j = 0; j < tptr->n; ++j)
        {
          size_t didx = j + tptr->start_idx;
          double r = data->R + data->altitude[didx];
          double thetaq = M_PI / 2.0 - data->qdlat[didx] * M_PI / 180.0;
          double phi = data->longitude[didx] * M_PI / 180.0;
          double val_int, val_ext;

          /* compute internal model M(r,thetaq,phi) */
          gsl_vector_set_zero(&v.vector);
          track_filter_matrix_row_int(tptr->t_eq, r, thetaq, phi, &v.vector, w);
          gsl_blas_ddot(&v.vector, params.c, &val_int);

          /* compute external model K(r,thetaq,phi) */
          gsl_vector_set_zero(&v.vector);
          track_filter_matrix_row_ext(tptr->t_eq, r, thetaq, phi, &v.vector, w);
          gsl_blas_ddot(&v.vector, params.c, &val_ext);

          fprintf(fp, "%f %f %f %f\n",
                  data->qdlat[didx],
                  tptr->Bx[j],
                  val_int,
                  val_ext);

          tptr->Bx[j] -= val_int + val_ext;
        }

      fprintf(fp, "\n\n");
      fflush(fp);
    }

  gsl_vector_free(params.c);
  gsl_vector_free(params.rhs);
  gsl_matrix_free(params.X);
  gsl_matrix_free(params.cov);

  fclose(fp);

  return s;
} /* track_filter() */

static int
track_filter_proc(const size_t track_idx, track_filter_params *params,
                  satdata_mag *data, track_workspace *w)
{
  int s = 0;
  track_data *tptr = &(w->tracks[track_idx]);

  /* exclude points below this QD latitude, except a few to stabilize filter */
  const double exclude_lat = TRACK_FILTER_EXCLUDE_LAT;

  const size_t ntot = tptr->n;
  gsl_matrix *X = params->X;
  gsl_matrix *cov = params->cov;
  gsl_vector *c = params->c;
  gsl_vector *rhs = params->rhs;
  const size_t p = c->size;
  size_t n = 0;        /* number of data in LS system */
  size_t i;
  double qd_min, qd_max, X_min, X_max;
  int found_lower = 0, found_upper = 0;

  for (i = 1; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;
      double qdprev = data->qdlat[didx - 1];
      double qdcur = data->qdlat[didx];

      if (fabs(qdprev) > exclude_lat && fabs(qdcur) <= exclude_lat)
        {
          assert(found_lower == 0);

          found_lower = 1;
          qd_min = qdcur;
          X_min = tptr->Bx[i];
        }
      else if (fabs(qdprev) < exclude_lat && fabs(qdcur) >= exclude_lat)
        {
          assert(found_upper == 0);

          found_upper = 1;
          qd_max = qdcur;
          X_max = tptr->Bx[i];
        }
    }

  if (!found_lower || !found_upper)
    {
      fprintf(stderr, "track_filter_proc: unable to find linear interpolation points\n");
      return -1;
    }

  /* build LS matrix and RHS */
  for (i = 0; i < ntot; ++i)
    {
      size_t didx = i + tptr->start_idx;
      double r = data->R + data->altitude[didx];
      double thetaq = M_PI / 2.0 - data->qdlat[didx] * M_PI / 180.0;
      double phi = data->longitude[didx] * M_PI / 180.0;
      double qdlat = data->qdlat[didx];
      double rhsval = tptr->Bx[i];
      gsl_vector_view v;
      size_t nc;

      /* exclude PEJ region */
      if (fabs(qdlat) > 60.0)
        continue;

      /* exclude points in the EEJ region */
      if (fabs(qdlat) < exclude_lat)
        {
          /* keep every 10th point in equatorial region */
          if (i % 10 == 0)
            {
              /* interpolate linearly across EEJ region */
              rhsval = interp1d(qd_min, qd_max, X_min, X_max, qdlat);
            }
          else
            continue;
        }

      /* build current row of matrix */
      v = gsl_matrix_row(X, n);
      nc = track_filter_matrix_row(tptr->t_eq, r, thetaq, phi, &v.vector, w);
      assert(nc == p);

      /* set RHS value */
      gsl_vector_set(rhs, n, rhsval);

      ++n;
    }

  /* solve X c = rhs */
  {
    double h = 0.17;
    gsl_multifit_linear_workspace *multifit_p =
      gsl_multifit_linear_alloc(n, p);
    double chisq, cnorm;
    gsl_matrix_view Xv = gsl_matrix_submatrix(X, 0, 0, n, p);
    gsl_vector_view bv = gsl_vector_subvector(rhs, 0, n);

#if 0
    /* L-curve analysis */
    for (h = 1.0e-3; h < 1.0; h += 0.01)
      {
        /* perform damped least squares fit */
        s = gsl_multifit_linear_ridge(h, &Xv.matrix, &bv.vector, c,
                                      cov, &chisq, multifit_p);

        cnorm = gsl_blas_dnrm2(c);

        printf("%e %e %e\n",
               h,
               log(sqrt(chisq)),
               log(h * cnorm));
      }
    exit(1);
#elif 0
    /* perform damped least squares fit */
    s = gsl_multifit_linear_ridge(h, &Xv.matrix, &bv.vector, c,
                                  cov, &chisq, multifit_p);
#endif

    gsl_multifit_linear_free(multifit_p);
  }

  return s;
} /* track_filter_proc() */

static size_t
track_filter_matrix_row_int(const double t, const double r, const double theta,
                            const double phi, gsl_vector *v, track_workspace *w)
{
  size_t cidx = 0;
  msynth_workspace *msynth_p = w->msynth_workspace_p;

  /* compute basis functions for M(r,thetaq,phi) */
  msynth_green(r, theta, phi, msynth_p);

  /* internal coefficients */

  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(1, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(1, 1, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(1, -1, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(2, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(2, 1, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(2, -1, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(3, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(4, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(5, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(6, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(7, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(8, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(9, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(10, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(11, 0, msynth_p)]);
  gsl_vector_set(v, cidx++, msynth_p->dX[msynth_nmidx(12, 0, msynth_p)]);

  return cidx;
} /* track_filter_matrix_row_int() */

static size_t
track_filter_matrix_row_ext(const double t, const double r, const double theta,
                            const double phi, gsl_vector *v, track_workspace *w)
{
  const size_t int_offset = 16;
  size_t cidx = 0;
  msynth_workspace *msynth_p = w->msynth_workspace_p;
  const double lat = M_PI / 2.0 - theta;
  double phi_sm, theta_sm, lat_sm;
  size_t i;
  time_t unix_time = satdata_epoch2timet(t);
  double fday = time2fday(unix_time);

  trans(GEO2SM, fday, phi, lat, &phi_sm, &lat_sm);
  theta_sm = M_PI / 2.0 - lat_sm;

  msynth_green_ext(r, theta_sm, phi_sm, msynth_p);

  for (i = 0; i < msynth_p->nnm_ext; ++i)
    {
      double phi_ss, lat_ss;
      double X, Y, Z;

      trans_vec(SM2GEO, fday, phi_sm, lat_sm, msynth_p->dX_ext[i],
                msynth_p->dY_ext[i], msynth_p->dZ_ext[i], &phi_ss,
                &lat_ss, &X, &Y, &Z);

      msynth_p->dX_ext[i] = X;
      msynth_p->dY_ext[i] = Y;
      msynth_p->dZ_ext[i] = Z;
    }

  gsl_vector_set(v, int_offset + cidx++, msynth_p->dX_ext[msynth_nmidx(1, 0, msynth_p)]);
  gsl_vector_set(v, int_offset + cidx++, msynth_p->dX_ext[msynth_nmidx(1, 1, msynth_p)]);
  gsl_vector_set(v, int_offset + cidx++, msynth_p->dX_ext[msynth_nmidx(2, 0, msynth_p)]);
  gsl_vector_set(v, int_offset + cidx++, msynth_p->dX_ext[msynth_nmidx(2, 1, msynth_p)]);

  return cidx;
} /* track_filter_matrix_row_ext() */

static size_t
track_filter_matrix_row(const double t, const double r, const double theta,
                        const double phi, gsl_vector *v, track_workspace *w)
{
  size_t nc = 0;

  nc += track_filter_matrix_row_int(t, r, theta, phi, v, w);
  nc += track_filter_matrix_row_ext(t, r, theta, phi, v, w);

  return nc;
} /* track_filter_matrix_row() */

#if 0
/*
track_filter()
  Filter vector residuals with exponential moving
average filter

Inputs: tau - filter time constant (s)
        w   - workspace
*/

int
track_filter(const double tau, track_workspace *w)
{
  int s = 0;
  size_t i, j;
  const double dt = 1.0; /* sample spacing in s */
  const double alpha = dt / (tau + dt);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      double prevX, prevY, prevZ;

      /* filter X residual */
      prevX = tptr->Bx[0];
      prevY = tptr->By[0];
      prevZ = tptr->Bz[0];
      for (j = 1; j < tptr->n; ++j)
        {
          tptr->Bx[j] = alpha * tptr->Bx[j] + (1.0 - alpha) * prevX;
          tptr->By[j] = alpha * tptr->By[j] + (1.0 - alpha) * prevY;
          tptr->Bz[j] = alpha * tptr->Bz[j] + (1.0 - alpha) * prevZ;
          prevX = tptr->Bx[j];
          prevY = tptr->By[j];
          prevZ = tptr->Bz[j];
        }
    }

  return s;
} /* track_filter() */
#endif
