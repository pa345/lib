/*
 * mag_eej.c
 *
 * Invert scalar magnetic latitude profile for EEJ current height-integrated
 * density
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>

#include <mainlib/ml_bsearch.h>
#include <mainlib/ml_coord.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_geo.h>
#include <mainlib/ml_apex.h>
#include <mainlib/ml_mageq.h>

#include "geod2geoc.h"
#include "mag.h"

static int mag_eej_matrix_row(const double r, const double theta, const double phi,
                              const double b[3], gsl_vector *v, const mag_eej_workspace *w);
static int eej_add_B(const double ri[3], const double rjk[3],
                     const double Jjk[3], double Bij[3]);
static int eej_ls(const gsl_matrix *X, const gsl_vector *y,
                  gsl_vector *c, gsl_matrix *cov, mag_eej_workspace *w);

/*
mag_eej_alloc()
  Allocate workspace for fitting line currents to magnetic data

Inputs: year      - year of data to be inverted (for computing QD transformations)
        ncurr     - number of line currents (placed symmetrically around magnetic equator)
        altitude  - altitude of line currents (km)
        qdlat_max - line currents will span [-qdlat_max,qdlat_max]

Return: workspace
*/

mag_eej_workspace *
mag_eej_alloc(const int year, const size_t ncurr,
              const double altitude, const double qdlat_max)
{
  mag_eej_workspace *w;
  const size_t n = 5000;   /* maximum number of F^(2) data to fit */
  const double qdlon_min = 0.0;
  const double qdlon_max = 359.0;
  const size_t nlon = 360;
  const double dlon = (qdlon_max - qdlon_min) / (nlon - 1.0);
  const double dlat = (2.0 * qdlat_max) / (ncurr - 1.0);
  const double epoch = (double) year;
  size_t j, k;

  w = calloc(1, sizeof(mag_eej_workspace));
  if (!w)
    return 0;

  w->qd_alt = altitude;
  w->nlon = nlon;
  w->p = ncurr;
  w->qdlat_max = qdlat_max;
  w->dqdlat = dlat;

  /* distance in km between arc currents */
  w->curr_dist_km = dlat * (40008.0) / 360.0;

  w->apex_workspace_p = apex_alloc();
  w->mageq_workspace_p = mageq_alloc();

  w->Jx = gsl_matrix_alloc(ncurr, nlon);
  w->Jy = gsl_matrix_alloc(ncurr, nlon);
  w->Jz = gsl_matrix_alloc(ncurr, nlon);
  w->mid_pos_x = gsl_matrix_alloc(ncurr, nlon);
  w->mid_pos_y = gsl_matrix_alloc(ncurr, nlon);
  w->mid_pos_z = gsl_matrix_alloc(ncurr, nlon);
  w->mid_lon = malloc(nlon * sizeof(double));

  w->X = gsl_matrix_alloc(n, w->p);
  w->cov = gsl_matrix_alloc(w->p, w->p);
  w->S = gsl_vector_alloc(w->p);
  w->rhs = gsl_vector_alloc(n);
  w->Xs = gsl_matrix_alloc(n, w->p);
  w->ys = gsl_vector_alloc(n);
  w->cs = gsl_vector_alloc(w->p);
  w->M = gsl_matrix_alloc(n, w->p);

  w->multifit_p = gsl_multifit_linear_alloc(n, w->p);
  w->nreg = 200;
  w->rho = gsl_vector_alloc(w->nreg);
  w->eta = gsl_vector_alloc(w->nreg);
  w->reg_param = gsl_vector_alloc(w->nreg);

  /* construct smoothing matrix L */
#if 1
  {
    const size_t k = 2;
    w->L = gsl_matrix_alloc(w->p - k, w->p);
    w->Ltau = gsl_vector_alloc(w->p - k);

    gsl_multifit_linear_Lk(w->p, k, w->L);
    gsl_multifit_linear_L_decomp(w->L, w->Ltau);
  }
#else
  w->L = gsl_matrix_alloc(w->p, w->p);

  {
    const size_t kmax = 2;
    gsl_vector *alpha = gsl_vector_alloc(kmax + 1);
    gsl_vector_set_all(alpha, 1.0);
    gsl_vector_set(alpha, 2, 10.0);

    gsl_multifit_linear_Lsobolev(w->p, kmax, alpha, w->L, w->multifit_p);

    gsl_vector_free(alpha);
  }
#endif

  /*
   * precompute Cartesian positions of each current segment for each
   * line current arc; current arcs follow lines of constant QD
   * latitude and segments are equally spaced in QD longitude. These
   * QD grid points are converted to spherical geocentric points and then
   * to Cartesian.
   *
   * XXX - 24 Jan 2014 - the apex_transform_inv() calls here produce
   * geographic latitudes which aren't as "smooth" as Stefan/GFZ modified
   * apxntrp.f calls - they tend to cluster closely together for
   * certain qdlat inputs as can be seen by plotting
   * (mid_lon(k), lat2-lat1) in the below loop
   */
  for (j = 0; j < ncurr; ++j)
    {
      double qdlat = -qdlat_max + j * dlat;

      for (k = 0; k < nlon; ++k)
        {
          double qdlon = qdlon_min + (k + 0.5) * dlon; /* midpoint lon */
          double qdlon1 = qdlon_min + k * dlon;        /* endpoint lons */
          double qdlon2 = qdlon_min + (k + 1.0) * dlon;
          double glat, glon; /* geodetic lat/lon */
          double glat1, glat2, glon1, glon2;
          double r1, r2, lat1, lat2;
          double lat, r;     /* geographic spherical lat/radius */
          double delta, az, delta_km;
          double Jt, Jp;     /* theta and phi components of J */
          double Jx, Jy, Jz; /* ECEF Cartesian current components */
          double mid_pos[3]; /* ECEF Cartesian position of segment midpt */

          /* compute geodetic lon/lat of qd point */
          apex_transform_inv_geodetic(epoch, qdlat, qdlon, w->qd_alt,
                                      &glat, &glon, w->apex_workspace_p);

          /* convert geodetic to geocentric spherical */
          geod2geoc(glat, w->qd_alt, &lat, &r);

          /* sanity check */
          {
            double tmplat, tmpalt;
            geo2geodetic(lat, glon, r, &tmplat, &tmpalt);
            gsl_test_rel(tmplat, glat, 1.0e-7, "latitude");
            gsl_test_rel(tmpalt, w->qd_alt, 1.0e-6, "altitude");
          }

          /*
           * now (r,lat,glon) contain the geocentric spherical coordinate
           * of the midpoint of segment k of current j; save ECEF
           * Cartesian position of the midpoints
           */
          sphere2cart(glon, lat, r, mid_pos);
          gsl_matrix_set(w->mid_pos_x, j, k, mid_pos[0]);
          gsl_matrix_set(w->mid_pos_y, j, k, mid_pos[1]);
          gsl_matrix_set(w->mid_pos_z, j, k, mid_pos[2]);

          /*
           * save segment longitude to check later if its within 30 deg;
           * this array is overwritten for each current j so its
           * an approximate position
           */
          w->mid_lon[k] = glon * 180.0 / M_PI;

          /* find geocentric coordinates of segment endpoints */
          apex_transform_inv_geodetic(epoch, qdlat, qdlon1, w->qd_alt,
                                      &glat1, &glon1, w->apex_workspace_p);
          apex_transform_inv_geodetic(epoch, qdlat, qdlon2, w->qd_alt,
                                      &glat2, &glon2, w->apex_workspace_p);

          /* convert geodetic to geocentric spherical */
          geod2geoc(glat1, w->qd_alt, &lat1, &r1);
          geod2geoc(glat2, w->qd_alt, &lat2, &r2);

          lat1 *= 180.0 / M_PI;
          glon1 *= 180.0 / M_PI;
          lat2 *= 180.0 / M_PI;
          glon2 *= 180.0 / M_PI;

          my_delaz(lat1, glon1, lat2, glon2, &delta, &az);

          delta_km = delta * M_PI / 180.0 * r;

          /* compute current components */
          Jt = -delta_km * cos(az * M_PI / 180.0);
          Jp = delta_km * sin(az * M_PI / 180.0);

          sphere2cart_vec(glon1*M_PI/180.0, lat1*M_PI/180.0,
                          -Jt, Jp, 0.0, &Jx, &Jy, &Jz);

          /* store current vector "green function" for later use */
          gsl_matrix_set(w->Jx, j, k, Jx);
          gsl_matrix_set(w->Jy, j, k, Jy);
          gsl_matrix_set(w->Jz, j, k, Jz);
        }
    }

  return w;
} /* mag_eej_alloc() */

void
mag_eej_free(mag_eej_workspace *w)
{
  if (w->apex_workspace_p)
    apex_free(w->apex_workspace_p);

  if (w->Jx)
    gsl_matrix_free(w->Jx);

  if (w->Jy)
    gsl_matrix_free(w->Jy);

  if (w->Jz)
    gsl_matrix_free(w->Jz);

  if (w->mid_pos_x)
    gsl_matrix_free(w->mid_pos_x);

  if (w->mid_pos_y)
    gsl_matrix_free(w->mid_pos_y);

  if (w->mid_pos_z)
    gsl_matrix_free(w->mid_pos_z);

  if (w->mid_lon)
    free(w->mid_lon);

  if (w->X)
    gsl_matrix_free(w->X);

  if (w->cov)
    gsl_matrix_free(w->cov);

  if (w->S)
    gsl_vector_free(w->S);

  if (w->rhs)
    gsl_vector_free(w->rhs);

  if (w->L)
    gsl_matrix_free(w->L);

  if (w->Ltau)
    gsl_vector_free(w->Ltau);

  if (w->Xs)
    gsl_matrix_free(w->Xs);

  if (w->ys)
    gsl_vector_free(w->ys);

  if (w->cs)
    gsl_vector_free(w->cs);

  if (w->M)
    gsl_matrix_free(w->M);

  if (w->rho)
    gsl_vector_free(w->rho);

  if (w->eta)
    gsl_vector_free(w->eta);

  if (w->reg_param)
    gsl_vector_free(w->reg_param);

  if (w->mageq_workspace_p)
    mageq_free(w->mageq_workspace_p);

  if (w->multifit_p)
    gsl_multifit_linear_free(w->multifit_p);

  free(w);
}

/*
mag_eej_proc()
  Invert F^(2) residual for EEJ current density

Inputs: track - satellite track
        J     - (output) where to store EEJ height integrated
                current density in A/m (array of size ncurr)
        w     - workspace

Notes:
1) On output, track->F2_fit is updated to contain the fit from the line
current model
*/

int
mag_eej_proc(mag_track *track, double *J, mag_eej_workspace *w)
{
  int s = 0;
  size_t i;
  size_t ndata = 0;
  gsl_matrix_view Xv;
  gsl_vector_view v;

  /* build LS matrix and rhs vector */
  for (i = 0; i < track->n; ++i)
    {
      double r = track->r[i];
      double theta = track->theta[i];
      double phi = track->phi[i];
      gsl_vector_view v;
      double bi[3]; /* unit field vector in NEC */

      /* ignore data outside [-qd_max,qd_max] */
      if (fabs(track->qdlat[i]) > w->qdlat_max)
        continue;

      /* compute unit magnetic field vector at this point */
      bi[0] = track->Bx_int[i] / track->F_int[i];
      bi[1] = track->By_int[i] / track->F_int[i];
      bi[2] = track->Bz_int[i] / track->F_int[i];

      /* add F^(2) to rhs vector */
      gsl_vector_set(w->rhs, ndata, track->F2[i]);

      /* construct this row of LS matrix */
      v = gsl_matrix_row(w->X, ndata);
      mag_eej_matrix_row(r, theta, phi, bi, &v.vector, w);

      ++ndata;
    }

  /* perform least squares fit */
  Xv = gsl_matrix_submatrix(w->X, 0, 0, ndata, w->p);
  v = gsl_vector_subvector(w->rhs, 0, ndata);
  s = eej_ls(&Xv.matrix, &v.vector, w->S, w->cov, w);
  if (s)
    return s;

  fprintf(stderr, "mag_eej_proc: final regularization factor = %g\n",
          gsl_vector_get(w->reg_param, w->reg_idx));
  fprintf(stderr, "mag_eej_proc: residual norm = %e\n", w->rnorm);
  fprintf(stderr, "mag_eej_proc: solution norm = %e\n", w->snorm);

  /* compute output current density */
  v = gsl_vector_view_array(J, w->p);
  gsl_vector_memcpy(&v.vector, w->S);
  gsl_vector_scale(&v.vector, 1.0 / w->curr_dist_km);

  /* compute and store F^(2) fit */
  for (i = 0; i < track->n; ++i)
    {
      double r = track->r[i];
      double theta = track->theta[i];
      double phi = track->phi[i];
      double bi[3]; /* unit field vector in NEC */

      /* compute unit magnetic field vector at this point */
      bi[0] = track->Bx_int[i] / track->F_int[i];
      bi[1] = track->By_int[i] / track->F_int[i];
      bi[2] = track->Bz_int[i] / track->F_int[i];

      /* construct this row of LS matrix */
      v = gsl_matrix_row(w->X, 0);
      mag_eej_matrix_row(r, theta, phi, bi, &v.vector, w);

      gsl_blas_ddot(&v.vector, w->S, &(track->F2_fit[i]));
    }

  return s;
} /* mag_eej() */

/*
mag_eej_vector_proc()
  Invert F^(2) residual for EEJ current density

Inputs: track - satellite track
        J     - (output) where to store EEJ height integrated
                current density in A/m (array of size ncurr)
        w     - workspace

Notes:
1) On output, track->F2_fit is updated to contain the fit from the line
current model
*/

int
mag_eej_vector_proc(mag_track *track, double *J, mag_eej_workspace *w)
{
  int s = 0;
  size_t i;
  magfit_workspace *magfit_p;
  magfit_parameters magfit_params = magfit_default_parameters();

  magfit_params.lat_min = -40.0;
  magfit_params.lat_max = 40.0;
  magfit_params.lat_spacing1d = 0.25;
#if 1
  magfit_params.flags = MAGFIT_FLG_FIT_Z;
#else
  magfit_params.flags = MAGFIT_FLG_FIT_X | MAGFIT_FLG_FIT_Z;
#endif
  magfit_params.flags |= MAGFIT_FLG_SECS_FIT_DF;
  magfit_p = magfit_alloc(magfit_secs1d, &magfit_params);

  /* build LS matrix and rhs vector */
  for (i = 0; i < track->n; ++i)
    {
      double B[3];

      /* ignore data outside [-qd_max,qd_max] */
      if (fabs(track->qdlat[i]) > w->qdlat_max)
        continue;

      B[0] = track->X2[i];
      B[1] = track->Y2[i];
      B[2] = track->Z2[i];

      magfit_add_datum(track->t[i], track->r[i], track->theta[i], track->phi[i],
                       track->qdlat[i], B, magfit_p);
    }

  /* fit EEJ model */
  magfit_fit(&(w->rnorm), &(w->snorm), magfit_p);

  fprintf(stderr, "mag_eej_vector_proc: final regularization factor = %g\n",
          gsl_vector_get(w->reg_param, w->reg_idx));
  fprintf(stderr, "mag_eej_vector_proc: residual norm = %e\n", w->rnorm);
  fprintf(stderr, "mag_eej_vector_proc: solution norm = %e\n", w->snorm);

  /* fill output current vector */
  {
    double epoch = satdata_epoch2year(track->t_eq);
    double r = R_EARTH_KM + 110.0;
    double qdlat, alon, alat;

    apex_transform(epoch, track->theta_eq, track->phi_eq, r, &alon, &alat,
                   &qdlat, NULL, NULL, NULL, w->apex_workspace_p);

    for (i = 0; i < w->p; ++i)
      {
        double K[3];
        double qdlat = -w->qdlat_max + i * w->dqdlat;
        double glat, glon, theta;

        apex_transform_inv(epoch, qdlat, alon, 110.0, &glat, &glon, w->apex_workspace_p);
        theta = M_PI / 2.0 - glat;

        magfit_eval_J(R_EARTH_KM + 110.0, theta, 0.0, K, magfit_p);
        J[i] = K[1] * 1.0e-3; /* convert to A/m */
      }
  }

  /* compute and store F^(2) fit */
  for (i = 0; i < track->n; ++i)
    {
      double B[3];

      magfit_eval_B(track->t[i], track->r[i], track->theta[i], track->phi[i],
                    B, magfit_p);

      track->X2_fit[i] = B[0];
      track->Y2_fit[i] = B[1];
      track->Z2_fit[i] = B[2];
    }

  magfit_free(magfit_p);

  return s;
}

static int
mag_eej_matrix_row(const double r, const double theta, const double phi,
                   const double b[3], gsl_vector *v, const mag_eej_workspace *w)
{
  const double max_dlon = 30.0; /* maximum longitude window in deg */
  const double lon_deg = phi * 180.0 / M_PI;
  const double lat_rad = M_PI / 2.0 - theta;
  size_t j, k;
  double ri[3];

  /* compute ECEF cartesian position vector of satellite position */
  sphere2cart(phi, lat_rad, r, ri);

  for (j = 0; j < w->p; ++j)
    {
      double Bij[3]; /* magnetic field at ri due to line j in ECEF */
      double Bij_NEC[3]; /* Bij in NEC */
      double Fij; /* B_{ij} . b_i */

      Bij[0] = Bij[1] = Bij[2] = 0.0;

      /* compute B_{ij} = sum_k dB_{ijk} (eq 14) */
      for (k = 0; k < w->nlon; ++k)
        {
          double rjk[3], Jjk[3];
          double dlon = wrap180(w->mid_lon[k] - lon_deg);

          /* check if this segment is within max_dlon of obs point */
          if (fabs(dlon) > max_dlon)
            continue;

          rjk[0] = gsl_matrix_get(w->mid_pos_x, j, k);
          rjk[1] = gsl_matrix_get(w->mid_pos_y, j, k);
          rjk[2] = gsl_matrix_get(w->mid_pos_z, j, k);

          Jjk[0] = gsl_matrix_get(w->Jx, j, k);
          Jjk[1] = gsl_matrix_get(w->Jy, j, k);
          Jjk[2] = gsl_matrix_get(w->Jz, j, k);

          eej_add_B(ri, rjk, Jjk, Bij);
        }

      /* convert B_{ij} to NEC coordinates */
      cart2sphere_vec(phi, lat_rad, Bij[0], Bij[1], Bij[2],
                      &Bij_NEC[0], &Bij_NEC[1], &Bij_NEC[2]);

      /* compute F_{ij} = B_{ij} . b_i (eq 15) */
      Fij = vec_dot(Bij_NEC, b);

      /* add to least squares matrix */
      gsl_vector_set(v, j, Fij);
    } /* for (j = 0; j < w->p; ++j) */

  return 0;
} /* mag_eej_matrix_row() */

/*
eej_add_B()
  Add contribution of a single line current segment to a running
vector sum using Biot-Savart law. See eq 13 of paper

Inputs: ri  - observer location (ECEF)
        rjk - midpoint of longitudinal segment k of arc current j (ECEF)
        Jjk - unit current vector of flow direction for this segment (ECEF)
        Bij - (output) updated with ECEF X,Y,Z components of magnetic field
              for this line segment at observer location ri (nT)

Notes:
1) on output, Bij = Bij + dBijk
*/

static int
eej_add_B(const double ri[3], const double rjk[3], const double Jjk[3],
          double Bij[3])
{
  /*
   * the factor of 1e9 gives Bij units of nT if current density Jjk
   * given in km
   */
  const double mu0 = (4.0 * M_PI) * 1.0e-7;
  const double C = 1.0e9 * mu0 / (4.0 * M_PI);

  double dx = ri[0] - rjk[0];
  double dy = ri[1] - rjk[1];
  double dz = ri[2] - rjk[2];
  double invr3 = pow(dx*dx + dy*dy + dz*dz, -1.5);

  Bij[0] += C * (Jjk[1]*dz - Jjk[2]*dy) * invr3;
  Bij[1] += C * (Jjk[2]*dx - Jjk[0]*dz) * invr3;
  Bij[2] += C * (Jjk[0]*dy - Jjk[1]*dx) * invr3;

  return GSL_SUCCESS;
} /* eej_add_B() */

/*
eej_ls()
  Solve current inversion least squares problem with ridge
regression
*/

static int
eej_ls(const gsl_matrix *X, const gsl_vector *y, gsl_vector *c,
       gsl_matrix *cov, mag_eej_workspace *w)
{
  int s = 0;
  const size_t n = X->size1;
  const size_t p = X->size2;
  const size_t m = w->L->size1;
  const size_t npm = n - (p - m);
  gsl_matrix_view Xs = gsl_matrix_submatrix(w->Xs, 0, 0, npm, m);
  gsl_vector_view ys = gsl_vector_subvector(w->ys, 0, npm);
  gsl_vector_view cs = gsl_vector_subvector(w->cs, 0, m);
  gsl_matrix_view M = (m < p) ? gsl_matrix_submatrix(w->M, 0, 0, n, p) : gsl_matrix_submatrix(w->M, 0, 0, m, p);
  double lambda, smax;

  /* convert (X, y) to standard form (X~,y~) */
  s = gsl_multifit_linear_stdform2(w->L, w->Ltau, X, y, &Xs.matrix, &ys.vector, &M.matrix, w->multifit_p);
  if (s)
    return s;

  /* compute SVD of X~ */
  s = gsl_multifit_linear_svd(&Xs.matrix, w->multifit_p);
  if (s)
    return s;

  /* compute L-curve and its corner of LS system (X~,y~) */
  s = gsl_multifit_linear_lcurve(&ys.vector, w->reg_param, w->rho, w->eta, w->multifit_p);
  if (s)
    return s;

  s = gsl_multifit_linear_lcorner(w->rho, w->eta, &(w->reg_idx));
  if (s)
    return s;

  /*
   * sometimes the corner finding routine fails if there are local corners;
   * therefore set a lower bound on lambda
   */
  smax = gsl_vector_get(w->multifit_p->S, 0);

  lambda = GSL_MAX(gsl_vector_get(w->reg_param, w->reg_idx), 1.0e-3*smax);

  w->reg_idx = bsearch_double(w->reg_param->data, lambda, 0, w->nreg - 1);
  lambda = gsl_vector_get(w->reg_param, w->reg_idx);

  /* solve LS system with optimal lambda */
  gsl_multifit_linear_solve(lambda, &Xs.matrix, &ys.vector, &cs.vector, &(w->rnorm), &(w->snorm), w->multifit_p);

  /* backtransform c~ to recover c */
  gsl_multifit_linear_genform2(w->L, w->Ltau, X, y, &cs.vector, &M.matrix, c, w->multifit_p);

  return s;
} /* eej_ls() */
