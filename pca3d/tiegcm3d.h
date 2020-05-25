/*
 * tiegcm3d.h
 */

#ifndef INCLUDED_tiegcm3d_h
#define INCLUDED_tiegcm3d_h

#include <time.h>
#include <netcdf.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>

/* time domain grid index function */
#define TIEGCM3D_IDX(it,ir,ilat,ilon,data)                    (CIDX4((it),(data->nt),(ir),(data->nr),(ilon),(data->nlon),(ilat),(data->nlat)))

/* frequency domain grid index function */
#define TIEGCM3D_FREQIDX(it,ifreq,ir,ilat,ilon,data,T,nfreq)  (CIDX5((it),(T),(ifreq),(nfreq),(ir),(data->nr),(ilon),(data->nlon),(ilat),(data->nlat)))

typedef struct
{
  time_t *t;    /* timestamps, size nt */
  double *year; /* year, size nt */
  double *doy;  /* day of year, size nt */
  double *ut;   /* UT (hours), size nt */
  double *r;    /* radius (m), size nr */
  double *glon; /* geodetic longitude (degrees), size nlon */
  double *glat; /* geodetic latitude (degrees), size nlat */

  double *Jr;   /* J_r grid, nt-by-nr-by-nlon-by-nlat (A/m^2) */
  double *Jt;   /* J_t grid, nt-by-nr-by-nlon-by-nlat (A/m^2) */
  double *Jp;   /* J_p grid, nt-by-nr-by-nlon-by-nlat (A/m^2) */

  double *work; /* temporary grid, nt-by-nr-nlon-by-nlat */

  size_t nt_max; /* maximum time snapshots in 1 TIEGCM file */

  size_t nt;    /* number of time snapshots stored */
  size_t nr;    /* length of 'r' */
  size_t nlon;  /* length of 'glon' */
  size_t nlat;  /* length of 'glat' */
  size_t ntot;  /* total data allocated */
} tiegcm3d_data;

/* FFT data based on J grids */
typedef struct
{
  size_t nt;       /* number of time steps */
  size_t nfreq;    /* number of frequencies */
  size_t nr;       /* number of r grid points */
  size_t nlon;     /* number of phi grid points */
  size_t nlat;     /* number of theta grid points */
  size_t T;        /* number of time window segments */
  double fs;       /* sampling frequency (samples/day) */
  double window_size;  /* window size in days */
  double window_shift; /* window shift in days */
  size_t nwindow;  /* number of samples in each window */
  gsl_vector *window; /* window function, size nwindow */
  time_t *t;       /* timestamps, size nt */
  double *r;       /* radius (km), size nr */
  double *glon;    /* geodetic longitude (degrees), size nlon */
  double *glat;    /* geodetic latitude (degrees), size nlat */
  double *Jr;      /* J_r grid, nt-by-nr-by-nlon-by-nlat */
  double *Jt;      /* J_t grid, nt-by-nr-by-nlon-by-nlat */
  double *Jp;      /* J_p grid, nt-by-nr-by-nlon-by-nlat */
  gsl_complex *Qr; /* J_r transform grid, nfreq-by-nr-by-nlat-by-nlon */
  gsl_complex *Qt; /* J_theta transform grid, nfreq-by-nr-by-nlat-by-nlon */
  gsl_complex *Qp; /* J_phi transform grid, nfreq-by-nr-by-nlat-by-nlon */
} tiegcm3d_fft_data;

/* tiegcm3d_alloc.c */
tiegcm3d_data *tiegcm3d_alloc(const size_t nt, const size_t nr, const size_t nlon, const size_t nlat);
tiegcm3d_data *tiegcm3d_realloc(const size_t nt, const size_t nr, const size_t nlon, const size_t nlat,
                                tiegcm3d_data *data);
void tiegcm3d_free(tiegcm3d_data *data);

/* tiegcm3d_print.c */
int tiegcm3d_print_time(const char *filename, const tiegcm3d_data *data, const int ir, const int ilat, const int ilon);

/* tiegcm3d_read.c */
tiegcm3d_data *tiegcm3d_read(const char *filename, tiegcm3d_data *data);
tiegcm3d_data *tiegcm3d_read_B(const char *filename, tiegcm3d_data *data);

#endif /* INCLUDED_tiegcm3d_h */
