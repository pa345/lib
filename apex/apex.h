/*
 * apex.h
 * Patrick Alken
 */

#ifndef INCLUDED_apex_h
#define INCLUDED_apex_h

#define APEX_NLAT    36
#define APEX_NLON    72
#define APEX_NALT    300

#define APEX_DATAFILE "/data/palken/lib/apex/apexsh.dat"

typedef struct
{
  size_t ntheta;
  size_t nphi;
  size_t nalt;

  int nepochs;
  float *epochgrid;

  int lwk;
  float *wk;
} apex_workspace;

/*
 * Prototypes
 */

apex_workspace *apex_alloc(int year);
void apex_free(apex_workspace *w);
int apex_makefile(const char *filename, apex_workspace *w);
int apex_readfile(const char *filename, int year, apex_workspace *w);
int apex_transform(double theta, double phi, double r,
                   double *apexlon, double *apexlat, double *qdlat,
                   double *E1, double *E2, double *E3, apex_workspace *w);
int apex_transform_geodetic(double theta, double phi, double alt,
                            double *apexlon, double *apexlat, double *qdlat,
                            double *E1, double *E2, double *E3, apex_workspace *w);
int apex_transform_inv(double qdlat, double qdlon, double alt,
                       double *glat, double *glon,
                       apex_workspace *w);
int apex_transform_inv_geodetic(double qdlat, double qdlon, double alt,
                                double *glat, double *glon,
                                apex_workspace *w);

#endif /* INCLUDED_apex_h */
