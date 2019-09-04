/*
 * tiegcm3d_print.c
 *
 * This module contains routines to print tiegcm data
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include <netcdf.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_interp.h>

#include "tiegcm3d.h"

/*
tiegcm3d_print_time()
  Print time series of Jr/Jt/Jp for a fixed location

Inputs: filename - output file
        data     - tiegcm data
        ir       - r grid point
        ilat     - latitude grid point
        ilon     - longitude grid point
*/

int
tiegcm3d_print_time(const char *filename, const tiegcm3d_data *data, const int ir, const int ilat, const int ilon)
{
  int s = 0;
  size_t it;
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "tiegcm3d_print_time: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Latitude: %.2f (deg)\n", data->glat[ilat]);
  fprintf(fp, "# Longitude: %.2f (deg)\n", data->glon[ilon]);
  fprintf(fp, "# Radius: %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - R_EARTH_KM);
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp, "# Field %zu: J_r (uA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: J_t (uA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: J_p (uA/m^2)\n", i++);

  for (it = 0; it < data->nt; ++it)
    {
      size_t idx = TIEGCM3D_IDX(it, ir, ilat, ilon, data);

      fprintf(fp, "%ld %16.4e %16.4e %16.4e\n",
              data->t[it],
              data->Jr[idx],
              data->Jt[idx],
              data->Jp[idx]);
    }

  fclose(fp);

  return s;
}
