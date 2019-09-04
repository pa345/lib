/*
 * swarmcdf.c
 *
 * This module contains routines to read/write Swarm EEF CDF data files
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cdf.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/time.h>

#include <mainlib/ml_satdata.h>

#include "swarmeef.h"

swarm_eef *
swarm_eef_alloc(size_t n)
{
  swarm_eef *data;

  data = calloc(1, sizeof(swarm_eef));
  if (!data)
    {
      fprintf(stderr, "swarm_eef_alloc: calloc failed: %s\n",
              strerror(errno));
      return 0;
    }

  swarm_eef_realloc(n, data);

  data->n = 0;

  return data;
} /* swarm_eef_alloc() */

swarm_eef *
swarm_eef_realloc(size_t n, swarm_eef *data)
{
  size_t i;

  if (!data)
    {
      return swarm_eef_alloc(n);
    }

  if (n <= data->ntot)
    return data; /* nothing to do */

  data->t = realloc(data->t, n * sizeof(double));
  data->latitude = realloc(data->latitude, n * sizeof(double));
  data->longitude = realloc(data->longitude, n * sizeof(double));
  data->EEF = realloc(data->EEF, n * sizeof(double));
  data->RelErr = realloc(data->RelErr, n * sizeof(double));
  data->Flags = realloc(data->Flags, n * sizeof(unsigned short));

  /* initialize newly allocated memory */
  for (i = data->ntot; i < n; ++i)
    {
      data->t[i] = 0.0;
      data->latitude[i] = 0.0;
      data->longitude[i] = 0.0;
      data->EEF[i] = 0.0;
      data->RelErr[i] = 0.0;
      data->Flags[i] = 0;
    }

  data->ntot = n;

  return data;
} /* swarm_eef_realloc() */

/*
swarm_eef_read()
  Read Swarm EEF CDF data file

Inputs: filename - CDF file
        data     - (input/output) where to store newly read data
                   If NULL, this data structure is allocated to store
                   new data. If not NULL, new data is appended to
                   pre-existing data

Return: pointer to data structure containing data
*/

swarm_eef *
swarm_eef_read(const char *filename, swarm_eef *data)
{
  CDFstatus s;
  CDFid input_id;
  long nrec;      /* number of records in file */
  long epoch_num;
  size_t cur_idx;
  size_t ntot, ncur;

  s = CDFopen(filename, &input_id);
  if (s != CDF_OK)
    {
      fprintf(stderr, "swarm_eef_read: error opening file %s: %s\n",
              filename, cdf_error(s));
      return 0;
    }

  epoch_num = CDFgetVarNum(input_id, "timestamp");
  if (epoch_num < 0)
    {
      fprintf(stderr, "swarm_eef_read: unable to read timestamp\n");
      return 0;
    }

  s = CDFgetzVarMaxWrittenRecNum(input_id, epoch_num, &nrec);
  if (s != CDF_OK)
    {
      fprintf(stderr, "swarm_eef_read: unable to read records written\n");
      return 0;
    }

  /* since records start at 0 */
  ++nrec;

  /* allocate 30 days of data at a time */
  ntot = data ? data->ntot : 0;
  ncur = data ? data->n : 0;
  if (ntot < ncur + nrec)
    data = swarm_eef_realloc(30 * 86400 + ntot, data);

  if (!data)
    return 0;

  /* now read the data into arrays */
  cur_idx = data->n;

  s += CDFgetVarAllRecordsByVarName(input_id, "timestamp", &(data->t[cur_idx]));
  s += CDFgetVarAllRecordsByVarName(input_id, "longitude", &(data->longitude[cur_idx]));
  s += CDFgetVarAllRecordsByVarName(input_id, "latitude", &(data->latitude[cur_idx]));
  s += CDFgetVarAllRecordsByVarName(input_id, "EEF", &(data->EEF[cur_idx]));
  s += CDFgetVarAllRecordsByVarName(input_id, "RelErr", &(data->RelErr[cur_idx]));
  s += CDFgetVarAllRecordsByVarName(input_id, "flags", &(data->Flags[cur_idx]));

  if (s != CDF_OK)
    {
      fprintf(stderr, "swarm_eef_read: error reading records\n");
    }

  data->n += nrec;

  CDFclose(input_id);

  return data;
} /* swarm_eef_read() */

swarm_eef *
swarm_eef_read_idx(const char *filename, const int verbose)
{
  swarm_eef *data = NULL;
  FILE *fp;
  char buffer[SWARMEEF_MAX_BUF];

  if (satdata_iscdf(filename))
    {
      data = swarm_eef_read(filename, NULL);
      return data;
    }

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "swarm_eef_read_idx: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  /* read in data files */
  while (fgets(buffer, SWARMEEF_MAX_BUF, fp) != 0)
    {
      struct timeval tv0, tv1;
      double dt;

      buffer[strlen(buffer) - 1] = '\0';

      if (verbose)
        {
          fprintf(stderr, "Reading %s...", buffer);
          gettimeofday(&tv0, NULL);
        }

      data = swarm_eef_read(buffer, data);
      if (!data)
        break;

      if (verbose)
        {
          gettimeofday(&tv1, NULL);
          dt = (tv1.tv_sec - tv0.tv_sec) + (tv1.tv_usec - tv0.tv_usec) * 1.0e-6;
          fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n, dt);
        }
    }

  fclose(fp);

  return data;
} /* swarm_eef_read_idx() */
