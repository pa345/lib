/*
 * invert_tmode.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

#include <common/common.h>

#include "invert_tmode.h"

/*
invert_tmode_alloc()
  Allocate invert_tmode_workspace

Inputs: nfreq  - number of frequency bands
        nmodes - number of modes per frequency band, length nfreq
        N      - length of time series for each mode

Return: pointer to workspace
*/

invert_tmode_workspace *
invert_tmode_alloc(const size_t nfreq, const size_t nmodes[], const size_t N)
{
  invert_tmode_workspace *w;
  size_t i;

  w = calloc(1, sizeof(invert_tmode_workspace));
  if (!w)
    return 0;

  w->nfreq = nfreq;
  w->N = N;

  w->nmodes = malloc(nfreq * sizeof(size_t));
  if (!w->nmodes)
    {
      invert_tmode_free(w);
      return 0;
    }

  w->modes = calloc(nfreq, sizeof(gsl_matrix_complex *));
  if (!w->modes)
    {
      invert_tmode_free(w);
      return 0;
    }

  for (i = 0; i < nfreq; ++i)
    {
      w->nmodes[i] = nmodes[i];
      w->modes[i] = gsl_matrix_complex_alloc(N, nmodes[i]);
      if (!w->modes[i])
        {
          invert_tmode_free(w);
          fprintf(stderr, "invert_tmode_alloc: failed to allocate matrix\n");
          return 0;
        }
    }

  return w;
}

void
invert_tmode_free(invert_tmode_workspace *w)
{
  if (w->nmodes)
    free(w->nmodes);

  if (w->modes)
    {
      size_t i;

      for (i = 0; i < w->nfreq; ++i)
        {
          if (w->modes[i])
            gsl_matrix_complex_free(w->modes[i]);
        }

      free(w->modes);
    }

  free(w);
}

/*
invert_tmode_read_ascii()
  Read ASCII file of temporal mode data

Inputs: filename - filename to read
        freq     - frequency band, in [0,nfreq-1]
        mode     - mode number, in [0,nmodes[freq]-1]
        w        - workspace
*/

int
invert_tmode_read_ascii(const char * filename, const size_t freq, const size_t mode, invert_tmode_workspace * w)
{
  if (freq >= w->nfreq)
    {
      GSL_ERROR ("invalid frequency band", GSL_EBADLEN);
    }
  else if (mode >= w->nmodes[freq])
    {
      GSL_ERROR ("invalid mode number", GSL_EBADLEN);
    }
  else
    {
      FILE *fp;
      gsl_vector_complex_view v = gsl_matrix_complex_column(w->modes[freq], mode);
      char buffer[2048];
      size_t n = 0;

      fp = fopen(filename, "r");
      if (!fp)
        {
          fprintf(stderr, "invert_tmode_read_ascii: unable to open %s: %s\n",
                  filename, strerror(errno));
          return -1;
        }

      while (fgets(buffer, sizeof(buffer), fp) != NULL)
        {
          int c;
          time_t t;
          double re, im;

          c = sscanf(buffer, "%ld %lf %lf", &t, &re, &im);
          if (c != 3)
            continue;

          gsl_vector_complex_set(&v.vector, n++, gsl_complex_rect(re, im));
        }

      if (n != w->N)
        {
          GSL_ERROR ("mode file has wrong length", GSL_EBADLEN);
        }

      fclose(fp);

      return GSL_SUCCESS;
    }
}

int
invert_tmode_write_binary(const char * filename, invert_tmode_workspace * w)
{
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "invert_tmode_write_binary: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  fwrite(&(w->nfreq), sizeof(size_t), 1, fp);
  fwrite(w->nmodes, sizeof(size_t), w->nfreq, fp);
  fwrite(&(w->N), sizeof(size_t), 1, fp);

  for (i = 0; i < w->nfreq; ++i)
    gsl_matrix_complex_fwrite(fp, w->modes[i]);

  fclose(fp);

  return GSL_SUCCESS;
}

invert_tmode_workspace *
invert_tmode_read_binary(const char * filename)
{
  FILE *fp;
  size_t i;
  invert_tmode_workspace * w;
  size_t nfreq, N;
  size_t * nmodes;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "invert_tmode_read_binary: unable to open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  fread(&nfreq, sizeof(size_t), 1, fp);

  nmodes = malloc(nfreq * sizeof(size_t));
  fread(nmodes, sizeof(size_t), nfreq, fp);

  fread(&N, sizeof(size_t), 1, fp);

  w = invert_tmode_alloc(nfreq, nmodes, N);

  for (i = 0; i < nfreq; ++i)
    gsl_matrix_complex_fread(fp, w->modes[i]);

  fclose(fp);

  return w;
}
