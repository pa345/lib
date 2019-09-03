/*
 * invert_tmode.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

#include <common/common.h>
#include <common/interp.h>
#include <satdata/satdata.h>

#include "invert_tmode.h"

/*
invert_tmode_alloc()
  Allocate invert_tmode_workspace

Inputs: nfreq  - number of frequency bands
        freqs  - frequences in cpd, length nfreq (can be NULL to be filled in later)
        nmodes - number of modes per frequency band, length nfreq
        N      - length of time series for each mode

Return: pointer to workspace
*/

invert_tmode_workspace *
invert_tmode_alloc(const size_t nfreq, const double freqs[], const size_t nmodes[], const size_t N)
{
  invert_tmode_workspace *w;
  size_t i;

  w = calloc(1, sizeof(invert_tmode_workspace));
  if (!w)
    return 0;

  w->nfreq = nfreq;
  w->N = N;

  w->freqs = malloc(nfreq * sizeof(double));
  w->nmodes = malloc(nfreq * sizeof(size_t));
  if (!w->freqs || !w->nmodes)
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
      if (freqs != NULL)
        w->freqs[i] = freqs[i];

      w->nmodes[i] = nmodes[i];

      w->modes[i] = gsl_matrix_complex_alloc(N, nmodes[i]);
      if (!w->modes[i])
        {
          invert_tmode_free(w);
          fprintf(stderr, "invert_tmode_alloc: failed to allocate matrix\n");
          return 0;
        }
    }

  w->t = malloc(N * sizeof(double));

  return w;
}

void
invert_tmode_free(invert_tmode_workspace *w)
{
  if (w->freqs)
    free(w->freqs);

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

  if (w->t)
    free(w->t);

  free(w);
}

/*
invert_tmode_get()
  Return value of temporal mode for a specified frequency band and mode number

Inputs: t    - timestamp (CDF_EPOCH)
        f    - frequency band in [0,nfreq-1]
        mode - mode number in frequency band f, in [0,nmodes[f]-1]
        w    - workspace

Return: alpha_{f,mode}(t)

Notes:
1) Temporal mode is interpolated to timestamp

2) Index acceleration lookup is not used, to keep function thread-safe
*/

gsl_complex
invert_tmode_get(const double t, const size_t f, const size_t mode, const invert_tmode_workspace * w)
{
  if (f >= w->nfreq)
    {
      GSL_ERROR_VAL ("invalid frequency band", GSL_EINVAL, GSL_COMPLEX_ZERO);
    }
  else if (mode >= w->nmodes[f])
    {
      GSL_ERROR_VAL ("invalid mode number", GSL_EINVAL, GSL_COMPLEX_ZERO);
    }
  else if (t < w->t[0] || t > w->t[w->N - 1])
    {
      fprintf(stderr, "t = %ld [%ld,%ld]\n", epoch2timet(t), epoch2timet(w->t[0]), epoch2timet(w->t[w->N-1]));
      GSL_ERROR_VAL ("t outside temporal mode range", GSL_EINVAL, GSL_COMPLEX_ZERO);
    }
  else
    {
      const size_t idx = gsl_interp_bsearch(w->t, t, 0, w->N - 1);
      gsl_complex z1 = gsl_matrix_complex_get(w->modes[f], idx, mode);
      gsl_complex z2 = gsl_matrix_complex_get(w->modes[f], idx + 1, mode);
      gsl_complex z;

      /* interpolate to time t */
      GSL_REAL(z) = interp1d(w->t[idx], w->t[idx + 1], GSL_REAL(z1), GSL_REAL(z2), t);
      GSL_IMAG(z) = interp1d(w->t[idx], w->t[idx + 1], GSL_IMAG(z1), GSL_IMAG(z2), t);

      return z;
    }
}

/*
invert_tmode_read_ascii()
  Read ASCII file of temporal mode data

Inputs: filename - filename to read
        ifreq    - frequency band, in [0,nfreq-1]
        mode     - mode number, in [0,nmodes[freq]-1]
        freq     - (output) frequency of mode in cpd
        w        - workspace
*/

int
invert_tmode_read_ascii(const char * filename, const size_t ifreq, const size_t mode, double * freq, invert_tmode_workspace * w)
{
  if (ifreq >= w->nfreq)
    {
      GSL_ERROR ("invalid frequency band", GSL_EBADLEN);
    }
  else if (mode >= w->nmodes[ifreq])
    {
      GSL_ERROR ("invalid mode number", GSL_EBADLEN);
    }
  else
    {
      FILE *fp;
      gsl_vector_complex_view v = gsl_matrix_complex_column(w->modes[ifreq], mode);
      char buffer[2048];
      int read_freq = 0;
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

          if (!read_freq)
            {
              if (*buffer == '#')
                {
                  c = sscanf(buffer, "# Frequency: %lf", freq);
                  if (c == 1)
                    read_freq = 1;
                }

              continue;
            }

          c = sscanf(buffer, "%ld %lf %lf", &t, &re, &im);
          if (c != 3)
            continue;

          w->t[n] = satdata_timet2epoch(t);
          gsl_vector_complex_set(&v.vector, n, gsl_complex_rect(re, im));

          if (n++ >= w->modes[ifreq]->size1)
            {
              fprintf(stderr, "invert_tmode_read_ascii: mode matrix size1 too small [%zu]\n", n);
              return -1;
            }
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
  fwrite(w->freqs, sizeof(double), w->nfreq, fp);
  fwrite(w->nmodes, sizeof(size_t), w->nfreq, fp);
  fwrite(&(w->N), sizeof(size_t), 1, fp);
  fwrite(w->t, sizeof(double), w->N, fp);

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
  double * freqs;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "invert_tmode_read_binary: unable to open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  fread(&nfreq, sizeof(size_t), 1, fp);

  freqs = malloc(nfreq * sizeof(double));
  fread(freqs, sizeof(double), nfreq, fp);

  nmodes = malloc(nfreq * sizeof(size_t));
  fread(nmodes, sizeof(size_t), nfreq, fp);

  fread(&N, sizeof(size_t), 1, fp);

  w = invert_tmode_alloc(nfreq, freqs, nmodes, N);

  fread(w->t, sizeof(double), N, fp);

  for (i = 0; i < nfreq; ++i)
    gsl_matrix_complex_fread(fp, w->modes[i]);

  fclose(fp);

  return w;
}

/*
invert_tmode_print()
  Print temporal modes to output directory
*/

int
invert_tmode_print(const char * dir_prefix, const invert_tmode_workspace * w)
{
  int s = 0;
  const double t0 = w->t[0] + 0.5*3.6e6;
  const double t1 = w->t[w->N - 1];
  const double dt = 1.0 * 3.6e6; /* time step in ms */
  char buf[2048];
  size_t f;
  
  for (f = 0; f < w->nfreq; ++f)
    {
      size_t mode;

      for (mode = 0; mode < w->nmodes[f]; ++mode)
        {
          double t;
          FILE *fp;
          size_t i;

          sprintf(buf, "%s/tmode_%02zu_%02zu.txt", dir_prefix, f + 1, mode + 1);
          fp = fopen(buf, "w");
          if (!fp)
            continue;

          i = 1;
          fprintf(fp, "# Frequency:   %f [cpd]\n", w->freqs[f]);
          fprintf(fp, "# Mode number: %zu\n", mode);
          fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
          fprintf(fp, "# Field %zu: real part of temporal mode\n", i++);
          fprintf(fp, "# Field %zu: imag part of temporal mode\n", i++);

          for (t = t0; t <= t1; t += dt)
            {
              gsl_complex z = invert_tmode_get(t, f, mode, w);
              fprintf(fp, "%ld %.12e %.12e\n", epoch2timet(t), GSL_REAL(z), GSL_IMAG(z));
            }

          fclose(fp);
        }
    }

  return s;
}
