/*
 * compare_plp.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

#include <gsl/gsl_math.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_bsearch.h>
#include <mainlib/ml_interp.h>

#include "ccorr.h"
#include "poltor.h"

int
write_corr(const char *outfile, const magdata *mdata)
{
  int s = 0;
  const double qdmin = -90.0;
  const double qdmax = 90.0;
  const double qdstep = 5.0;
  const size_t nqd = (size_t) ((qdmax - qdmin) / qdstep);
  ccorr_workspace *ccorr_p;
  size_t i, j;
  FILE *fp;

  fp = fopen(outfile, "w");
  if (!fp)
    {
      fprintf(stderr, "write_corr: unable to open %s: %s\n",
              outfile, strerror(errno));
      return -1;
    }

  ccorr_p = ccorr_alloc(qdmin, qdmax, nqd);

  for (i = 0; i < mdata->n; ++i)
    {
      if (mdata->flags[i] & MAGDATA_FLG_TRACK_START)
        {
          for (j = 0; j < ccorr_p->nbins; ++j)
            {
              double qdlat, r;
              size_t n;

              ccorr_xval(j, &qdlat, ccorr_p);
              r = ccorr_r(qdlat, ccorr_p);
              n = ccorr_n(qdlat, ccorr_p);

              fprintf(fp, "%f %f %zu\n", qdlat, r, n);
            }

          fprintf(fp, "\n\n");
          fflush(fp);

          ccorr_reset(ccorr_p);
        }

      if (mdata->ne[i] != 0.0)
        {
          double B[4];

          magdata_residual(i, B, mdata);
          ccorr_add(mdata->qdlat[i], B[0], mdata->ne[i], ccorr_p);
        }
    }

  ccorr_free(ccorr_p);
  fclose(fp);

  return s;
} /* write_corr() */

int
fill_dens(magdata *mdata, const satdata_lp *dens_data)
{
  int s = 0;
  size_t i;

  /* initialize ne array */
  for (i = 0; i < mdata->n; ++i)
    mdata->ne[i] = 0.0;

  for (i = 0; i < mdata->n; ++i)
    {
      double t = mdata->t[i];
      size_t idx;
      double dt, ne;
      double t0, t1;

      if (t < dens_data->t[0] || t > dens_data->t[dens_data->n - 1])
        continue;

      idx = bsearch_double(dens_data->t, t, 0, dens_data->n);

      t0 = dens_data->t[idx];
      t1 = dens_data->t[idx + 1];

      dt = t1 - t0;
      dt /= 1000.0; /* convert to s */

      /* check for data gaps in PLP measurements */
      if (fabs(dt) > 30.0)
        continue;

      ne = interp1d(t0, t1, dens_data->ne[idx], dens_data->ne[idx + 1], t);

      mdata->ne[i] = ne;
    }

  return s;
}

int
main(int argc, char *argv[])
{
  magdata *mdata = NULL;
  satdata_lp *lp_data = NULL;
  int c;

  while ((c = getopt(argc, argv, "i:c:s:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading magnetic data from %s...", optarg);
            mdata = magdata_read(optarg, NULL);
            fprintf(stderr, "done\n");
            break;

          case 'c':
            fprintf(stderr, "main: reading CHAMP PLP data from %s...", optarg);
            lp_data = satdata_champ_plp_read_idx(optarg, 0);
            fprintf(stderr, "done\n");
            break;

          case 's':
            fprintf(stderr, "main: reading Swarm LP data from %s...", optarg);
            lp_data = satdata_swarm_lp_read_idx(optarg, 0);
            fprintf(stderr, "done\n");
            break;

          default:
            break;
        }
    }

  if (!mdata || !lp_data)
    {
      fprintf(stderr, "Usage: %s <-i mag_input_file> [-c champ_lp_idx_file] [-s swarm_lp_idx_file]\n", argv[0]);
      exit(1);
    }

  fill_dens(mdata, lp_data);

  {
    char *outfile = "corr.dat";

    fprintf(stderr, "main: writing correlation data to %s...", outfile);
    write_corr(outfile, mdata);
    fprintf(stderr, "done\n");
  }

  {
    char *outfile = "plp.dat";

    fprintf(stderr, "main: writing data to %s...", outfile);
    magdata_print(outfile, mdata);
    fprintf(stderr, "done\n");
  }

  return 0;
} /* main() */
