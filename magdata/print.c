/*
 * print.c
 *
 * Print contents of magdata file
 *
 * ./print <-i magdata_file>
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include <common/common.h>
#include <msynth/msynth.h>

#include "magdata.h"

/*
magdata_print2()
  Output data points in ASCII format, only printing
time, position and magnetic measurements for points selected
for main field modeling

Inputs: filename - where to store data
        data     - data

Return: success/error
*/

int
magdata_print2(const char *filename, const magdata *data)
{
  int s = 0;
  size_t i;
  FILE *fp;
  msynth_workspace *msynth_p;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "magdata_print2: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  msynth_p = msynth_read(MSYNTH_BOUMME_FILE);
  msynth_set(1, 15, msynth_p);

  if (data->name[0] != '\0')
    fprintf(fp, "# Name: %s\n", data->name);

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (seconds since 01-01-1970 00:00:00 UTC\n", i++);
  fprintf(fp, "# Field %zu: geocentric radius (km)\n", i++);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: B_x NEC (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y NEC (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z NEC (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_x minus core field (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y minus core field (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z minus core field (nT)\n", i++);
  fprintf(fp, "# Field %zu: d/dt B_x (nT/year)\n", i++);
  fprintf(fp, "# Field %zu: d/dt B_y (nT/year)\n", i++);
  fprintf(fp, "# Field %zu: d/dt B_z (nT/year)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double B_core[4], B[4];

      if (data->flags[i] & MAGDATA_FLG_DISCARD)
        continue;

      if (!(data->flags[i] & MAGDATA_FLG_FIT_MF))
        continue;

#if 0
      if (!(data->flags[i] & MAGDATA_FLG_Z))
        continue;
#else
      if (!(data->flags[i] & MAGDATA_FLG_DZDT))
        continue;
#endif

      /* subtract external POMME */
      magdata_residual(i, B, data);

      msynth_eval(satdata_epoch2year(data->t[i]), data->r[i], data->theta[i], data->phi[i], B_core, msynth_p);

      if ((i > 0) && (data->flags[i] & MAGDATA_FLG_TRACK_START))
        fprintf(fp, "\n");

      fprintf(fp, "%ld %.4f %.3f %.3f %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
              unix_time,
              data->r[i],
              wrap180(data->phi[i] * 180.0 / M_PI),
              90.0 - data->theta[i] * 180.0 / M_PI,
              data->Bx_nec[i],
              data->By_nec[i],
              data->Bz_nec[i],
              B[0] - B_core[0],
              B[1] - B_core[1],
              B[2] - B_core[2],
              data->dXdt_nec[i],
              data->dYdt_nec[i],
              data->dZdt_nec[i]);
    }

  msynth_free(msynth_p);

  fclose(fp);

  return s;
}

int
main(int argc, char *argv[])
{
  magdata *data = NULL;
  char *outfile = "data.dat";
  struct timeval tv0, tv1;
  int c;

  while ((c = getopt(argc, argv, "i:o:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data = magdata_read(optarg, NULL);
            if (!data)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      fprintf(stderr, "Usage: %s <-i magdata_file> [-o output_file]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: writing output to %s...", outfile);
#if 0
  magdata_print(outfile, data);
#else
  magdata_print2(outfile, data);
#endif
  fprintf(stderr, "done\n");

  magdata_free(data);

  return 0;
} /* main() */
