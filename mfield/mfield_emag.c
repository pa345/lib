/*
 * mfield_emag.c
 *
 * Read EMAG2 grid and store in magdata format
 *
 * The result is an output file in magdata format containing all data point to
 * be used in the modeling. All data points will have a MAGDATA_FLG_FIT_xxx flag
 * set, and other flags will vary.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <zlib.h>

#include <gsl/gsl_math.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_geo.h>
#include <mainlib/ml_magdata.h>

/* from: zcat EMAG2_V2.xyz.gz | wc -l */
#define EMAG_N           52601353

typedef struct
{
  double lon[EMAG_N];
  double lat[EMAG_N];
  double F[EMAG_N];
  size_t n;
} emag_data;

emag_data *
read_emag(const char *filename)
{
  gzFile fp;
  emag_data *data;
  char buf[2048];
  size_t n = 0;

  fp = gzopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "read_emag: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  data = malloc(sizeof(emag_data));

  while (gzgets(fp, buf, 2048) != NULL)
    {
      int c;
      double lon, lat, F;

      c = sscanf(buf, "%lf %lf %lf",
                 &lon,
                 &lat,
                 &F);
      if (c < 3)
        continue;

      data->lon[n] = lon;
      data->lat[n] = lat;
      data->F[n] = F;

      if (++n > EMAG_N)
        {
          fprintf(stderr, "read_emag: error: EMAG_N too small\n");
          break;
        }
    }

  gzclose(fp);

  assert(n == EMAG_N);
  data->n = n;

  return data;
}

magdata *
copy_data(const emag_data *data)
{
  const double R = R_EARTH_KM;
  const size_t n = data->n;
  const double h = 4.0; /* geodetic height above ellipsoid */
  magdata *mdata = magdata_alloc(n, R);
  size_t i;
  magdata_datum datum;
  int s;

  mdata->global_flags = MAGDATA_GLOBFLG_SCALAR_GRID;

  magdata_datum_init(&datum);

  for (i = 0; i < n; ++i)
    {
      double r, geoc_lat;
      double latd = data->lat[i] * M_PI / 180.0;

      /* convert geodetic to geocentric */
      geodetic2geo(latd, h, &geoc_lat, &r);

      /* discard points exactly at the north pole, since we cannot
       * compute Legendre functions here for MF modeling */
      if (fabs(geoc_lat) == M_PI / 2.0)
        continue;

      datum.r = r;
      datum.theta = M_PI / 2.0 - geoc_lat;
      datum.phi = data->lon[i] * M_PI / 180.0;
      datum.F = data->F[i];
      datum.flags = MAGDATA_FLG_F | MAGDATA_FLG_FIT_MF;

      s = magdata_add(&datum, mdata);
      if (s)
        {
          fprintf(stderr, "copy_data: error adding data point to magdata\n");
          break;
        }
    }

  return mdata;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --emag_file   | -i emag_grid_file_gz  - gzipped EMAG grid file\n");
  fprintf(stderr, "\t --output_file | -o output_file        - output file\n");
}

int
main(int argc, char *argv[])
{
  const char *datamap_file = "datamap.dat";
  char *output_file = NULL;
  emag_data *data = NULL;
  magdata *mdata;
  struct timeval tv0, tv1;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "emag_file", required_argument, NULL, 'i' },
          { "output_file", required_argument, NULL, 'o' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:o:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = read_emag(optarg);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu points read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));
            break;

          case 'o':
            output_file = optarg;
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: converting to magdata format...");
  gettimeofday(&tv0, NULL);
  mdata = copy_data(data);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  magdata_init(mdata);

  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

#if 0
  fprintf(stderr, "main: writing data map to %s...", datamap_file);
  magdata_map(datamap_file, mdata);
  fprintf(stderr, "done\n");
#endif

  if (output_file)
    {
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");
    }

  magdata_free(mdata);

  return 0;
}
