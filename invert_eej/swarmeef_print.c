/*
 * swarmeef_print()
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <mainlib/ml_satdata.h>

#include "swarmeef.h"

int
print_data(const swarm_eef *data)
{
  size_t i;

  i = 1;
  printf("# Field %zu: decimal year (UT)\n", i++);
  printf("# Field %zu: latitude (degrees)\n", i++);
  printf("# Field %zu: longitude (degrees)\n", i++);
  printf("# Field %zu: EEF (mV/m)\n", i++);
  printf("# Field %zu: RelErr\n", i++);
  printf("# Field %zu: Flags\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      printf("%f %f %f %f %f %d\n",
             satdata_epoch2year(data->t[i]),
             data->latitude[i],
             data->longitude[i],
             data->EEF[i] * 1.0e3,
             data->RelErr[i],
             data->Flags[i]);
    }

  return 0;
}

int
main(int argc, char *argv[])
{
  int c;
  char *infile = NULL;
  swarm_eef *data;

  while ((c = getopt(argc, argv, "i:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s <-i eef_idx_file>\n", argv[0]);
            exit(1);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i eef_idx_file>\n", argv[0]);
      exit(1);
    }

  data = swarm_eef_read_idx(infile, 1);
  if (!data)
    {
      fprintf(stderr, "main: swarm_eef_read_idx failed\n");
      exit(1);
    }

  print_data(data);

  return 0;
}
