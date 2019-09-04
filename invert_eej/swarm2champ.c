#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mainlib/ml_satdata.h>

int
main()
{
  const char *infile = "champ2.cdf";
  const char *outfile = "swarm2.cdf";
  satdata_mag *data;

  fprintf(stderr, "reading %s...", infile);
  data = satdata_champ_mag_read(infile, NULL);
  fprintf(stderr, "done (%zu points read)\n", data->n);

  fprintf(stderr, "writing %s...", outfile);
  satdata_swarm_write(1, outfile, data);
  fprintf(stderr, "done\n");

  satdata_mag_free(data);

  return 0;
}
