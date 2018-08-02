#include "euler_calc.h"

/*
stage2_euler()
  Perform Euler angle calculation for unflagged satellite data

Inputs: res_file - filename for residuals
        data     - satellite data
        track_p  - track workspace

Return: success/error

Notes:
1) All unflagged points in 'data' are input to the Euler angle calculation
*/

int
stage2_euler(const char *res_file, const satdata_mag * data, const track_workspace * track_p)
{
  int s = 0;
  size_t nflagged = satdata_nflagged(data);
  size_t n = data->n - nflagged;
  size_t i, j;
  euler_calc_workspace *euler_p = euler_calc_alloc(n);
  struct timeval tv0, tv1;

  /* add unflagged data to euler workspace */
  fprintf(stderr, "stage2_euler: adding data for Euler angles...");
  gettimeofday(&tv0, NULL);

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;

      if (tptr->flags)
        continue;

      for (j = start_idx; j <= end_idx; ++j)
        {
          double B_VFM[3], B_model[4];

          if (data->flags[j])
            continue;

          if (fabs(data->qdlat[j]) > 55.0)
            continue;

          B_VFM[0] = SATDATA_VEC_X(data->B_VFM, j);
          B_VFM[1] = SATDATA_VEC_Y(data->B_VFM, j);
          B_VFM[2] = SATDATA_VEC_Z(data->B_VFM, j);

          satdata_mag_model(j, B_model, data);

          euler_calc_add(data->t[j], data->qdlat[j], B_VFM, B_model, &(data->q[4*j]), euler_p);
        }
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  euler_calc_proc(euler_p);

  if (res_file)
    {
      fprintf(stderr, "stage2_euler: writing Euler residuals to %s...", res_file);
      euler_calc_print_residuals(res_file, euler_p);
      fprintf(stderr, "done\n");
    }

  euler_calc_free(euler_p);

  return s;
}
