#include "fluxcal.h"

/*
stage2_calibrate()
  Calculate scalar calibration parameters for all unflagged satellite data

Inputs: data    - satellite data
        track_p - track workspace
        c       - (output) scalar calibration parameters (indexed by FLUXCAL_IDX_xxx)
        rms     - (output) scalar residual rms after calibration (nT)

Return: success/error

Notes:
1) All unflagged points in 'data' are input to the scalar calibration

2) Satellite dataset is left unmodified (calibration is not applied to data)
*/

int
stage2_calibrate(const char *res_file, const satdata_mag * data, const track_workspace * track_p,
                 gsl_vector * c, double *rms)
{
  int s = 0;
  size_t nflagged = satdata_nflagged(data);
  size_t n = data->n - nflagged;
  size_t i, j;
  fluxcal_workspace *fluxcal_p = fluxcal_alloc(n);
  struct timeval tv0, tv1;

  /* add unflagged data to fluxcal workspace */
  fprintf(stderr, "stage2_calibrate: adding data for scalar calibration...");
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

          fluxcal_add_datum(data->t[j], B_VFM, B_model[3], fluxcal_p);
        }
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* set initial values of calibration parameters */
  gsl_vector_set(c, FLUXCAL_IDX_SX, 1.0);
  gsl_vector_set(c, FLUXCAL_IDX_SY, 1.0);
  gsl_vector_set(c, FLUXCAL_IDX_SZ, 1.0);
  gsl_vector_set(c, FLUXCAL_IDX_OX, 0.0);
  gsl_vector_set(c, FLUXCAL_IDX_OY, 0.0);
  gsl_vector_set(c, FLUXCAL_IDX_OZ, 0.0);
  gsl_vector_set(c, FLUXCAL_IDX_U1, 0.0);
  gsl_vector_set(c, FLUXCAL_IDX_U2, 0.0);
  gsl_vector_set(c, FLUXCAL_IDX_U3, 0.0);

  fprintf(stderr, "stage2_calibrate: performing scalar calibration...");
  gettimeofday(&tv0, NULL);
  fluxcal_proc(c, fluxcal_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  *rms = fluxcal_rms(fluxcal_p);

  if (res_file)
    {
      fprintf(stderr, "stage2_calibrate: writing scalar calibration residuals to %s...", res_file);
      fluxcal_print_residuals(res_file, c, fluxcal_p);
      fprintf(stderr, "done\n");
    }

  fluxcal_free(fluxcal_p);

  return s;
}
