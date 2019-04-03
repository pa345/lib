#include <att/att_calc.h>

/*
stage2_align()
  Calculate alignment parameters between VFM and spacecraft-fixed frame

Inputs: data    - satellite data
        track_p - track workspace
        c       - (output) alignment parameters (indexed by FLUXCAL_IDX_xxx)

Return: success/error

Notes:
1) All unflagged points in 'data' are input to the alignment calculation
2) Satellite dataset is left unmodified (alignment is not applied to data)
*/

int
stage2_align(const satdata_mag * data, const track_workspace * track_p,
             gsl_vector * c, att_calc_workspace * w)
{
  int s = 0;
  size_t i, j;
  att_workspace * att_quat_p = att_alloc(att_quaternion);
  struct timeval tv0, tv1;

  /* add unflagged data to fluxcal workspace */
  fprintf(stderr, "stage2_align: adding data for alignment...");
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
          gsl_vector_view q = gsl_vector_view_array(&(data->q[4 * j]), 4);
          double *B_VFM = &(data->B_VFM[3 * j]);
          double B_model[4];
          gsl_vector_view v = gsl_vector_view_array(B_model, 3);

          if (data->flags[j])
            continue;

          if (fabs(data->qdlat[j]) > 50.0)
            continue;

          satdata_mag_model(j, B_model, data);

          /* rotate B_model to spacecraft frame */
          att_rotate_inverse(&q.vector, &v.vector, &v.vector, att_quat_p);

          att_calc_add(B_VFM, B_model, w);
        }
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, %zu data added)\n", time_diff(tv0, tv1), w->n);

  fprintf(stderr, "stage2_align: performing alignment...");
  gettimeofday(&tv0, NULL);
  att_calc_proc(w);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  gsl_vector_memcpy(c, w->c);

  att_free(att_quat_p);

  return s;
}

/*
stage2_align_vfm2nec()
  Compute B_nec vector by applying alignment and attitude rotation
to data

Inputs: c    - alignment parameters, computed by stage2_align()
        data - satellite data
        w    - workspace

Return: success/error

Notes:
1) Only data which is not time-flagged is computed, in case the
alignment parameters are computed for different time periods
*/

static int
stage2_align_vfm2nec(const gsl_vector * c, satdata_mag *data, att_calc_workspace * w)
{
  att_workspace * att_quat_p = att_alloc(att_quaternion);
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      gsl_vector_const_view B_VFM = gsl_vector_const_view_array(&(data->B_VFM[3 * i]), 3);
      gsl_vector_view B_NEC = gsl_vector_view_array(&(data->B[3 * i]), 3);
      gsl_vector_view q = gsl_vector_view_array(&(data->q[4 * i]), 4);

      /* ignore time flagged data */
      if (data->flags[i] & SATDATA_FLG_TIME)
        continue;

      att_rotate(c, &B_VFM.vector, &B_NEC.vector, w->att_workspace_p);
      att_rotate(&q.vector, &B_NEC.vector, &B_NEC.vector, att_quat_p);
    }

  att_free(att_quat_p);

  return 0;
}
