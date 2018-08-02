#include <gsl/gsl_filter.h>
#include <gsl/gsl_movstat.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

/*
stage2_correct_spikes()
  Detect and correct spikes in VFM time series data

Methodology:

1. Rotate main field model into VFM frame with given quaternions
2. Compute VFM residuals: r = B_VFM - B_model
3. Impulse rejection filter each component of VFM residual r
4. Detect outliers in new residuals: r - filter

Inputs: header      - 1 to write file header
        fp          - output file pointer
        window_size - window size for impulse rejection filter
        nsigma      - sigma multiplier for identifying spikes
        start_idx   - start index for processing
        end_idx     - end index for processing
        nspikes     - (output) number of spikes found in each component
        data        - satellite data

Notes:
1) On output, data->B_VFM is updated to remove detected spikes
2) On output, data->F is updated with new |B_VFM|
3) The number of spikes in [start_idx,end_idx] is added to nspikes, so
   the nspikes array should be initialized prior to this function
*/

static int
stage2_correct_spikes(const int header, FILE *fp, const size_t window_size, const double nsigma[3],
                      const size_t start_idx, const size_t end_idx,
                      size_t nspikes[3], double min_spike[3], double max_spike[3], satdata_mag * data)
{
  const size_t n = end_idx - start_idx + 1;
  const double min_threshold = 3.0; /* minimum required amplitude for spike (nT) */
  gsl_vector *input[3], *output[3];
  gsl_vector_int *ioutlier[3];
  gsl_filter_impulse_workspace *impulse_p;
  gsl_vector *xmedian[3], *xsigma[3];
  size_t noutlier[3];
  size_t i, j;

  if (header && fp)
    {
      i = 1;
      fprintf(fp, "# Field %zu: timestamp\n", i++);
      fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: original VFM X residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original VFM Y residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original VFM Z residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original F residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: median VFM X (nT)\n", i++);
      fprintf(fp, "# Field %zu: median VFM Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: median VFM Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: sigma VFM X (nT)\n", i++);
      fprintf(fp, "# Field %zu: sigma VFM Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: sigma VFM Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: upper VFM X limit (nT)\n", i++);
      fprintf(fp, "# Field %zu: lower VFM X limit (nT)\n", i++);
      fprintf(fp, "# Field %zu: upper VFM Y limit (nT)\n", i++);
      fprintf(fp, "# Field %zu: lower VFM Y limit (nT)\n", i++);
      fprintf(fp, "# Field %zu: upper VFM Z limit (nT)\n", i++);
      fprintf(fp, "# Field %zu: lower VFM Z limit (nT)\n", i++);
      fprintf(fp, "# Field %zu: outlier found in VFM X (nT)\n", i++);
      fprintf(fp, "# Field %zu: outlier found in VFM Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: outlier found in VFM Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected VFM X residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected VFM Y residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected VFM Z residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected F residual (nT)\n", i++);
      return 0;
    }

  impulse_p = gsl_filter_impulse_alloc(window_size);

  for (i = 0; i < 3; ++i)
    {
      input[i] = gsl_vector_alloc(n);
      output[i] = gsl_vector_alloc(n);
      xmedian[i] = gsl_vector_alloc(n);
      xsigma[i] = gsl_vector_alloc(n);
      ioutlier[i] = gsl_vector_int_calloc(n);
    }

  /* compute vector residuals in VFM frame with a-priori model */
  for (i = 0; i < n; ++i)
    {
      size_t didx = i + start_idx;
      double *q = &(data->q[4 * didx]);
      double *B_VFM = &(data->B_VFM[3 * didx]);
      double B_model[4], B_model_VFM[3];

      satdata_mag_model(didx, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      /* compute residuals in VFM frame */
      gsl_vector_set(input[0], i, B_VFM[0] - B_model_VFM[0]);
      gsl_vector_set(input[1], i, B_VFM[1] - B_model_VFM[1]);
      gsl_vector_set(input[2], i, B_VFM[2] - B_model_VFM[2]);
    }

  /* apply impulse detection filter to VFM residuals */

  for (j = 0; j < 3; ++j)
    {
      gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_QN, nsigma[j], input[j], output[j],
                         xmedian[j], xsigma[j], &noutlier[j], ioutlier[j], impulse_p);
    }

  for (i = 0; i < n; ++i)
    {
      size_t didx = i + start_idx;
      double *q = &(data->q[4 * didx]);
      double *B_VFM = &(data->B_VFM[3 * didx]);
      double B_model[4], B_model_VFM[3], B_VFM_orig[4];

      satdata_mag_model(didx, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      for (j = 0; j < 3; ++j)
        {
          int outlier = gsl_vector_int_get(ioutlier[j], i);

          /* save original B_VFM for later printing */
          B_VFM_orig[j] = B_VFM[j];

          if (outlier)
            {
              /* correct B_VFM[j] by spike amplitude */
              double diff = gsl_vector_get(input[j], i) - gsl_vector_get(output[j], i);
              double absdiff = fabs(diff);

              /* only correct spike if it is larger than a threshold value */
              /* only correct spikes for X/Y at low/mid-latitudes; high-latitudes is not reliable */
              if ((j == 2 || fabs(data->qdlat[didx]) < 55.0) && (absdiff > min_threshold))
                {
                  B_VFM[j] -= diff;

                  if (absdiff < min_spike[j])
                    min_spike[j] = absdiff;
                  if (absdiff > max_spike[j])
                    max_spike[j] = absdiff;
                }
              else
                {
                  --noutlier[j];
                  gsl_vector_int_set(ioutlier[j], i, 0);
                }
            }
        }

      /* save original F and recompute new F */
      B_VFM_orig[3] = data->F[didx];
      data->F[didx] = gsl_hypot3(B_VFM[0], B_VFM[1], B_VFM[2]);

      if (fp)
        {
          /* separate tracks with newlines */
          if (didx > 0 && data->flags[didx] & SATDATA_FLG_TRACK_START)
            fprintf(fp, "\n\n");

          fprintf(fp, "%ld %8.4f %8.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %d %d %d %10.4f %10.4f %10.4f %10.4f\n",
                  satdata_epoch2timet(data->t[didx]),
                  data->latitude[didx],
                  data->qdlat[didx],
                  gsl_vector_get(input[0], i),
                  gsl_vector_get(input[1], i),
                  gsl_vector_get(input[2], i),
                  B_VFM_orig[3] - B_model[3],
                  gsl_vector_get(xmedian[0], i),
                  gsl_vector_get(xmedian[1], i),
                  gsl_vector_get(xmedian[2], i),
                  gsl_vector_get(xsigma[0], i),
                  gsl_vector_get(xsigma[1], i),
                  gsl_vector_get(xsigma[2], i),
                  gsl_vector_get(xmedian[0], i) + nsigma[0] * gsl_vector_get(xsigma[0], i),
                  gsl_vector_get(xmedian[0], i) - nsigma[0] * gsl_vector_get(xsigma[0], i),
                  gsl_vector_get(xmedian[1], i) + nsigma[1] * gsl_vector_get(xsigma[1], i),
                  gsl_vector_get(xmedian[1], i) - nsigma[1] * gsl_vector_get(xsigma[1], i),
                  gsl_vector_get(xmedian[2], i) + nsigma[2] * gsl_vector_get(xsigma[2], i),
                  gsl_vector_get(xmedian[2], i) - nsigma[2] * gsl_vector_get(xsigma[2], i),
                  gsl_vector_int_get(ioutlier[0], i),
                  gsl_vector_int_get(ioutlier[1], i),
                  gsl_vector_int_get(ioutlier[2], i),
                  B_VFM[0] - B_model_VFM[0],
                  B_VFM[1] - B_model_VFM[1],
                  B_VFM[2] - B_model_VFM[2],
                  data->F[didx] - B_model[3]);
        }
    }

  for (j = 0; j < 3; ++j)
    nspikes[j] += noutlier[j];

  gsl_filter_impulse_free(impulse_p);

  for (i = 0; i < 3; ++i)
    {
      gsl_vector_free(input[i]);
      gsl_vector_free(output[i]);
      gsl_vector_free(xmedian[i]);
      gsl_vector_free(xsigma[i]);
    }

  return 0;
}
