#include <dsp/gsl_dsp_filter.h>

/*
correct_track()
  Correct a track for possible (multiple) data jumps in the VFM frame

Inputs: tol    - tolerance to detect a jump (nT)
        B      - input with jumps
                 B[0] = filtered VFM X minus model (nT)
                 B[1] = filtered VFM Y minus model (nT)
                 B[2] = filtered VFM Z minus model (nT)
        offset - (input/output) matrix of size n-by-3
                 on output, offsets to add to B to remove jumps:
                   y = x + offset (jump-free)
        data   - satellite data
*/

static int
correct_track(const double tol[3], const gsl_vector * B[3], gsl_matrix * offset, satdata_mag * data)
{
  int s = 0;
  const double n = offset->size1;
  double alpha[3] = { 0.0, 0.0, 0.0 };       /* current offsets to add to vector components */
  double alpha_reset[3] = { 0.0, 0.0, 0.0 }; /* offsets computed after first track */
  int found_reset = 0;                       /* set to 1 after first track is processed */
  size_t j, k, l;
  double *off[3];
  const size_t sigma_window = 5;             /* number of samples for standard deviation */
  double *x = malloc(sigma_window * sizeof(double));

  for (j = 0; j < n - 1; ++j)
    {
      int dir = satdata_satdir(j, data->n, data->latitude);
      int dir_next = satdata_satdir(j + 1, data->n, data->latitude);
      double Bj[3], Bjp1[3];
      double off_next[3];
      double sample[3], sample_next[3];
      double delta[3], sigma[3];

      for (k = 0; k < 3; ++k)
        {
          Bj[k] = gsl_vector_get(B[k], j);
          Bjp1[k] = gsl_vector_get(B[k], j + 1);

          off[k] = gsl_matrix_ptr(offset, j, k);
          off_next[k] = gsl_matrix_get(offset, j + 1, k);

          *off[k] += alpha[k];
        }

#if 0
      /* compute std deviation of previous sigma_window samples */
      if (j <= 1)
        {
          sigma[0] = sigma[1] = sigma[2] = 1.0;
        }
      else
        {
          size_t nsamp = (j < sigma_window) ? j : sigma_window;

          for (k = 0; k < 3; ++k)
            {
              /* add the offsets to B(j-nsamp:j-1) before computing sigma */
              for (l = 0; l < nsamp; ++l)
                {
                  x[l] = gsl_vector_get(B[k], l + j - nsamp) + gsl_matrix_get(offset, l + j - nsamp, k);
                  /*x[l] = gsl_vector_get(B[k], l + j - nsamp);*/
                }

              sigma[k] = gsl_stats_sd(x, 1, nsamp);
            }
        }
#endif

#if 1
      if (j > 1)
        {
          double dt = (data->t[j] - data->t[j - 1]) / 1.0e3; /* dt in sec */

          if (fabs(dt) > 4.0)
            {
              continue;
            }
        }
#endif

      /* compute difference between adjacent samples and check for jump */
      for (k = 0; k < 3; ++k)
        {
          sample[k] = Bj[k] + *off[k] - alpha[k];
          sample_next[k] = Bjp1[k] + off_next[k];

          delta[k] = sample_next[k] - sample[k];

#if 0
          if (k == 2 && (fabs(delta[k]) > 15.0*sigma[k]))
            {
              printf("%ld %f %f %f %f %f\n",
                     satdata_epoch2timet(data->t[j]),
                     data->qdlat[j],
                     sigma[0],
                     sigma[1],
                     sigma[2],
                     fabs(delta[k]));
            }
#endif

          if (fabs(delta[k]) > tol[k])
            {
              int found_jump = 1;

              /* don't look for X/Y jumps in polar regions */
              if (k < 2 && fabs(data->qdlat[j]) > 55.0)
                found_jump = 0;

              if (found_jump)
                {
                  /* jump detected */
                  alpha[k] -= delta[k];
                }
            }
        }

#if 1
      /* reset alpha counter after each track */
      if (dir == 1 && dir != dir_next)
        {
          if (!found_reset)
            {
              /* store alpha values after first track - all future tracks will
               * be reset to these values */
              for (k = 0; k < 3; ++k)
                alpha_reset[k] = alpha[k];

              found_reset = 1;
            }
          else
            {
              /* reset alpha values - if I don't do this I find a significant
               * drift in the offset baseline over time */
              for (k = 0; k < 3; ++k)
                alpha[k] = alpha_reset[k];
            }
        }
#endif
    }

  /* account for the last sample */
  for (k = 0; k < 3; ++k)
    {
      off[k] = gsl_matrix_ptr(offset, n - 1, k);
      *off[k] += alpha[k];
    }

  free(x);

  return s;
}

/*
correct_jumps()
  Detect and correct jumps in VFM time series data

Methodology:

1. Rotate main field model into VFM frame with given quaternions
2. Compute VFM residuals: r = B_VFM - B_model
3. Median filter each component of VFM residual r
4. Detect outliers in new residuals: r - median_filter

Inputs: filename    - data file
        window_size - window size for median filter
        njumps      - (output) number of jumps found in each component
        data        - satellite data

Notes:
1) On output, data->B_VFM is updated to remove detected jumps
2) data->F is not updated
*/

static int
correct_jumps(const char *filename, const size_t window_size, size_t njumps[3], satdata_mag * data)
{
  FILE *fp = NULL;
  const size_t n = data->n;
  double tol[3];
  gsl_vector *input[3], *output[3];
  gsl_dsp_filter_median_workspace *median_p = gsl_dsp_filter_median_alloc(window_size);
  gsl_matrix *offset = gsl_matrix_calloc(n, 3);
  size_t i, j;

  tol[0] = 1000.0;
  tol[1] = 1000.0;
  tol[2] = 4.0;

  for (i = 0; i < 3; ++i)
    {
      njumps[i] = 0;
      input[i] = gsl_vector_alloc(n);
      output[i] = gsl_vector_alloc(n);
    }

  if (filename)
    {
      fp = fopen(filename, "w");

      i = 1;
      fprintf(fp, "# Field %zu: timestamp\n", i++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: original VFM X (nT)\n", i++);
      fprintf(fp, "# Field %zu: original VFM Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: original VFM Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: original F(nT)\n", i++);
      fprintf(fp, "# Field %zu: modeled VFM X (nT)\n", i++);
      fprintf(fp, "# Field %zu: modeled VFM Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: modeled VFM Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: modeled F(nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected VFM X (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected VFM Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected VFM Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected F (nT)\n", i++);
      fprintf(fp, "# Field %zu: satellite direction (+1 north, -1 south)\n", i++);
    }

  for (i = 0; i < n; ++i)
    {
      double *q = &(data->q[4 * i]);
      double *B_VFM = &(data->B_VFM[3 * i]);
      double B_model[4], B_model_VFM[3];

      satdata_mag_model(i, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      /* compute residuals in VFM frame */
      gsl_vector_set(input[0], i, B_VFM[0] - B_model_VFM[0]);
      gsl_vector_set(input[1], i, B_VFM[1] - B_model_VFM[1]);
      gsl_vector_set(input[2], i, B_VFM[2] - B_model_VFM[2]);
    }

  /* apply median filter to VFM residuals */
  for (j = 0; j < 3; ++j)
    gsl_dsp_filter_median(input[j], output[j], median_p);

  correct_track(tol, (const gsl_vector **) output, offset, data);

  for (i = 0; i < n; ++i)
    {
      double *q = &(data->q[4 * i]);
      double *B_VFM = &(data->B_VFM[3 * i]);
      double B_model[4], B_model_VFM[3];
      int dir = satdata_satdir(i, data->n, data->latitude);
      double F_orig = data->F[i];

      satdata_mag_model(i, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      /* compute residuals in VFM frame */
      gsl_vector_set(input[0], i, B_VFM[0] - B_model_VFM[0]);
      gsl_vector_set(input[1], i, B_VFM[1] - B_model_VFM[1]);
      gsl_vector_set(input[2], i, B_VFM[2] - B_model_VFM[2]);

      /* correct B_VFM */
      B_VFM[0] += gsl_matrix_get(offset, i, 0);
      B_VFM[1] += gsl_matrix_get(offset, i, 1);
      B_VFM[2] += gsl_matrix_get(offset, i, 2);

      /* recompute F */
      data->F[i] = gsl_hypot3(B_VFM[0], B_VFM[1], B_VFM[2]);

      if (fp)
        {
          fprintf(fp, "%ld %f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %d\n",
                  satdata_epoch2timet(data->t[i]),
                  data->qdlat[i],
                  B_VFM[0] - gsl_matrix_get(offset, i, 0),
                  B_VFM[1] - gsl_matrix_get(offset, i, 1),
                  B_VFM[2] - gsl_matrix_get(offset, i, 2),
                  F_orig,
                  B_model_VFM[0],
                  B_model_VFM[1],
                  B_model_VFM[2],
                  B_model[3],
                  B_VFM[0],
                  B_VFM[1],
                  B_VFM[2],
                  data->F[i],
                  dir);

          if (i < n - 1)
            {
              int dir_next = satdata_satdir(i + 1, data->n, data->latitude);
              if (dir != dir_next)
                {
                  fprintf(fp, "\n\n");
                }
            }
        }
    }

  gsl_dsp_filter_median_free(median_p);

  for (i = 0; i < 3; ++i)
    {
      gsl_vector_free(input[i]);
      gsl_vector_free(output[i]);
    }

  if (fp)
    fclose(fp);

  return 0;
}

/*
correct_jumps2()
  Detect and correct jumps in VFM time series data

Methodology:

1. Rotate main field model into VFM frame with given quaternions
2. Compute VFM residuals: r = B_VFM - B_model
3. Median filter each component of VFM residual r
4. Detect outliers in new residuals: r - median_filter

Inputs: filename  - data file
        H         - number of previous samples to include in moving MAD window
        njumps    - (output) number of jumps found in each component
        data      - satellite data

Notes:
1) On output, data->B_VFM is updated to remove detected jumps
2) data->F is not updated
*/

static int
correct_jumps2(const char *filename, const size_t H, size_t njumps[3], satdata_mag * data)
{
  int s = 0;
  FILE *fp = NULL;
  const size_t n = data->n;
  gsl_dsp_movmad_workspace * movmad_p = gsl_dsp_movmad_alloc2(H, 0);
  gsl_vector *input[3], *median[3], *mad[3];
  size_t i;

  if (filename)
    {
      fp = fopen(filename, "w");

      i = 1;
      fprintf(fp, "# Field %zu: timestamp\n", i++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: original VFM X residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original VFM Y residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original VFM Z residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: X window median (nT)\n", i++);
      fprintf(fp, "# Field %zu: Y window median (nT)\n", i++);
      fprintf(fp, "# Field %zu: Z window median (nT)\n", i++);
      fprintf(fp, "# Field %zu: X window MAD (nT)\n", i++);
      fprintf(fp, "# Field %zu: Y window MAD (nT)\n", i++);
      fprintf(fp, "# Field %zu: Z window MAD (nT)\n", i++);
      fprintf(fp, "# Field %zu: X first difference (nT)\n", i++);
      fprintf(fp, "# Field %zu: Y first difference (nT)\n", i++);
      fprintf(fp, "# Field %zu: Z first difference (nT)\n", i++);
      fprintf(fp, "# Field %zu: found jump in X (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: found jump in Y (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: found jump in Z (0 or 1)\n", i++);
      fprintf(fp, "# Field %zu: satellite direction (+1 north, -1 south)\n", i++);
    }

  for (i = 0; i < 3; ++i)
    {
      njumps[i] = 0;
      input[i] = gsl_vector_alloc(n);
      median[i] = gsl_vector_alloc(n);
      mad[i] = gsl_vector_alloc(n);
    }

  /* compute residuals in VFM frame (data - model) for each vector component */
  for (i = 0; i < n; ++i)
    {
      double *q = &(data->q[4 * i]);
      double *B_VFM = &(data->B_VFM[3 * i]);
      double B_model[4], B_model_VFM[3];

      satdata_mag_model(i, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      /* compute residuals in VFM frame */
      gsl_vector_set(input[0], i, B_VFM[0] - B_model_VFM[0]);
      gsl_vector_set(input[1], i, B_VFM[1] - B_model_VFM[1]);
      gsl_vector_set(input[2], i, B_VFM[2] - B_model_VFM[2]);
    }

  /* apply moving MAD to VFM residuals */
  for (i = 0; i < 3; ++i)
    gsl_dsp_movmad(input[i], median[i], mad[i], movmad_p);

  for (i = 0; i < n; ++i)
    {
      double *q = &(data->q[4 * i]);
      double *B_VFM = &(data->B_VFM[3 * i]);
      double B_model[4], B_model_VFM[3];
      int dir = satdata_satdir(i, data->n, data->latitude);
      double F_orig = data->F[i];
      double diff[3];
      int found_jump[3] = { 0, 0, 0 };
      size_t j;

      satdata_mag_model(i, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      /* compute first differences of VFM residuals */
      if (i < n - 1)
        {
          for (j = 0; j < 3; ++j)
            diff[j] = gsl_vector_get(input[j], i + 1) - gsl_vector_get(input[j], i);
        }
      else
        {
          for (j = 0; j < 3; ++j)
            diff[j] = gsl_vector_get(input[j], i) - gsl_vector_get(input[j], i - 1);
        }

      if (fabs(diff[2]) > 3.0 * gsl_vector_get(mad[2], i))
        found_jump[2] = 1;

      if (fp)
        {
          fprintf(fp, "%ld %f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %d %d %d %d\n",
                  satdata_epoch2timet(data->t[i]),
                  data->qdlat[i],
                  gsl_vector_get(input[0], i),
                  gsl_vector_get(input[1], i),
                  gsl_vector_get(input[2], i),
                  gsl_vector_get(median[0], i),
                  gsl_vector_get(median[1], i),
                  gsl_vector_get(median[2], i),
                  gsl_vector_get(mad[0], i),
                  gsl_vector_get(mad[1], i),
                  gsl_vector_get(mad[2], i),
                  diff[0],
                  diff[1],
                  diff[2],
                  found_jump[0],
                  found_jump[1],
                  found_jump[2],
                  dir);

          if (i < n - 1)
            {
              int dir_next = satdata_satdir(i + 1, data->n, data->latitude);
              if (dir != dir_next)
                {
                  fprintf(fp, "\n\n");
                }
            }
        }
    }

  if (fp)
    fclose(fp);

  for (i = 0; i < 3; ++i)
    {
      gsl_vector_free(input[i]);
      gsl_vector_free(median[i]);
      gsl_vector_free(mad[i]);
    }

  gsl_dsp_movmad_free(movmad_p);

  return s;
}

static int
stage2_correct_jumps(const char *jump_file, satdata_mag *data)
{
  int s = 0;
  size_t njumps[3];
  size_t i;

  fprintf(stderr, "\n");

  fprintf(stderr, "\t stage2_correct_jumps: fixing data jumps...");
#if 0
  correct_jumps(jump_file, 3, njumps, data);
#else
  correct_jumps2(jump_file, 5, njumps, data);
#endif
  fprintf(stderr, "done (data written to %s, %zu X jumps, %zu Y jumps, %zu Z jumps)\n",
          jump_file,
          njumps[0],
          njumps[1],
          njumps[2]);

  fprintf(stderr, "\t stage2_correct_jumps: recomputing F...");

  for (i = 0; i < data->n; ++i)
    {
      double *B_VFM = &(data->B_VFM[3 * i]);
      data->F[i] = gsl_hypot3(B_VFM[0], B_VFM[1], B_VFM[2]);
    }

  fprintf(stderr, "done\n");

  return s;
}

/* compute y = [x(2)-x(1) x(3)-x(2) ... x(n-1)-x(n-2) x(n)-x(n-1) x(n)-x(n-1)]; in-place y = x allowed */
static int
compute_diff(const gsl_vector * x, gsl_vector * y)
{
  const size_t n = x->size;
  double xi = gsl_vector_get(x, 0);
  size_t i;

  for (i = 0; i < n - 1; ++i)
    {
      double xip1 = gsl_vector_get(x, i + 1);

      gsl_vector_set(y, i, xip1 - xi);
      xi = xip1;
    }

  gsl_vector_set(y, n - 1, gsl_vector_get(y, n - 2));

  return 0;
}

static int
stage2_correct_jumps2(const int header, FILE *fp, const size_t start_idx, const size_t end_idx,
                      satdata_mag * data)
{
  const size_t n = end_idx - start_idx + 1;
  const double min_jump_threshold = 2.0; /* minimum required jump amplitude betweeen adjacent samples (nT) */
  const size_t K = 11;
  const double nsigma = 5.0;
  size_t njump[3];
  gsl_vector *input[3], *output_rmf[3], *output_fd[3], *output_hampel[3];
  gsl_vector *xmedian[3], *xsigma[3];
  gsl_vector_int *ioutlier[3];
  gsl_filter_rmf_workspace *rmf_p = gsl_filter_rmf_alloc(K);
  gsl_filter_hampel_workspace *hampel_p = gsl_filter_hampel_alloc(K);
  double cumjump[3] = { 0.0, 0.0, 0.0 }; /* cumulative jump sums */
  size_t i;

  if (header && fp)
    {
      i = 1;
      fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: original X VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original Y VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: original Z VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: RM-filtered X VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: RM-filtered Y VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: RM-filtered Z VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: first differenced RM-filtered X (nT)\n", i++);
      fprintf(fp, "# Field %zu: first differenced RM-filtered Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: first differenced RM-filtered Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: Hampel filtered X (nT)\n", i++);
      fprintf(fp, "# Field %zu: Hampel filtered Y (nT)\n", i++);
      fprintf(fp, "# Field %zu: Hampel filtered Z (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected X VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected Y VFM residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected Z VFM residual (nT)\n", i++);
      return 0;
    }

  for (i = 0; i < 3; ++i)
    {
      input[i] = gsl_vector_alloc(n);
      output_rmf[i] = gsl_vector_alloc(n);
      output_fd[i] = gsl_vector_alloc(n);
      output_hampel[i] = gsl_vector_alloc(n);
      xmedian[i] = gsl_vector_alloc(n);
      xsigma[i] = gsl_vector_alloc(n);
      ioutlier[i] = gsl_vector_int_alloc(n);
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

  /*
   * 1) apply recursive median filter to VFM residuals - to smooth noise and preserve jumps
   * 2) compute first differences of RM-filtered signal
   * 3) apply Hampel filter to first differences signal to locate spikes (indicating edges of jumps)
   */

  for (i = 0; i < 3; ++i)
    {
      gsl_filter_rmf(input[i], output_rmf[i], rmf_p);
      compute_diff(output_rmf[i], output_fd[i]);
      gsl_filter_hampel(nsigma, output_fd[i], output_hampel[i], xmedian[i], xsigma[i], &njump[i], ioutlier[i], hampel_p);
    }

  /* now correct B_VFM for detected jumps */
  for (i = 0; i < n; ++i)
    {
      size_t didx = i + start_idx;
      double *q = &(data->q[4 * didx]);
      double *B_VFM = &(data->B_VFM[3 * didx]);
      double B_model[4], B_model_VFM[3];
      size_t j;

      satdata_mag_model(didx, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_VFM);

      for (j = 0; j < 3; ++j)
        {
          int outlier = gsl_vector_int_get(ioutlier[j], i);

          if (outlier)
            {
              /* compute jump amplitude (difference between FD signal and Hampel-filtered FD signal) */
              double diff = gsl_vector_get(output_fd[j], i) - gsl_vector_get(output_hampel[j], i);
              double absdiff = fabs(diff);

              /* if jump amplitude satisfies minimum criteria, add it to running jump total */
              if (absdiff > min_jump_threshold)
                cumjump[j] += diff;
            }

          /* correct B_VFM data by subtracting total sum of all previous jumps found */
          B_VFM[j] -= cumjump[j];
        }

      if (fp)
        {
          fprintf(fp, "%ld %8.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
                  satdata_epoch2timet(data->t[didx]),
                  data->qdlat[didx],
                  gsl_vector_get(input[0], i),
                  gsl_vector_get(input[1], i),
                  gsl_vector_get(input[2], i),
                  gsl_vector_get(output_rmf[0], i),
                  gsl_vector_get(output_rmf[1], i),
                  gsl_vector_get(output_rmf[2], i),
                  gsl_vector_get(output_fd[0], i),
                  gsl_vector_get(output_fd[1], i),
                  gsl_vector_get(output_fd[2], i),
                  gsl_vector_get(output_hampel[0], i),
                  gsl_vector_get(output_hampel[1], i),
                  gsl_vector_get(output_hampel[2], i),
                  B_VFM[0] - B_model_VFM[0],
                  B_VFM[1] - B_model_VFM[1],
                  B_VFM[2] - B_model_VFM[2]);
        }
    }

  if (fp)
    fprintf(fp, "\n\n");

  for (i = 0; i < 3; ++i)
    {
      gsl_vector_free(input[i]);
      gsl_vector_free(output_rmf[i]);
      gsl_vector_free(output_fd[i]);
      gsl_vector_free(output_hampel[i]);
      gsl_vector_free(xmedian[i]);
      gsl_vector_free(xsigma[i]);
      gsl_vector_int_free(ioutlier[i]);
    }

  gsl_filter_rmf_free(rmf_p);
  gsl_filter_hampel_free(hampel_p);

  return 0;
}
