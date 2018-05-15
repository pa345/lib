/*
 * common.c
 */

/*
count_windows()
  Count number of time segment windows in FFT analysis. So
for a sliding 2 day window with a 1 day overlap, set

window_size = 2
window_shift = 1

Inputs: nsamples     - total number of samples in time series
        fs           - sampling frequency in 1/days
        window_size  - number of days per window
        window_shift - number of days to advance/slide forward

Return: total number of windows
*/

static size_t
count_windows(const size_t nsamples, const double fs,
              const double window_size, const double window_shift)
{
  const size_t nwindow = (size_t) (window_size * fs);   /* number of samples per window */
  const size_t nforward = (size_t) (window_shift * fs); /* number of samples to slide forward */
  size_t T = 0;                                         /* number of windows */
  size_t end_idx = nwindow - 1;                         /* first window contains samples [0, nwindow - 1] */

  while (end_idx < nsamples)
    {
      ++T;
      end_idx += nforward;
    }

  return T;
}
