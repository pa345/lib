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

/*
calc_nlm_complex()
  Calculate number of complex coefficients, including negative and positive m values

Inputs: lmin - minimum SH degree
        lmax - maximum SH degree
        mmax - maximum SH order

Return: number of SH coefficients
*/

static size_t
calc_nlm_complex(const size_t lmin, const size_t lmax, const size_t mmax)
{
  size_t nlm = (mmax + lmin + 1) * (mmax - lmin + 1) + (2 * mmax + 1) * (lmax - mmax);
  return nlm;
}

static size_t
lmidx_complex(const size_t l, const int m, const size_t mmax)
{
  size_t idx;

  if (l <= mmax)
    {
      size_t base = (l - 1) * (l + 1);
      int offset = m + (int) l;
      idx = base + (size_t) offset;
    }
  else
    {
      size_t base1 = mmax * (mmax + 2);
      size_t base2 = (l - mmax - 1) * (2 * mmax + 1);
      int offset = m + (int) mmax;
      idx = base1 + base2 + (size_t) offset;
    }

  return idx;
}
