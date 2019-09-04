/*
 * crustal1.c
 *
 * Read NGDC-720 coefficients, use a cosine taper to force spectrum
 * to zero around n=300
 *
 * Usage: ./crustal1 [args]
 * [-b NGDC720_file]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_msynth.h>

/*
taper_gnm()
  Taper high degree gnm coefficients.
We use a cosine taper set to 1 at nmin, and going to 0 at
nmax + 1. Don't make it go to 0 at nmax or the spectrum will
suddenly drop many orders of magnitude between nmax-1 and nmax
*/

static int
taper_gnm(const size_t nmin, const size_t nmax, msynth_workspace *w)
{
  const double fac = M_PI / (2.0 * (nmax + 1 - nmin));
  size_t n;

  for (n = 1; n <= w->nmax; ++n)
    {
      int M = (int) n;
      int m;
      double wn = 1.0;

      /* compute taper weight */
      if (n > nmin)
        {
          double val = cos(fac * (n - nmin));
          wn = val * val;
        }

      for (m = -M; m <= M; ++m)
        {
          size_t cidx = msynth_nmidx(n, m, w);

          w->c[cidx] *= wn;
        }
    }

  return 0;
}

int
main(int argc, char *argv[])
{
  msynth_workspace *w1;
  msynth_workspace *w2;
  char *coeffile = "new_coef.txt";
  char *outfile_orig = "spectrum_orig.s";
  char *outfile = "spectrum.s";
  int c;
  const size_t nmin = 100;
  const size_t nmax = 300;
  double epoch = -1.0;

  fprintf(stderr, "main: reading %s...", MSYNTH_NGDC720_FILE);
  w1 = msynth_ngdc720_read(MSYNTH_NGDC720_FILE);
  fprintf(stderr, "done\n");

  w2 = msynth_copy(w1);

  if (epoch < 0.0)
    epoch = w1->epochs[0];

  fprintf(stderr, "main: epoch = %g\n", epoch);
  fprintf(stderr, "main: spectrum nmax = %zu\n", w1->eval_nmax);

  fprintf(stderr, "main: tapering spectrum...");
  taper_gnm(nmin, nmax, w2);
  w2->nmax = nmax;
  w2->eval_nmax = nmax;
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: printing original spectrum to %s...", outfile_orig);
  msynth_print_spectrum(outfile_orig, epoch, w1);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: printing spectrum to %s...", outfile);
  msynth_print_spectrum(outfile, epoch, w2);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing new coefficients to %s...", coeffile);
  msynth_crust_write(coeffile, w2);
  fprintf(stderr, "done\n");

  msynth_free(w1);
  msynth_free(w2);

  return 0;
}
