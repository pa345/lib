#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>

#include "yuen.c"

int
run_trial(const double trim, const size_t nx, const size_t ny, double *x, double *y, int * status1, int * status2, gsl_rng * r)
{
  const double alpha = 0.001;
  yuen_stats_type ystats;
  medtest_stats_type mstats;
  size_t i;

  for (i = 0; i < nx; ++i)
    x[i] = gsl_ran_gaussian(r, 2.0);

  for (i = 0; i < ny; ++i)
    y[i] = gsl_ran_gaussian(r, 2.0) + 7.0;

  *status1 = yuen(trim, alpha, x, 1, nx, y, 1, ny, &ystats);
  *status2 = medtest(alpha, x, 1, nx, y, 1, ny, &mstats);

#if 0
  printf("==== YUEN TRIMMED MEAN TEST ====\n");

  printf("CI = [%g,%g] p-value = %g\n", ystats.lower, ystats.upper, ystats.p);
  printf("mean1 = %g mean2 = %g\n", ystats.mean1, ystats.mean2);

  if (status == GSL_SUCCESS)
    printf("datasets have EQUAL trimmed means\n");
  else
    printf("datasets have UNEQUAL trimmed means\n");

  printf("==== MEDIAN/QN TEST ====\n");

  status = medtest(alpha, x, 1, nx, y, 1, ny, &mstats);

  printf("CI = [%g,%g] p-value = %g\n", mstats.lower, mstats.upper, mstats.p);
  printf("median1 = %g median2 = %g\n", mstats.median1, mstats.median2);

  if (status == GSL_SUCCESS)
    printf("datasets have EQUAL medians\n");
  else
    printf("datasets have UNEQUAL medians\n");
#endif

  return 0;
}

int
main(int argc, char *argv[])
{
  int status1, status2;
  const size_t ntrial = 10000;
  const size_t nx = 10;
  const size_t ny = 10;
  double trim = 0.2;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
  double *x = malloc(nx * sizeof(double));
  double *y = malloc(ny * sizeof(double));
  size_t i;
  size_t nyuen = 0, nmedian = 0;

  if (argc > 1)
    trim = atof(argv[1]);

  gsl_rng_set(r, time(NULL));

  for (i = 0; i < ntrial; ++i)
    {
      run_trial(trim, nx, ny, x, y, &status1, &status2, r);

      if (status1 == GSL_SUCCESS)
        ++nyuen;

      if (status2 == GSL_SUCCESS)
        ++nmedian;
    }

  fprintf(stderr, "Yuen:      successes: %zu/%zu (%.1f%%)\n", nyuen, ntrial, (double)nyuen / (double)ntrial * 100.0);
  fprintf(stderr, "Median/Qn: successes: %zu/%zu (%.1f%%)\n", nmedian, ntrial, (double)nmedian / (double)ntrial * 100.0);

  gsl_rng_free(r);
  free(x);
  free(y);

  return 0;
}
