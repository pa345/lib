#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

typedef struct
{
  double p;
  double lower;
  double upper;
  double mean1;
  double mean2;
} yuen_stats_type;

typedef struct
{
  double p;
  double lower;
  double upper;
  double median1;
  double median2;
  double tstat;
  double se;
} medtest_stats_type;

/* Winsorized variance */
double
winvar_from_sorted_data(const double trim, const double sorted_data[],
                        const size_t stride, const size_t n)
{
  if (trim >= 0.5)
    {
      /* entire data set is trimmed */
      return (0.0);
    }
  else
    {
      size_t ilow = (size_t) floor(trim * n);
      size_t ihigh = n - ilow - 1;
      double xlow = sorted_data[ilow * stride];
      double xhigh = sorted_data[ihigh * stride];
      double mean = 0.0;
      double M2 = 0.0;
      double k = 0.0;
      size_t i;

      /* process samples on left tail (equal to xlow) */
      for (i = 0; i < ilow; ++i)
        {
          double delta = xlow - mean;
          k += 1.0;
          mean += delta / k;
          M2 += delta * (xlow - mean);
        }

      /* process samples in middle (not modified) */
      for (i = ilow; i <= ihigh; ++i)
        {
          double delta = sorted_data[i * stride] - mean;
          k += 1.0;
          mean += delta / k;
          M2 += delta * (sorted_data[i * stride] - mean);
        }

      /* process samples on right tail (equal to xhigh) */
      for (i = ihigh + 1; i < n; ++i)
        {
          double delta = xhigh - mean;
          k += 1.0;
          mean += delta / k;
          M2 += delta * (xhigh - mean);
        }

      return (M2 / (k - 1.0));
    }
}

int
yuen_from_sorted_data(const double trim, const double alpha,
                      const double xsorted_data[], const size_t xstride, const size_t xlen,
                      const double ysorted_data[], const size_t ystride, const size_t ylen,
                      yuen_stats_type * stats)
{
  if (trim >= 0.5)
    {
      GSL_ERROR("trim factor must be less than 0.5", GSL_EDOM);
    }
  else
    {
      double h1 = (double) xlen - 2.0 * floor(trim * xlen);
      double h2 = (double) ylen - 2.0 * floor(trim * ylen);
      double q1 = (xlen - 1.0) * winvar_from_sorted_data(trim, xsorted_data, xstride, xlen) / (h1 * (h1 - 1.0)); /* squared standard error of x */
      double q2 = (ylen - 1.0) * winvar_from_sorted_data(trim, ysorted_data, ystride, ylen) / (h2 * (h2 - 1.0)); /* squared standard error of y */
      double df = (q1 + q2) * (q1 + q2) / (q1*q1/(h1 - 1.0) + q2*q2/(h2 - 1.0));                                 /* degrees of freedom */
      double mean1 = gsl_stats_trmean_from_sorted_data(trim, xsorted_data, xstride, xlen);                       /* trmean(x) */
      double mean2 = gsl_stats_trmean_from_sorted_data(trim, ysorted_data, ystride, ylen);                       /* trmean(y) */
      double difference = mean1 - mean2;                                                                         /* numerator of test statistic (difference of trimmed means) */
      double se = sqrt(q1 + q2);                                                                                 /* denominator of test statistic (combined standard error) */
      double ystat = fabs(difference / se);                                                                      /* Yuen test statistic */
      double pval = 2.0 * (1.0 - gsl_cdf_tdist_P(ystat, df));                                                    /* p-value */

      if (stats != NULL)
        {
          double spread = gsl_cdf_tdist_Pinv(1.0 - 0.5*alpha, df) * se;

          stats->lower = difference - spread; /* lower confidence bound */
          stats->upper = difference + spread; /* upper confidence bound */
          stats->p = pval;
          stats->mean1 = mean1;
          stats->mean2 = mean2;
        }

      if (pval <= alpha)
        return GSL_FAILURE;
      else
        return GSL_SUCCESS;
    }
}

int
yuen(const double trim, const double alpha,
     double xdata[], const size_t xstride, const size_t xlen,
     double ydata[], const size_t ystride, const size_t ylen,
     yuen_stats_type * stats)
{
  gsl_sort(xdata, xstride, xlen);
  gsl_sort(ydata, ystride, ylen);

  return yuen_from_sorted_data(trim, alpha,
                               xdata, xstride, xlen,
                               ydata, ystride, ylen, stats);
}

int
medtest_from_sorted_data(const double alpha,
                         const double xsorted_data[], const size_t xstride, const size_t xlen,
                         const double ysorted_data[], const size_t ystride, const size_t ylen,
                         medtest_stats_type * stats)
{
  double *work = malloc(3 * GSL_MAX(xlen, ylen) * sizeof(double));
  int *work_int = malloc(5 * GSL_MAX(xlen, ylen) * sizeof(int));
  double Q1 = gsl_stats_Qn_from_sorted_data(xsorted_data, xstride, xlen, work, work_int); /* Q_n(x) */
  double Q2 = gsl_stats_Qn_from_sorted_data(ysorted_data, ystride, ylen, work, work_int); /* Q_n(y) */
  double median1 = gsl_stats_median_from_sorted_data(xsorted_data, xstride, xlen);        /* median(x) */
  double median2 = gsl_stats_median_from_sorted_data(ysorted_data, ystride, ylen);        /* median(y) */
  double difference = median1 - median2;                                                  /* numerator of test statistic (difference of medians) */
  double se = sqrt(M_PI_2 * (Q1 * Q1 / xlen + Q2 * Q2 / ylen));                           /* denominator of test statistic (combined standard error) */
  double tstat = fabs(difference / se);                                                   /* test statistic */
  double pval = 2.0 * (1.0 - gsl_cdf_ugaussian_P(tstat));                                 /* p-value */

  if (stats != NULL)
    {
      double spread = gsl_cdf_ugaussian_Pinv(1.0 - 0.5*alpha) * se;

      stats->lower = difference - spread; /* lower confidence bound */
      stats->upper = difference + spread; /* upper confidence bound */
      stats->p = pval;
      stats->median1 = median1;
      stats->median2 = median2;
      stats->tstat = tstat;
      stats->se = se;
    }

  free(work);
  free(work_int);

  if (pval <= alpha)
    return GSL_FAILURE;
  else
    return GSL_SUCCESS;
}

int
medtest(const double alpha,
        double xdata[], const size_t xstride, const size_t xlen,
        double ydata[], const size_t ystride, const size_t ylen,
        medtest_stats_type * stats)
{
  gsl_sort(xdata, xstride, xlen);
  gsl_sort(ydata, ystride, ylen);

  return medtest_from_sorted_data(alpha, xdata, xstride, xlen,
                                  ydata, ystride, ylen, stats);
}

/* t-test for medians, assuming same variance in both windows */
int
medtest2_from_sorted_data(const double alpha,
                          const double xsorted_data[], const size_t xstride, const size_t xlen,
                          const double ysorted_data[], const size_t ystride, const size_t ylen,
                          medtest_stats_type * stats)
{
  double *z = malloc((xlen + ylen) * sizeof(double));
  double *work = malloc(3 * (xlen + ylen) * sizeof(double));
  int *work_int = malloc(5 * (xlen + ylen) * sizeof(int));
  double median_x = gsl_stats_median_from_sorted_data(xsorted_data, xstride, xlen);   /* median(x) */
  double median_y = gsl_stats_median_from_sorted_data(ysorted_data, ystride, ylen);   /* median(y) */
  double difference = median_x - median_y;                                            /* numerator of test statistic (difference of medians) */
  double se;                                                                          /* denominator of test statistic (combined standard error) */
  double tstat;                                                                       /* test statistic */
  double pval;                                                                        /* p-value */
  double Qn;                                                                          /* Q_n of combined window x and y */
  size_t idx = 0;
  size_t i;

  for (i = 0; i < xlen; ++i)
    z[idx++] = xsorted_data[i * xstride] - median_x;

  for (i = 0; i < ylen; ++i)
    z[idx++] = ysorted_data[i * ystride] - median_y;

  gsl_sort(z, 1, xlen + ylen);
  Qn = gsl_stats_Qn_from_sorted_data(z, 1, xlen + ylen, work, work_int);

  se = sqrt(M_PI_2 * Qn * Qn * (1.0 / (double)xlen + 1.0 / (double)ylen));
  tstat = fabs(difference / se);
  pval = 2.0 * (1.0 - gsl_cdf_ugaussian_P(tstat));

  if (stats != NULL)
    {
      double spread = gsl_cdf_ugaussian_Pinv(1.0 - 0.5*alpha) * se;

      stats->lower = difference - spread; /* lower confidence bound */
      stats->upper = difference + spread; /* upper confidence bound */
      stats->p = pval;
      stats->median1 = median_x;
      stats->median2 = median_y;
      stats->tstat = tstat;
      stats->se = se;
    }

  free(z);
  free(work);
  free(work_int);

  if (pval <= alpha)
    return GSL_FAILURE;
  else
    return GSL_SUCCESS;
}

int
medtest2(const double alpha,
         double xdata[], const size_t xstride, const size_t xlen,
         double ydata[], const size_t ystride, const size_t ylen,
         medtest_stats_type * stats)
{
  gsl_sort(xdata, xstride, xlen);
  gsl_sort(ydata, ystride, ylen);

  return medtest2_from_sorted_data(alpha, xdata, xstride, xlen,
                                   ydata, ystride, ylen, stats);
}
