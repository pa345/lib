/*
 * window.c
 *
 * Window functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_blas.h>

/*
apply_ps1()

Inputs: in  - input vector (may be NULL, in which case output will contain
              window function)
        out - input multiplied by window function (if input is NULL, will
              contain only window function)

Return: success/error
*/

int
apply_ps1(const gsl_vector *in, gsl_vector *out)
{
  const size_t N = out->size;

  if (in != NULL && N != in->size)
    {
      GSL_ERROR("input and output vectors must be same length", GSL_EBADLEN);
    }
  else
    {
      const double d = sqrt(2.0 / (double)N);
      const double p[11] = { 5.3476939016920851e-11, 2.2654256220146656e-9,
                             7.8075102004229667e-8, 2.1373409644281953e-6,
                             4.5094847544714943e-5, 7.0498957221483167e-4,
                             7.7412693304064753e-3, 5.5280627452077586e-2,
                             2.2753754228751827e-1, 4.3433904277546202e-1,
                             2.2902051859068017e-1 };
      size_t i;

      for (i = 0; i < N; ++i)
        {
          double xt = (2.0*i-1.0)/(double)N - 1.0;
          double u = (1.0 - xt) * (1.0 + xt);
          double wi = ((((((((((p[0]*u+p[1])*u+p[2])*u+p[3])*u+p[4])*u+p[5])*u+p[6])*u+p[7])*u+p[8])*u+p[9])*u+p[10])*d;
          double *vi = gsl_vector_ptr(out, i);

          if (in != NULL)
            *vi = wi * gsl_vector_get(in, i);
          else
            *vi = wi;
        }

      return GSL_SUCCESS;
    }
}

int
apply_hamming(const gsl_vector *in, gsl_vector *out)
{
  const size_t N = out->size;

  if (in != NULL && N != in->size)
    {
      GSL_ERROR("input and output vectors must be same length", GSL_EBADLEN);
    }
  else
    {
      const double alpha = 0.53836;
      const double beta = 0.46164;
      size_t i;

      for (i = 0; i < N; ++i)
        {
          double ratio = 2.0 * M_PI * i / ((double)N - 1.0);
          double wi = alpha - beta*cos(ratio);
          double *vi = gsl_vector_ptr(out, i);

          if (in != NULL)
            *vi = wi * gsl_vector_get(in, i);
          else
            *vi = wi;
        }

      return GSL_SUCCESS;
    }
}

int
apply_kaiser(const gsl_vector *in, gsl_vector *out)
{
  const size_t N = out->size;

  if (in != NULL && N != in->size)
    {
      GSL_ERROR("input and output vectors must be same length", GSL_EBADLEN);
    }
  else
    {
      const double alpha = 3.0;
      const double beta = M_PI * alpha;
      const double d = gsl_sf_bessel_I0(beta);
      size_t i;

      for (i = 0; i < N; ++i)
        {
          double term1 = 2.0*i / (N - 1.0) - 1.0;
          double term2 = 1.0 - term1*term1;
          double wi = gsl_sf_bessel_I0(beta * sqrt(term2)) / d;
          double *vi = gsl_vector_ptr(out, i);

          if (in != NULL)
            *vi = wi * gsl_vector_get(in, i);
          else
            *vi = wi;
        }

      return GSL_SUCCESS;
    }
}

/* modified sin^2 window from Gary; constructed so that sum(window) = N/2 */
int
apply_modsinsq(const gsl_vector *in, gsl_vector *out)
{
  const size_t N = out->size;

  if (in != NULL && N != in->size)
    {
      GSL_ERROR("input and output vectors must be same length", GSL_EBADLEN);
    }
  else
    {
      size_t i;
      gsl_vector * a = gsl_vector_alloc(N);

      for (i = 0; i < N / 2; ++i)
        {
          double xi = (i + 1.0) * M_PI / (double) N;
          double sxi = sin(xi);
          double wi = sxi * sxi;

          gsl_vector_set(out, i, wi);
          gsl_vector_set(out, N - i - 1, wi);
        }

      for (i = 0; i < N / 2; ++i)
        {
          double wi = gsl_vector_get(out, i);
          double ai = gsl_vector_get(out, N / 2 + i);

          gsl_vector_set(a, i, ai);
          gsl_vector_set(a, N / 2 + i, wi);
        }

      for (i = 0; i < N / 2; ++i)
        {
          double * wi = gsl_vector_ptr(out, i);
          double * wNi = gsl_vector_ptr(out, N - i - 1);
          double ai = gsl_vector_get(a, i);
          double scale = *wi + ai;

          *wi /= scale;
          *wNi /= scale;
        }

      if (in != NULL)
        {
          for (i = 0; i < N; ++i)
            {
              double in_i = gsl_vector_get(in, i);
              double * out_i = gsl_vector_ptr(out, i);
              *out_i *= in_i;
            }
        }

      gsl_vector_free(a);

      return GSL_SUCCESS;
    }
}
