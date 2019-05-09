/*
 * task2.c
 *
 * WMM Simulations Task 2
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <bspline2/gsl_bspline2.h>
#include <msynth/msynth.h>

/* subtract 3 years from timestamps, so that range (2014.5, 2017.5) becomes (2011.5, 2014.5), similar to WMM modeling timeframe */
static int
task2_times(mfield_data_workspace * w)
{
  size_t i, j;

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_ExistVector(mptr->flags[j]) ||
              (mptr->flags[j] & MAGDATA_FLG_Z))
            {
              long year, month, day, hour, min, sec, msec;

              EPOCHbreakdown(mptr->t[j], &year, &month, &day, &hour, &min, &sec, &msec);

              year = year - 3;

              mptr->t[j] = computeEPOCH(year, month, day, hour, min, sec, msec);
              mptr->ts[j] = mptr->t[j];
            }
          else
            {
              mptr->flags[j] |= MAGDATA_FLG_DISCARD;
            }
        }
    }

  return 0;
}

/* replace B_nec with CHAOS values */
static int
task2_synth(mfield_data_workspace * w)
{
  msynth_workspace * msynth_p = msynth_shc_read(MSYNTH_CHAOS_FILE);
  size_t i, j;

  msynth_set(1, 15, msynth_p);

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* if no vector measurement, do nothing */
          if (MAGDATA_ExistVector(mptr->flags[j]) ||
              (mptr->flags[j] & MAGDATA_FLG_Z))
            {
              double r = mptr->r[j];
              double theta = mptr->theta[j];
              double phi = mptr->phi[j];
              const double tyr = epoch2year(mptr->t[j]);
              double B[4];

              msynth_eval(tyr, r, theta, phi, B, msynth_p);

              mptr->Bx_nec[j] = B[0];
              mptr->By_nec[j] = B[1];
              mptr->Bz_nec[j] = B[2];
              mptr->F[j] = B[3];

              mptr->Bx_model[j] = 0.0;
              mptr->By_model[j] = 0.0;
              mptr->Bz_model[j] = 0.0;
            }
        }
    }
    

  msynth_free(msynth_p);

  return 0;
}

/*
task2_add_noise()
  For WMM simulation studies - add random Gaussian noise to vector
data (b1) in NEC frame

Inputs: sigma - standard deviation of noise (nT);
                if negative, not used
        w     - workspace
*/

static int
task2_add_b1(const double sigma, mfield_data_workspace * w)
{
  int s = 0;
  size_t i, j;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* if no vector measurement, do nothing */
          if (MAGDATA_ExistVector(mptr->flags[j]) ||
              (mptr->flags[j] & MAGDATA_FLG_Z))
            {
              if (sigma > 0.0)
                {
                  mptr->Bx_nec[j] += gsl_ran_gaussian(r, sigma);
                  mptr->By_nec[j] += gsl_ran_gaussian(r, sigma);
                  mptr->Bz_nec[j] += gsl_ran_gaussian(r, sigma);

                  /* recompute scalar field measurement */
                  mptr->F[j] = gsl_hypot3(mptr->Bx_nec[j], mptr->By_nec[j], mptr->Bz_nec[j]);
                }
            }
        }
    }

  gsl_rng_free(r);

  return s;
}

static int
task2_add_b2(mfield_data_workspace * w)
{
  const size_t ncontrol = 22; /* number of B-spline coefficients */
  const size_t k = 4;         /* B-spline order */
  gsl_bspline2_workspace * bspline_r = gsl_bspline2_alloc_ncontrol(k, ncontrol);
  gsl_bspline2_workspace * bspline_t = gsl_bspline2_alloc_ncontrol(k, ncontrol);
  gsl_bspline2_workspace * bspline_p = gsl_bspline2_alloc_ncontrol(k, ncontrol);
  gsl_vector * coef_r = gsl_vector_alloc(ncontrol);
  gsl_vector * coef_t = gsl_vector_alloc(ncontrol);
  gsl_vector * coef_p = gsl_vector_alloc(ncontrol);
  size_t i, j;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  FILE *fp;

  gsl_bspline2_init_uniform(0.0, 90.0, bspline_r);
  gsl_bspline2_init_uniform(0.0, 90.0, bspline_t);
  gsl_bspline2_init_uniform(0.0, 90.0, bspline_p);

  fp = fopen("/data/palken/matlab/Arnaud/coef_r.txt", "r");
  gsl_vector_fscanf(fp, coef_r);
  fclose(fp);

  fp = fopen("/data/palken/matlab/Arnaud/coef_t.txt", "r");
  gsl_vector_fscanf(fp, coef_t);
  fclose(fp);

  fp = fopen("/data/palken/matlab/Arnaud/coef_p.txt", "r");
  gsl_vector_fscanf(fp, coef_p);
  fclose(fp);

  for (i = 0; i < w->nsources; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* if no vector measurement, do nothing */
          if (MAGDATA_ExistVector(mptr->flags[j]) ||
              (mptr->flags[j] & MAGDATA_FLG_Z))
            {
              double sigma_r, sigma_t, sigma_p;
              double qdlat = fabs(mptr->qdlat[j]);

              gsl_bspline2_eval(qdlat, coef_r, &sigma_r, bspline_r);
              gsl_bspline2_eval(qdlat, coef_t, &sigma_t, bspline_t);
              gsl_bspline2_eval(qdlat, coef_p, &sigma_p, bspline_p);

              mptr->Bx_nec[j] += gsl_ran_gaussian(r, sigma_t);
              mptr->By_nec[j] += gsl_ran_gaussian(r, sigma_p);
              mptr->Bz_nec[j] += gsl_ran_gaussian(r, sigma_r);

              /* recompute scalar field measurement */
              mptr->F[j] = gsl_hypot3(mptr->Bx_nec[j], mptr->By_nec[j], mptr->Bz_nec[j]);
            }
        }
    }

  gsl_bspline2_free(bspline_r);
  gsl_bspline2_free(bspline_t);
  gsl_bspline2_free(bspline_p);
  gsl_vector_free(coef_r);
  gsl_vector_free(coef_t);
  gsl_vector_free(coef_p);
  gsl_rng_free(r);

  return 0;
}

static int
task2_add_noise(const double sigma, mfield_data_workspace * w)
{
  int s = 0;

  fprintf(stderr, "task2_add_noise: recalculating timestamps...");
  task2_times(w);
  fprintf(stderr, "done\n");

  fprintf(stderr, "task2_add_noise: synthing CHAOS values...");
  task2_synth(w);
  fprintf(stderr, "done\n");

  if (sigma > 0.0)
    task2_add_b1(sigma, w);

  task2_add_b2(w);

  return s;
}
