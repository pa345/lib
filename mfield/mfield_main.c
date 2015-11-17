/*
 * mfield_main.c
 *
 * Usage:
 * ./mfield [flags]
 *
 * Flags:
 *   -c coef_output_file
 *   -n max_iterations
 *   -e epoch_decimal_year
 *   -v lambda_sv
 *   -a lambda_sa
 *   -p euler_period_days
 *   -r residual_file
 *   -l Lcurve_data_file
 *
 * After each iteration, the file 'res.#.dat' is written
 * where # is the iteration number. This file contains the
 * residuals of a sample of the DMSP dataset.
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>

#include <satdata/satdata.h>

#include "common.h"
#include "euler.h"
#include "oct.h"
#include "magdata.h"
#include "mfield.h"
#include "mfield_test.h"
#include "msynth.h"
#include "track.h"

/* maximum spherical harmonic degree (internal) */
#define NMAX_MF              5
#define NMAX_SV              5
#define NMAX_SA              12

#define MAX_BUFFER           2048

#define MFIELD_CALC_UNCERTAINTIES   0

/* n m gnm dgnm ddgnm */
static mfield_test_coeff test_gnm[] = {
  /*
   * need to start somewhat close to IGRF for main field, otherwise
   * Euler angles won't converge
   */
  { 1, 0, -30000.0, -25.0, 1.2 },
  { 1, -1, 5000.0, 30.0, -3.2 },
  { 2, 1, 3000.0, -4.0, 4.2 },
  { 3, -2, 250.0, -3.0, -6.2 },

  { 0, 0, 0.0, 0.0, 0.0 }
};

/* replace with synthetic data for testing */
int
replace_synthetic_data(mfield_workspace *w)
{
  const size_t nmax = w->nmax_mf;
  const double t0 = w->epoch;
  size_t i, j;
  size_t plm_size = gsl_sf_legendre_array_n(nmax);
  double *Plm = malloc(plm_size * sizeof(double));
  double *dPlm = malloc(plm_size * sizeof(double));
  size_t c_size = w->p;
  gsl_vector *g = gsl_vector_calloc(c_size);
  mfield_test_coeff *gptr;

#if MFIELD_FIT_EULER
  /* Euler angles */
  const double alpha = -13.1 * M_PI / 180.0 * 1.0;
  const double beta = -5.2 * M_PI / 180.0 * 1.0;
  const double gamma = 3.4 * M_PI / 180.0 * 1.0;
#endif

  /* initialize coefficients g from test_gnm[] */
  for (gptr = &test_gnm[0]; gptr->n != 0; ++gptr)
    {
      size_t n = gptr->n;
      int m = gptr->m;
      size_t cidx = mfield_coeff_nmidx(n, m);

      mfield_set_mf(g, cidx, gptr->gnm, w);
      mfield_set_sv(g, cidx, gptr->dgnm, w);
      mfield_set_sa(g, cidx, gptr->ddgnm, w);
    }

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);

      for (j = 0; j < mptr->n; ++j)
        {
          double t = satdata_epoch2year(mptr->t[j]);
          double t1 = t - t0;
          double t2 = 0.5 * t1 * t1;
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          double sint = sin(theta);
          double cost = cos(theta);
          double X = 0.0, Y = 0.0, Z = 0.0;
          double ratio = w->R / r;
          size_t n;
          int m;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                          cost, Plm, dPlm);

          for (n = 1; n <= nmax; ++n)
            {
              int ni = (int) n;
              double rterm = pow(ratio, n + 2.0);

              for (m = 0; m <= ni; ++m)
                {
                  double c = cos(m * phi);
                  double s = sin(m * phi);
                  size_t pidx = gsl_sf_legendre_array_index(n, m);
                  size_t cidx = mfield_coeff_nmidx(n, m);
                  double gnm, hnm = 0.0;

                  gnm = mfield_get_mf(g, cidx, w) +
                        mfield_get_sv(g, cidx, w) * t1 +
                        mfield_get_sa(g, cidx, w) * t2;

                  if (m > 0)
                    {
                      cidx = mfield_coeff_nmidx(n, -m);
                      hnm = mfield_get_mf(g, cidx, w) +
                            mfield_get_sv(g, cidx, w) * t1 +
                            mfield_get_sa(g, cidx, w) * t2;
                    }

                  X += rterm * (gnm * c + hnm * s) * dPlm[pidx];
                  Y += rterm / sint * m * (gnm * s - hnm * c) * Plm[pidx];
                  Z -= (n + 1.0) * rterm *
                       (gnm * c + hnm * s) * Plm[pidx];
                }
            }

          mptr->Bx_nec[j] = X;
          mptr->By_nec[j] = Y;
          mptr->Bz_nec[j] = Z;
          mptr->F[j] = gsl_hypot3(X, Y, Z);

          /* no crustal/external field for synthetic data */
          mptr->Bx_model[j] = 0.0;
          mptr->By_model[j] = 0.0;
          mptr->Bz_model[j] = 0.0;

#if MFIELD_FIT_EULER
          /* rotate NEC vector to VFM frame */
          {
            double *q = &(mptr->q[4*j]);
            double B_nec[3], B_vfm[3];

            B_nec[0] = X;
            B_nec[1] = Y;
            B_nec[2] = Z;

            euler_nec2vfm(EULER_FLG_ZYX, alpha, beta, gamma, q, B_nec, B_vfm);

            mptr->Bx_vfm[j] = B_vfm[0];
            mptr->By_vfm[j] = B_vfm[1];
            mptr->Bz_vfm[j] = B_vfm[2];
          }
#endif
        }
    }

  free(Plm);
  free(dPlm);
  gsl_vector_free(g);

  return 0;
} /* replace_synthetic_data() */

#if 0

/*
parse_input()
  Read in index file for a given satellite, perform data
selection and downsampling
*/

satdata_mag *
parse_input(const size_t sat_idx)
{
  satdata_mag *data = NULL;
  const char *idxfile = index_files[sat_idx];
  size_t i;

  fprintf(stderr, "parse_input: reading %s...", idxfile);

  if (sat_idx >= IDX_SWA && sat_idx <= IDX_SWC)
    data = satdata_swarm_read_idx(idxfile, 0);
  else if (sat_idx == IDX_CHAMP)
    data = satdata_champ_read_idx(idxfile, 0);
  else if (sat_idx == IDX_OERSTED)
    data = satdata_oersted_read_idx(idxfile, 0);
  else if (sat_idx >= IDX_F15 && sat_idx <= IDX_F18)
    data = satdata_dmsp_read_idx(idxfile, 0);

  if (!data)
    return NULL;

  fprintf(stderr, "done (%zu points read)\n", data->n);

  if (sat_idx >= IDX_SWA && sat_idx <= IDX_SWC)
    {
      size_t nrms;
      track_workspace *track_p = track_alloc();
      double thresh[] = { 20.0, 25.0, 15.0, 15.0 };

      satdata_swarm_filter_instrument(1, data);

      /* filter by track rms */

      track_init(data, NULL, track_p);

      nrms = track_flag_rms("swarm_rms.dat", thresh, data, track_p);
      fprintf(stderr, "parse_input: flagged (%zu/%zu) (%.1f%%) points due to high rms\n",
              nrms, data->n, (double) nrms / (double) data->n * 100.0);

      track_free(track_p);

      satdata_filter_wmm(1, data);
    }

  fprintf(stderr, "parse_input: downsampling data by factor %d...", DOWNSAMPLE);

  for (i = 0; i < data->n; ++i)
    {
      if (i % DOWNSAMPLE != 0)
        data->flags[i] |= SATDATA_FLG_OUTLIER;
    }

  fprintf(stderr, "done\n");

  /* flag local time */
  if (!(sat_idx >= IDX_F15 && sat_idx <= IDX_F18))
    {
      size_t nlt;
      const double lt_min = 5.0;
      double lt_max = 22.0;
      const double euler_lt_min = 6.0;
      const double euler_lt_max = 18.0;

      /* in the first half of 2013, Oersted is in a ~10am/10pm orbit */
      if (sat_idx == IDX_OERSTED)
        lt_max = 20.0;

      fprintf(stderr, "parse_input: flagging points inside LT window [%g,%g], euler [%g,%g]...",
              lt_min, lt_max, euler_lt_min, euler_lt_max);

      nlt = flag_local_time(lt_min, lt_max, euler_lt_min, euler_lt_max, data);

      fprintf(stderr, "done (%zu/%zu data flagged)\n", nlt, data->n);
    }

  {
    size_t nflagged = satdata_nflagged(data);
    fprintf(stderr, "parse_input: total flagged points: %zu/%zu (%.1f%%) (%zu remaining)\n",
            nflagged, data->n, (double)nflagged / (double)data->n * 100.0,
            data->n - nflagged);
  }

  return data;
} /* parse_input() */

#endif /* 0 */

/*
initial_guess()
  Construct initial guess for main field coefficients. These
are based on the relevant IGRF coefficients, extrapolated forward
to the desired epoch using the SV coefficients. Initial SA coefficients
are set to 0.
*/

int
initial_guess(gsl_vector *c, mfield_workspace *w)
{
  msynth_workspace *msynth_p = msynth_igrf_read(MSYNTH_IGRF_FILE);
  const size_t nmax = GSL_MIN(w->nmax_mf, msynth_p->nmax);
  const double t = w->epoch;                       /* desired epoch */
  const double t0 = msynth_get_epoch(t, msynth_p); /* IGRF epoch */
  const double dt = t - t0;
  size_t n;
  int m;

  gsl_vector_set_zero(c);

  for (n = 1; n <= nmax; ++n)
    {
      int ni = (int) n;

      for (m = -ni; m <= ni; ++m)
        {
          size_t midx = msynth_nmidx(n, m, msynth_p);
          size_t cidx = mfield_coeff_nmidx(n, m);
          double gnm = msynth_get_mf(t, midx, msynth_p);
          double dgnm = msynth_get_sv(t, midx, msynth_p);

          /*
           * use SV prediction to update main field coefficients for new
           * epoch
           */
          mfield_set_mf(c, cidx, gnm + dt * dgnm, w);
          mfield_set_sv(c, cidx, dgnm, w);
          mfield_set_sa(c, cidx, 0.0, w);
        }
    }

  msynth_free(msynth_p);

  return 0;
} /* initial_guess() */

int
print_spectrum(const char *filename, mfield_workspace *w)
{
  size_t n;
  FILE *fp = fopen(filename, "w");

  fprintf(stderr, "print_spectrum: writing spectrum to %s...", filename);
  for (n = 1; n <= w->nmax_mf; ++n)
    {
      double gn = mfield_spectrum(n, w);
      double dgn = mfield_spectrum_sv(n, w);
      double ddgn = mfield_spectrum_sa(n, w);
      fprintf(fp, "%zu %.12e %.12e %.12e\n", n, gn, dgn, ddgn);
    }
  fprintf(stderr, "done\n");

  fclose(fp);

  return 0;
} /* print_spectrum() */

/*
print_residuals()
  Output residuals for each satellite, using data stored
in w->data_workspace_p and coefficients stored in w->c
*/

int
print_residuals(const char *filename, mfield_workspace *w)
{
  size_t i, j, k;

  for (i = 0; i < w->nsat; ++i)
    {
      const size_t nbins = 500;
      const double a = -100.0;
      const double b = 100.0;
      magdata *mptr = mfield_data_ptr(i, w->data_workspace_p);
      FILE *fp_res, *fp_hist;
      char fileres[2048];
      char filehist[2048];
      gsl_histogram *hf, *hz;
      gsl_rstat_workspace *rm_f = gsl_rstat_alloc();

      hf = gsl_histogram_alloc(nbins);
      hz = gsl_histogram_alloc(nbins);
      gsl_histogram_set_ranges_uniform(hf, a, b);
      gsl_histogram_set_ranges_uniform(hz, a, b);

      if (mptr->n == 0)
        continue;

      sprintf(fileres, "%s.sat%zu", filename, i);
      fp_res = fopen(fileres, "w");
      if (!fp_res)
        {
          fprintf(stderr, "print_residuals: unable to open %s: %s\n",
                  fileres, strerror(errno));
          return GSL_FAILURE;
        }

      sprintf(filehist, "%s.sat%zu.hist", filename, i);
      fp_hist = fopen(filehist, "w");
      if (!fp_hist)
        {
          fprintf(stderr, "print_residuals: unable to open %s: %s\n",
                  filehist, strerror(errno));
          return GSL_FAILURE;
        }

      k = 1;
      fprintf(fp_res, "# Field %zu: time (years)\n", k++);
      fprintf(fp_res, "# Field %zu: local time (hours)\n", k++);
      fprintf(fp_res, "# Field %zu: altitude (km)\n", k++);
      fprintf(fp_res, "# Field %zu: longitude (deg)\n", k++);
      fprintf(fp_res, "# Field %zu: latitude (deg)\n", k++);
      fprintf(fp_res, "# Field %zu: QD latitude (deg)\n", k++);
      fprintf(fp_res, "# Field %zu: scalar residual (nT)\n", k++);
      fprintf(fp_res, "# Field %zu: X residual (nT)\n", k++);
      fprintf(fp_res, "# Field %zu: Y residual (nT)\n", k++);
      fprintf(fp_res, "# Field %zu: Z residual (nT)\n", k++);
      fprintf(fp_res, "# Field %zu: NEC X component (nT)\n", k++);
      fprintf(fp_res, "# Field %zu: NEC Y component (nT)\n", k++);
      fprintf(fp_res, "# Field %zu: NEC Z component (nT)\n", k++);
      fprintf(fp_res, "# Field %zu: satellite direction (+1 north -1 south)\n", k++);
      fprintf(fp_res, "# Field %zu: scalar data used in MF fitting (1 or 0)\n", k++);
      fprintf(fp_res, "# Field %zu: vector data used in MF fitting (1 or 0)\n", k++);
      fprintf(fp_res, "# Field %zu: vector data used in Euler angle fitting (1 or 0)\n", k++);
      fprintf(fp_res, "# Field %zu: along-track gradient available (1 or 0)\n", k++);

      fprintf(stderr, "Writing residuals to %s...", fileres);

      for (j = 0; j < mptr->n; ++j)
        {
          double t = satdata_epoch2year(mptr->t[j]);
          double r = mptr->r[j];
          double theta = mptr->theta[j];
          double phi = mptr->phi[j];
          time_t unix_time = satdata_epoch2timet(mptr->t[j]);
          double lt = get_localtime(unix_time, phi);
          double B_nec[3], B_int[4], B_ext[4], B_model[4];
          double res[4];

          mfield_eval(mptr->t[j], r, theta, phi, B_int, w);
          mfield_eval_ext(mptr->t[j], r, theta, phi, B_ext, w);

          /* add apriori external/crustal fields to computed external field */
          B_ext[0] += mptr->Bx_model[j];
          B_ext[1] += mptr->By_model[j];
          B_ext[2] += mptr->Bz_model[j];

#if MFIELD_FIT_EULER
          if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
            {
              size_t euler_idx = mfield_euler_idx(i, mptr->t[j], w);
              double alpha = gsl_vector_get(w->c, euler_idx);
              double beta = gsl_vector_get(w->c, euler_idx + 1);
              double gamma = gsl_vector_get(w->c, euler_idx + 2);
              double *q = &(mptr->q[4*j]);

              B_nec[0] = mptr->Bx_vfm[j];
              B_nec[1] = mptr->By_vfm[j];
              B_nec[2] = mptr->Bz_vfm[j];

              /* rotate to NEC with computed Euler angles */
              euler_vfm2nec(EULER_FLG_ZYX, alpha, beta, gamma, q, B_nec, B_nec);
            }
          else
#endif
            {
              B_nec[0] = mptr->Bx_nec[j];
              B_nec[1] = mptr->By_nec[j];
              B_nec[2] = mptr->Bz_nec[j];
            }

          /* compute total modeled field */
          for (k = 0; k < 3; ++k)
            B_model[k] = B_int[k] + B_ext[k];

          B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

          /* compute vector residuals in NEC frame */
          for (k = 0; k < 3; ++k)
            res[k] = B_nec[k] - B_model[k];

          /* compute scalar residual */
          res[3] = mptr->F[j] - B_model[3];

          gsl_histogram_increment(hf, res[3]);
          gsl_histogram_increment(hz, res[2]);

          if (fabs(res[3]) < b)
            gsl_rstat_add(res[3], rm_f);

          fprintf(fp_res, "%.12f %f %f %f %f %f %.5e %.5e %.5e %.5e %.4e %.4e %.4e %d %d %d %d %d\n",
                  t,
                  lt,
                  mptr->r[j] - 6371.2,
                  wrap180(phi * 180.0 / M_PI),
                  90.0 - theta * 180.0 / M_PI,
                  mptr->qdlat[j],
                  res[3],
                  res[0],
                  res[1],
                  res[2],
                  B_nec[0],
                  B_nec[1],
                  B_nec[2],
                  mptr->satdir[j],
                  ((mptr->flags[j] & MAGDATA_FLG_FIT_MF) && (mptr->flags[j] & MAGDATA_FLG_F)) ? 1 : 0,
                  ((mptr->flags[j] & MAGDATA_FLG_FIT_MF) && (mptr->flags[j] & MAGDATA_FLG_Z)) ? 1 : 0,
                  mptr->flags[j] & MAGDATA_FLG_FIT_EULER ? 1 : 0,
                  mptr->flags[j] & MAGDATA_FLG_DZ_NS ? 1 : 0);
        }

      fprintf(stderr, "done\n");

      fprintf(stderr, "\n=== Scalar Residual Histogram statistics ===\n");
      fprintf(stderr, "\t mean of residuals = %.2f (%.2f) [nT]\n",
              gsl_histogram_mean(hf),
              gsl_rstat_mean(rm_f));
      fprintf(stderr, "\t stddev of residuals = %.2f (%.2f) [nT]\n",
              gsl_histogram_sigma(hf),
              gsl_rstat_sd(rm_f));
      fprintf(stderr, "\t scaling by sum of all bins = %.0f\n",
              gsl_histogram_sum(hf));
      gsl_histogram_scale(hf, 1.0 / gsl_histogram_sum(hf));

      fprintf(stderr, "\n=== B_z Residual Histogram statistics ===\n");
      fprintf(stderr, "\t mean of residuals = %.2f [nT]\n",
              gsl_histogram_mean(hz));
      fprintf(stderr, "\t stddev of residuals = %.2f [nT]\n",
              gsl_histogram_sigma(hz));
      fprintf(stderr, "\t scaling by sum of all bins = %.0f\n",
              gsl_histogram_sum(hz));
      gsl_histogram_scale(hz, 1.0 / gsl_histogram_sum(hz));

      fprintf(stderr, "Writing histogram to %s...",
              filehist);
      for (k = 0; k < nbins; ++k)
        {
          double lower, upper;
          gsl_histogram_get_range(hf, k, &lower, &upper);
          fprintf(fp_hist, "%g %g %f %f\n",
                  lower, upper,
                  gsl_histogram_get(hf, k),
                  gsl_histogram_get(hz, k));
        }
      fprintf(stderr, "done\n");

      gsl_histogram_free(hf);
      gsl_histogram_free(hz);
      gsl_rstat_free(rm_f);

      fclose(fp_res);
      fclose(fp_hist);
    }

  return GSL_SUCCESS;
} /* print_residuals() */

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options] sat1.dat sat2.dat ...\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --maxit | -n num_iterations     - number of robust iterations\n");
  fprintf(stderr, "\t --output_file | -o file         - coefficient output file\n");
  fprintf(stderr, "\t --epoch | -e epoch              - model epoch in decimal years\n");
  fprintf(stderr, "\t --lambda_sv | -v lambda_sv      - SV damping parameter\n");
  fprintf(stderr, "\t --lambda_sa | -a lambda_sa      - SA damping parameter\n");
  fprintf(stderr, "\t --euler | -p period             - Euler bin size in days\n");
  fprintf(stderr, "\t --residual_file | -r file       - residual output file\n");
  fprintf(stderr, "\t --lcurve_file | -l file         - L-curve data file\n");
  fprintf(stderr, "\t --tmin | -b min_time            - minimum data period time in decimal years\n");
  fprintf(stderr, "\t --tmax | -c max_time            - maximum data period time in decimal years\n");
} /* print_help() */

int
main(int argc, char *argv[])
{
  double epoch = MFIELD_EPOCH;
  double R = MFIELD_RE_KM;
  char *outfile = NULL;
  char *resfile = NULL;
  char *Lfile = NULL;
  char *datamap_file = "datamap.dat";
  mfield_workspace *mfield_workspace_p;
  mfield_parameters mfield_params;
  mfield_data_workspace *mfield_data_p;
  gsl_vector *coeffs; /* model coefficients */
  size_t iter = 0;
  size_t maxit = 1;
  double lambda_sv = 0.0;     /* coefficient damping */
  double lambda_sa = 0.0;
  double euler_period = 30.0; /* set to 0 for single set of angles */
  double tmin = -1.0;         /* minimum time for data in years */
  double tmax = -1.0;         /* maximum time for data in years */
  int nsat = 0;               /* number of satellites */
  struct timeval tv0, tv1;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "residual_file", required_argument, NULL, 'r' },
          { "output_file", required_argument, NULL, 'o' },
          { "epoch", required_argument, NULL, 'e' },
          { "lambda_sv", required_argument, NULL, 'v' },
          { "lambda_sa", required_argument, NULL, 'a' },
          { "lcurve_file", required_argument, NULL, 'l' },
          { "maxit", required_argument, NULL, 'n' },
          { "tmin", required_argument, NULL, 'b' },
          { "tmax", required_argument, NULL, 'c' },
          { "euler", required_argument, NULL, 'p' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:c:e:l:n:o:p:r:v:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'b':
            tmin = atof(optarg);
            break;

          case 'c':
            tmax = atof(optarg);
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'n':
            maxit = (size_t) atoi(optarg);
            break;

          case 'e':
            epoch = atof(optarg);
            break;

          case 'v':
            lambda_sv = atof(optarg);
            break;

          case 'a':
            lambda_sa = atof(optarg);
            break;

          case 'p':
            euler_period = atof(optarg);
            break;

          case 'r':
            resfile = optarg;
            break;

          case 'l':
            Lfile = optarg;
            break;

          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  nsat = argc - optind;
  if (nsat == 0)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: epoch = %.2f\n", epoch);
  fprintf(stderr, "main: radius = %g [km]\n", R);
  fprintf(stderr, "main: MF nmax = %d\n", NMAX_MF);
#if MFIELD_FIT_SECVAR
  fprintf(stderr, "main: SV nmax = %d\n", NMAX_SV);
  fprintf(stderr, "main: SV damping = %g\n", lambda_sv);
#endif
#if MFIELD_FIT_SECACC
  fprintf(stderr, "main: SA nmax = %d\n", NMAX_SA);
  fprintf(stderr, "main: SA damping = %g\n", lambda_sa);
#endif
  fprintf(stderr, "main: euler period = %g [days]\n", euler_period);
  fprintf(stderr, "main: tmin = %g\n", tmin);
  fprintf(stderr, "main: tmax = %g\n", tmax);
  fprintf(stderr, "main: number of robust iterations = %zu\n", maxit);
  fprintf(stderr, "main: number of satellites = %d\n", nsat);
  if (outfile)
    fprintf(stderr, "main: output coefficient file = %s\n", outfile);
  if (resfile)
    fprintf(stderr, "main: residual output file = %s\n", resfile);
  if (Lfile)
    fprintf(stderr, "main: L-curve output file = %s\n", Lfile);

  /* allocate data workspace */
  mfield_data_p = mfield_data_alloc(nsat, epoch);

  {
    int satnum = 0;

    while (optind < argc)
      {
        magdata **mdata = &(mfield_data_p->mdata[satnum]);

        assert(satnum++ < nsat);

        fprintf(stderr, "main: reading %s...", argv[optind]);
        gettimeofday(&tv0, NULL);
        *mdata = magdata_read(argv[optind], NULL);
        gettimeofday(&tv1, NULL);

        if (!(*mdata))
          exit(1);

        fprintf(stderr, "done (%zu data total, %g seconds)\n",
                (*mdata)->n, time_diff(tv0, tv1));

        magdata_init(*mdata);
        magdata_calc(*mdata);

        ++optind;
      }
  }

  {
    size_t nflag;

    /* flag any datapoints outside of [tmin,tmax] */
    fprintf(stderr, "main: flagging points outside of time [%g,%g]...", tmin, tmax);
    nflag = mfield_data_filter_time(tmin, tmax, mfield_data_p);
    fprintf(stderr, "done (%zu data flagged)\n", nflag);

#if !MFIELD_FIT_EULER
    fprintf(stderr, "main: flagging Euler-only data points...");
    nflag = mfield_data_filter_euler(mfield_data_p);
    fprintf(stderr, "done (%zu data flagged)\n", nflag);
#endif
  }

  fprintf(stderr, "main: data epoch = %.2f\n", mfield_data_epoch(mfield_data_p));

  /* print spatial coverage maps for each satellite */
  mfield_data_map(datamap_file, mfield_data_p);

  mfield_params.epoch = epoch;
  mfield_params.R = R;
  mfield_params.nmax_mf = NMAX_MF;
  mfield_params.nmax_sv = NMAX_SV;
  mfield_params.nmax_sa = NMAX_SA;
  mfield_params.nsat = nsat;
  mfield_params.euler_period = euler_period;
  mfield_params.mfield_data_p = mfield_data_p;

  /* allocate mfield workspace */
  mfield_workspace_p = mfield_alloc(&mfield_params);

#if MFIELD_SYNTH_DATA
  {
    fprintf(stderr, "main: replacing with synthetic data...");
    gettimeofday(&tv0, NULL);
    replace_synthetic_data(mfield_workspace_p);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
  }
#endif

  /* initialize model parameters */
  mfield_init(mfield_workspace_p);

  /* coefficient damping parameters */
  mfield_set_damping(lambda_sv, lambda_sa, mfield_workspace_p);

  /* construct initial guess vector from IGRF */
  coeffs = gsl_vector_alloc(mfield_workspace_p->p);
  fprintf(stderr, "main: constructing initial coefficient vector...");
  initial_guess(coeffs, mfield_workspace_p);
  fprintf(stderr, "done\n");

  while (iter++ < maxit)
    {
      fprintf(stderr, "main: ROBUST ITERATION %zu/%zu\n", iter, maxit);

      mfield_calc_nonlinear(coeffs, mfield_workspace_p);

      /* reset workspace for a new iteration */
      mfield_reset(mfield_workspace_p);
    }

#if MFIELD_CALC_UNCERTAINTIES
  /* calculate errors in coefficients */
  {
    struct timeval tv0, tv1;

    fprintf(stderr, "main: calculating coefficient uncertainties...");
    gettimeofday(&tv0, NULL);
    mfield_calc_uncertainties(mfield_workspace_p);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
  }
#endif

  /* L-curve data */
  if (Lfile)
    {
      gsl_vector *res_f = gsl_multifit_fdfridge_residual(mfield_workspace_p->fdf_s);
      double xnorm = gsl_blas_dnrm2(mfield_workspace_p->c);
      double fnorm = gsl_blas_dnrm2(res_f);
      FILE *fp = fopen(Lfile, "a");

      if (!fp)
        {
          fprintf(stderr, "main: unable to open %s: %s\n",
                  Lfile, strerror(errno));
        }
      else
        {
          fprintf(fp, "%.12e %.12e %f %f\n",
                  log(fnorm),
                  log(xnorm),
                  lambda_sv,
                  lambda_sa);

          fclose(fp);
        }
    }

  if (outfile)
    {
      fprintf(stderr, "main: writing coefficients to %s...", outfile);
      mfield_write(outfile, mfield_workspace_p);
      fprintf(stderr, "done\n");
    }

  if (resfile)
    print_residuals(resfile, mfield_workspace_p);

  print_spectrum("mfield.s", mfield_workspace_p);

  /* print coefficients */
  {
    size_t n;
    int m;

#if MFIELD_FIT_EXTFIELD
    char *ext_file = "coeffs.ext";
    FILE *fp = fopen(ext_file, "w");
    gsl_vector_view ev = gsl_vector_subvector(coeffs, mfield_workspace_p->ext_offset,
                                              mfield_workspace_p->next);

    fprintf(stderr, "main: printing external coefficients to %s...", ext_file);
    gsl_vector_fprintf(fp, &ev.vector, "%g");
    fprintf(stderr, "done\n");

    fclose(fp);
#endif

#if MFIELD_FIT_EULER
    /* print Euler angles */
    for (n = 0; n < mfield_workspace_p->nsat; ++n)
      {
        magdata *mptr = mfield_data_ptr(n, mfield_workspace_p->data_workspace_p);

        if (mptr->global_flags & MAGDATA_GLOBFLG_EULER)
          {
            double t0 = mfield_workspace_p->data_workspace_p->t0[n];
            size_t euler_idx = mfield_euler_idx(n, t0, mfield_workspace_p);
            double alpha = gsl_vector_get(coeffs, euler_idx);
            double beta = gsl_vector_get(coeffs, euler_idx + 1);
            double gamma = gsl_vector_get(coeffs, euler_idx + 2);
            char filename[2048];

            fprintf(stderr, "main: satellite %zu: alpha = %f beta = %f gamma = %f [deg]\n",
                    n,
                    wrap180(alpha * 180.0 / M_PI),
                    wrap180(beta * 180.0 / M_PI),
                    wrap180(gamma * 180.0 / M_PI));

            sprintf(filename, "euler.%zu", n);
            fprintf(stderr, "main: satellite %zu: printing Euler angles to %s...", n, filename);
            mfield_euler_print(filename, n, mfield_workspace_p);
            fprintf(stderr, "done\n");
          }
      }
#endif

    fprintf(stderr, "main: printing internal coefficients up to degree 3\n");
    for (n = 1; n <= GSL_MIN(3, NMAX_MF); ++n)
      {
        int ni = (int) n;
        for (m = -ni; m <= ni; ++m)
          {
            int mabs = abs(m);
            size_t cidx = mfield_coeff_nmidx(n, m);
            char c = (m < 0) ? 'h' : 'g';
            double gnm = mfield_get_mf(coeffs, cidx, mfield_workspace_p);
            double dgnm = mfield_get_sv(coeffs, cidx, mfield_workspace_p);
            double ddgnm = mfield_get_sa(coeffs, cidx, mfield_workspace_p);

            fprintf(stderr, "%c(%d,%d) = %12g (%12g,%12g)\n", c, ni, mabs,
                    gnm, dgnm, ddgnm);
          }
      }
  }

  mfield_free(mfield_workspace_p);
  mfield_data_free(mfield_data_p);
  gsl_vector_free(coeffs);

  return 0;
} /* main() */
