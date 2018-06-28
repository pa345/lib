#include <math.h>
#include <getopt.h>

/*
pde_theta()
  Compute angle theta corresponding to an index idx \in [0, NTHETA - 1]
in radians
*/
static inline double
pde_theta(size_t idx, pde_workspace *w)
{
  return (w->theta_min + idx * w->dtheta);
}

/* sin(theta) */
static inline double
pde_sint(size_t idx, pde_workspace *w)
{
  double theta = pde_theta(idx, w);

  return sin(theta);
}

/*
pde_r_m()
  Compute radius r corresponding to an index idx \in [0, NR - 1] in
meters
*/

static inline double
pde_r_m(size_t idx, pde_workspace *w)
{
  return (w->rmin + idx * w->dr);
} /* pde_r_m() */

/*
pde_r()
  Compute radius r corresponding to an index idx \in [0, NR - 1] in
dimensionless units
*/

static inline double
pde_r(size_t idx, pde_workspace *w)
{
  return (pde_r_m(idx, w) / w->r_s);
} /* pde_r() */

/*
pde_r_km()
  Compute radius r corresponding to an index idx \in [0, NR - 1] in km
*/
static inline double
pde_r_km(size_t idx, pde_workspace *w)
{
  return pde_r_m(idx, w) / 1000.0;
} /* pde_r_km() */

static inline double
pde_dr(pde_workspace *w)
{
  return w->dr / w->r_s;
}

static inline double
pde_dr_sq(pde_workspace *w)
{
  return (w->dr * w->dr / w->r_s / w->r_s);
}

static inline double
pde_dtheta(pde_workspace *w)
{
  return (w->dtheta);
}

static inline double
pde_dtheta_sq(pde_workspace *w)
{
  return (w->dtheta * w->dtheta);
}

/*
pde_ridx()
  Determine index r in [0, NR - 1] corresponding to given radius

Inputs: r - altitude (relative to earth surface) in km

*/
static inline size_t
pde_ridx(double r, pde_workspace *w)
{
  size_t idx;

  idx = (size_t) (((r + R_EARTH_KM) * 1.0e3 - w->rmin) / w->dr);

  return idx;
}

/*
pde_thidx()
  Determine index theta in [0, NTHETA - 1] corresponding to given theta

Inputs: theta - colatitude in radians

*/
static inline size_t
pde_thidx(double theta, pde_workspace *w)
{
  size_t idx;

  idx = (size_t) ((theta - w->theta_min) / w->dtheta);

  return idx;
}

inline static double
E_phi(size_t i, size_t j, pde_workspace *w)
{
  double r = pde_r_km(i, w);
  double theta = pde_theta(j, w);
  double Ep;

  Ep = R_BOTTOM * w->E_phi0 / r / sin(theta) / w->E_s;

  return Ep;
}
