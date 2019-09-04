/*
 * mag_log_grad.c
 *
 * Contains routines which log various steps of the inversion
 * for E/W gradients
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_interp.h>

#include "mag.h"
#include "pde.h"

#include "pde_common.c"

/*
mag_log_B2_grad()
  Log computed B^(2) residuals to a file

Inputs: header - 1 = print header, 0 = don't
        w      - workspace

Notes
1) track is stored in w->track
*/

int
mag_log_B2_grad(const int header, const mag_workspace *w)
{
  int s = 0;
  const size_t downsample = 5; /* downsample to keep file size reasonable */
  size_t i;
  mag_track *track = (mag_track *) &(w->track);

  if (header)
    {
      /* print header information */
      i = 1;
      log_proc(w->log_B2_grad, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: timestamp at gradient point (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: QD latitude (degrees)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: QD latitude at gradient point (degrees)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: X^(1) (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: Y^(1) (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: Z^(1) (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: F^(1) (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: X^(1) at gradient point (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: Y^(1) at gradient point (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: Z^(1) at gradient point (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: F^(1) at gradient point (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: X^(2) (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: Y^(2) (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: Z^(2) (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: F^(2) (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: internal Sq model X (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: internal Sq model Y (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: internal Sq model Z (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: external Sq model X (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: external Sq model Y (nT)\n", i++);
      log_proc(w->log_B2_grad, "# Field %zu: external Sq model Z (nT)\n", i++);
      return s;
    }

  for (i = 0; i < track->n; i += downsample)
    {
      time_t unix_time = satdata_epoch2timet(track->t[i]);
      time_t unix_time2 = satdata_epoch2timet(track->t_grad[i]);

      log_proc(w->log_B2_grad, "%ld %ld %8.4f %8.4f %8.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
               unix_time,
               unix_time2,
               track->qdlat[i],
               track->qdlat_grad[i],
               track->X1[i],
               track->Y1[i],
               track->Z1[i],
               track->F1[i],
               track->X1_grad[i],
               track->Y1_grad[i],
               track->Z1_grad[i],
               track->F1_grad[i],
               track->X2[i],
               track->Y2[i],
               track->Z2[i],
               track->F2[i],
               track->X_Sq_int[i],
               track->Y_Sq_int[i],
               track->Z_Sq_int[i],
               track->X_Sq_ext[i],
               track->Y_Sq_ext[i],
               track->Z_Sq_ext[i]);
    }

  log_proc(w->log_B2_grad, "\n\n");

  return s;
}
