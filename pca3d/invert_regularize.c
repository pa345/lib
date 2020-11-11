/*
 * invert_regularize.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include "invert.h"
#include "io.h"
#include "pca3d.h"

/*
invert_regularize()
  Build regularization matrix

Inputs: L - (output) diagonal regularization matrix
        w - workspace
*/

int
invert_regularize(gsl_vector * L, invert_workspace * w)
{
  invert_tmode_workspace * tmode_p = w->tmode_workspace_p;
  invert_smode_workspace * smode_p = w->smode_workspace_p;
  size_t nf, ns, nt;
  char buf[1024];

  for (nf = 0; nf < w->nfreq; ++nf)
    {
      gsl_vector *S;

      /* read singular values for this band */
      sprintf(buf, "%s_%zu", PCA3D_STAGE3B_SVAL_DAT, nf + 1);
      S = pca3d_read_vector(buf);

      for (ns = 0; ns < smode_p->nmodes[nf]; ++ns)
        {
          double Sj = gsl_vector_get(S, ns);
          double Sjinv = 1.0 / Sj;

          for (nt = 0; nt < tmode_p->nmodes[nf]; ++nt)
            {
              size_t idx = invert_coeff_idx(nf, nt, ns, w);
              gsl_vector_set(L, idx, Sjinv);
              gsl_vector_set(L, w->p_complex + idx, Sjinv);
            }
        }

      gsl_vector_free(S);
    }

  return 0;
}
