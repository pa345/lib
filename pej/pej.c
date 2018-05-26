/*
 * pej.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pej.h"

pej_workspace *
pej_alloc()
{
  pej_workspace *w;

  w = calloc(1, sizeof(pej_workspace));
  if (!w)
    return NULL;

  return w;
}

void
pej_free(pej_workspace *w)
{
  free(w);
}
