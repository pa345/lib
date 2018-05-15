/*
 * attitude.h
 */

#ifndef INCLUDED_attitude_h
#define INCLUDED_attitude_h

#include <satdata/satdata.h>

#include "track.h"

/*
 * Prototypes
 */

int attitude_correct(const char *filename, const satdata_mag * data, const track_workspace * track_p);

#endif /* INCLUDED_attitude_h */
