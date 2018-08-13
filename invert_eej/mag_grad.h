/*
 * mag_grad.h
 */

#ifndef INCLUDED_mag_grad_h
#define INCLUDED_mag_grad_h

/*
 * Prototypes
 */

int mag_grad_proc(const mag_params *params, track_workspace *track_p, satdata_mag *data,
                  track_workspace *track_p2, satdata_mag *data2, mag_workspace *w);

#endif /* INCLUDED_mag_grad_h */
