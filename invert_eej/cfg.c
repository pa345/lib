/*
 * cfg.c
 *
 * This module contains routines to read the configuration file and
 * store the values of the parameters
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <libconfig.h>

#include "cfg.h"

static int cfg_read(const char *config_filename, cfg_workspace *w);

cfg_parameters cfg_params;

/* configuration directives and their locations in cfg_params */
static cfg_struct cfg_dirs[] = {
  { "input_file_prev", &(cfg_params.prev_input_filename), CFG_STRING|CFG_OPTIONAL },
  { "input_file", &(cfg_params.input_filename), CFG_STRING|CFG_OPTIONAL },
  { "output_file", &(cfg_params.output_filename), CFG_STRING|CFG_OPTIONAL },
  { "log_dir", &(cfg_params.log_dir), CFG_STRING|CFG_OPTIONAL },
  { "test_output_dir", &(cfg_params.test_output_dir), CFG_STRING|CFG_OPTIONAL },

  { "kp_data_file", &(cfg_params.kp_data_file), CFG_STRING|CFG_OPTIONAL },
  { "f107_data_file", &(cfg_params.f107_data_file), CFG_STRING|CFG_OPTIONAL },
  { "dst_data_file", &(cfg_params.dst_data_file), CFG_STRING|CFG_OPTIONAL },
  { "igrf_cof_file", &(cfg_params.igrf_cof_file), CFG_STRING|CFG_OPTIONAL },
  { "lith_cof_file", &(cfg_params.lith_cof_file), CFG_STRING|CFG_OPTIONAL },
  { "core_cof_file", &(cfg_params.core_cof_file), CFG_STRING|CFG_OPTIONAL },

  { "min_local_time", &(cfg_params.lt_min), CFG_DOUBLE|CFG_OPTIONAL },
  { "max_local_time", &(cfg_params.lt_max), CFG_DOUBLE|CFG_OPTIONAL },
  { "max_latitude_gap", &(cfg_params.max_lat_gap), CFG_DOUBLE|CFG_OPTIONAL },
  { "max_time_gap", &(cfg_params.max_time_gap), CFG_DOUBLE|CFG_OPTIONAL },
  { "dst_num_days", &(cfg_params.dst_num_days), CFG_INT|CFG_OPTIONAL },
  { "current_shell_altitude", &(cfg_params.curr_altitude), CFG_DOUBLE|CFG_OPTIONAL },
  { "ncurr", &(cfg_params.ncurr), CFG_INT|CFG_OPTIONAL },

  { "prof_lat_min", &(cfg_params.prof_lat_min), CFG_DOUBLE|CFG_OPTIONAL },
  { "prof_lat_max", &(cfg_params.prof_lat_max), CFG_DOUBLE|CFG_OPTIONAL },
  { "prof_lat_step", &(cfg_params.prof_lat_step), CFG_DOUBLE|CFG_OPTIONAL },

  { "nr", &(cfg_params.nr), CFG_INT|CFG_OPTIONAL },
  { "ntheta", &(cfg_params.ntheta), CFG_INT|CFG_OPTIONAL },

  { "r_earth", &(cfg_params.r_earth), CFG_DOUBLE|CFG_OPTIONAL },

  { "alt_min", &(cfg_params.alt_min), CFG_DOUBLE|CFG_OPTIONAL },
  { "alt_max", &(cfg_params.alt_max), CFG_DOUBLE|CFG_OPTIONAL },
  { "theta_min", &(cfg_params.theta_min), CFG_DOUBLE|CFG_OPTIONAL },
  { "theta_max", &(cfg_params.theta_max), CFG_DOUBLE|CFG_OPTIONAL },

  { "kp_max", &(cfg_params.kp_max), CFG_DOUBLE|CFG_OPTIONAL },
  { "min_longitude", &(cfg_params.lon_min), CFG_DOUBLE|CFG_OPTIONAL },
  { "max_longitude", &(cfg_params.lon_max), CFG_DOUBLE|CFG_OPTIONAL },
  { "track_qdmax", &(cfg_params.track_qdmax), CFG_DOUBLE|CFG_OPTIONAL },

  { "sq_nmax_int", &(cfg_params.sq_nmax_int), CFG_INT|CFG_OPTIONAL },
  { "sq_mmax_int", &(cfg_params.sq_mmax_int), CFG_INT|CFG_OPTIONAL },
  { "sq_nmax_ext", &(cfg_params.sq_nmax_ext), CFG_INT|CFG_OPTIONAL },
  { "sq_mmax_ext", &(cfg_params.sq_mmax_ext), CFG_INT|CFG_OPTIONAL },

  { "sq_qdmin", &(cfg_params.sq_qdmin), CFG_DOUBLE|CFG_OPTIONAL },
  { "sq_qdmax", &(cfg_params.sq_qdmax), CFG_DOUBLE|CFG_OPTIONAL },

  { "calc_field_models", &(cfg_params.calc_field_models), CFG_INT|CFG_OPTIONAL },
  { "main_nmax_int", &(cfg_params.main_nmax_int), CFG_INT|CFG_OPTIONAL },
  { "crust_nmax_int", &(cfg_params.crust_nmax_int), CFG_INT|CFG_OPTIONAL },

  { 0, 0, 0 }
};

cfg_workspace *
cfg_alloc(const char *config_filename)
{
  cfg_workspace *w;
  int s;
  cfg_struct *cptr;

  w = calloc(1, sizeof(cfg_workspace));
  if (!w)
    {
      fprintf(stderr, "cfg_alloc: calloc failed: %s\n", strerror(errno));
      return 0;
    }

  /* initialize cfg_params structure */
  for (cptr = cfg_dirs; cptr->name; cptr++)
    {
      if (cptr->flags & CFG_STRING)
        *(char **)cptr->location = NULL;
      else if (cptr->flags & CFG_INT)
        *(int *)cptr->location = -9999;
      else if (cptr->flags & CFG_DOUBLE)
        *(double *)cptr->location = -9999.0;
    }

  s = cfg_read(config_filename, w);
  if (s)
    {
      fprintf(stderr, "cfg_alloc: cfg_read failed\n");
      cfg_free(w);
      return 0;
    }

  return w;
} /* cfg_alloc() */

void
cfg_free(cfg_workspace *w)
{
  config_destroy(&(w->config));

  free(w);
} /* cfg_free() */

/*
cfg_read()
  Read configuration file and store parameters into w->params
*/

static int
cfg_read(const char *config_filename, cfg_workspace *w)
{
  int s = 0;
  cfg_struct *cptr;

  config_init(&(w->config));

  s = config_read_file(&(w->config), config_filename);
  if (s != CONFIG_TRUE)
    {
       fprintf(stderr, "cfg_read: %s:%d - %s\n",
               config_filename,
               config_error_line(&(w->config)),
               config_error_text(&(w->config)));

      return 1;
    }

  for (cptr = cfg_dirs; cptr->name; cptr++)
    {
      if (cptr->flags & CFG_STRING)
        {
          s = config_lookup_string(&(w->config),
                                   cptr->name,
                                   (const char **) cptr->location);
        }
      else if (cptr->flags & CFG_INT)
        {
          s = config_lookup_int(&(w->config),
                                cptr->name,
                                (int *) cptr->location);
        }
      else
        {
          s = config_lookup_float(&(w->config),
                                  cptr->name,
                                  (double *) cptr->location);
        }

      if (s != CONFIG_TRUE && !(cptr->flags & CFG_OPTIONAL))
        {
          fprintf(stderr, "cfg_read: error: %s not specified in config file %s\n",
                  cptr->name,
                  config_filename);
          return 1;
        }
    }

  /* save a copy into w->params */
  w->params = cfg_params;

  return 0;
} /* cfg_read() */
