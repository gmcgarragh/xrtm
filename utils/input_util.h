/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef INPUT_UTIL_H
#define INPUT_UTIL_H

#include <xrtm_fd.h>
#include <xrtm_model.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
     int n_layers;
     int n_derivs;
     int fd_method;
     char **coef_files;
     char ***coef_files_l;
} misc_data;


void misc_init(xrtm_data *d, misc_data *md, int fd_method);
void misc_free(misc_data *md);

int xrtm_sread_input(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, char *s, int use_aligned, int use_array2, int use_dash, int use_equal);
int xrtm_fread_input(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, const char *filename, int use_aligned, int use_array2, int use_dash, int use_equal);
int xrtm_swrite_input(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, char *s, char nl, int use_aligned, int use_array2, int use_dash, int use_equal);
int xrtm_fwrite_input_fn(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, const char *filename, char nl, int use_aligned, int use_array2, int use_dash, int use_equal);
int xrtm_fwrite_input_fp(xrtm_data *d, xrtm_fd_data *fd, misc_data *md, FILE *fp, char nl, int use_aligned, int use_array2, int use_dash, int use_equal);


#ifdef __cplusplus
}
#endif

#endif /* INPUT_UTIL_H */
