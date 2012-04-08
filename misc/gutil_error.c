/******************************************************************************%
**
**    Copyright (C) 1998-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"



/*******************************************************************************
 *
 ******************************************************************************/
void eprintf(const char *fmt, ...) {

     va_list ap;

     va_start(ap, fmt);

     vfprintf(stderr, fmt, ap);

     va_end(ap);
}



/*******************************************************************************
 *
 ******************************************************************************/
void check_arg_count(int i, int argc, int n, const char *s) {

     if (i + n >= argc) {
          fprintf(stderr, "ERROR: not enough arguments for %s\n", s);
          exit(1);
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
int strtoi_errmsg(const char *string, const char *name, int *value) {

     char *endptr;

     errno = 0;

     if (string == NULL) {
          fprintf(stderr, "ERROR: invalid value for %s: NULL\n", name);
          return -1;
     }

     if (string == '\0') {
          fprintf(stderr, "ERROR: invalid value for %s: empty string\n", name);
          return -1;
     }

     *value = strtol(string, &endptr, 0);

     if (errno == ERANGE) {
          if (*value == INT_MIN) {
               fprintf(stderr, "ERROR: invalid value for %s, too small: %s\n", name, string);
               return -1;
          }
          else {
               fprintf(stderr, "ERROR: invalid value for %s, too large: %s\n", name, string);
               return -1;
          }
     }

     if (*endptr != '\0') {
          fprintf(stderr, "ERROR: invalid value for %s: %s\n", name, string);
          return -1;
     }
/*
     if (sscanf(string, "%d", value) != 1) {
          fprintf(stderr, "ERROR: invalid value for %s: %s\n", name, string);
          return -1;
     }
*/
     return 0;
}



int strtoi_errmsg_exit(const char *string, const char *name) {

     int value;

     if (strtoi_errmsg(string, name, &value))
          exit(1);

     return value;
}



long strtol_errmsg(const char *string, const char *name, long *value) {

     char *endptr;

     errno = 0;

     if (string == NULL) {
          fprintf(stderr, "ERROR: invalid value for %s: NULL\n", name);
          exit(1);
     }

     if (string == '\0') {
          fprintf(stderr, "ERROR: invalid value for %s: empty string\n", name);
          exit(1);
     }

     *value = strtol(string, &endptr, 0);

     if (errno == ERANGE) {
          if (*value == LONG_MIN) {
               fprintf(stderr, "ERROR: invalid value for %s, too small: %s\n", name, string);
               exit(1);
          }
          else {
               fprintf(stderr, "ERROR: invalid value for %s, too large: %s\n", name, string);
               exit(1);
          }
     }

     if (*endptr != '\0') {
          fprintf(stderr, "ERROR: invalid value for %s: %s\n", name, string);
          exit(1);
     }
/*
     if (sscanf(string, "%ld", &value) != 1) {
          fprintf(stderr, "ERROR: invalid value for %s: %s\n", name, string);
          exit(1);
     }
*/
     return 0;
}



long strtol_errmsg_exit(const char *string, const char *name) {

     long value;

     if (strtol_errmsg(string, name, &value))
          exit(1);

     return value;
}



double strtod_errmsg(const char *string, const char *name, double *value) {

     char *endptr;

     errno = 0;

     if (string == NULL) {
          fprintf(stderr, "ERROR: invalid value for %s: NULL\n", name);
          exit(1);
     }

     if (string == '\0') {
          fprintf(stderr, "ERROR: invalid value for %s: empty string\n", name);
          exit(1);
     }

     *value = strtod(string, &endptr);

     if (errno == ERANGE) {
          if (*value == HUGE_VAL) {
               fprintf(stderr, "ERROR: invalid value for %s, too small: %s\n", name, string);
               exit(1);
          }
          else {
               fprintf(stderr, "ERROR: invalid value for %s, too large: %s\n", name, string);
               exit(1);
          }
     }

     if (*endptr != '\0') {
          fprintf(stderr, "ERROR: invalid value for %s: %s\n", name, string);
          exit(1);
     }
/*
     if (sscanf(string, "%lf", &value) != 1) {
          fprintf(stderr, "ERROR: invalid value for %s: %s\n", name, string);
          exit(1);
     }
*/
     return 0;
}



double strtod_errmsg_exit(const char *string, const char *name) {

     double value;

     if (strtod_errmsg(string, name, &value))
          exit(1);

     return value;
}



/*******************************************************************************
 *
 ******************************************************************************/
int system_errmsg(const char *cmd, int flag) {

     int r;

     r = system(cmd);

     if (r < 0) {
          fprintf(stderr, "ERROR: system() function failed\n");
          return -1;
     }
     else
     if (r > 0 && flag) {
          fprintf(stderr, "ERROR: system(), command failed: %s\n", cmd);
          return -1;
     }

     return 0;
}



void system_errmsg_exit(const char *cmd, int flag) {

     int r;

     r = system(cmd);

     if (r < 0) {
          fprintf(stderr, "ERROR: system() function failed\n");
          exit(1);
     }
     else
     if (r > 0 && flag) {
          fprintf(stderr, "ERROR: system(), command failed: %s\n", cmd);
          exit(1);
     }
}
