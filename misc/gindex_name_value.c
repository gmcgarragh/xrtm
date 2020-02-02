/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"

#include "gindex_name_value.h"


/*******************************************************************************
 *
 ******************************************************************************/
int max_name_length(const char **names, int n) {

     int i;

     int length;
     int length_max = 0;

     for (i = 0; i < n; ++i) {
          length = strlen(names[i]);
          if (length > length_max)
               length_max = length;
     }

     return length_max;
}



/*******************************************************************************
 *
 ******************************************************************************/
const char *index_to_name(int index, const char **names, int n, const char *desc) {

     if (index < 0 || index >= n) {
          fprintf(stderr, "ERROR: invalid %s index: %d\n", desc, index);
          return NULL;
     }

     return names[index];
}



int name_to_index(const char *name, const char **names, int n, const char *desc) {

     int i;

     for (i = 0; i < n; ++i) {
          if (strcmp(names[i], name) == 0) {
               return i;
          }
     }

     fprintf(stderr, "ERROR: invalid %s name: %s\n", desc, name);

     return -1;
}



/*******************************************************************************
 *
 ******************************************************************************/
long index_to_value(int index, long *values, int n, const char *desc) {

     if (index < 0 || index >= n) {
          fprintf(stderr, "ERROR: invalid %s index: %d\n", desc, index);
          return -1;
     }

     return values[index];
}



int value_to_index(long value, long *values, int n, const char *desc) {

     int i;

     for (i = 0; i < n; ++i) {
          if (value == values[i])
               return i;
     }

     fprintf(stderr, "ERROR: invalid %s value: %ld\n", desc, value);

     return -1;
}



int mask_to_index_list(long mask, long *values, int n_values, int *indexes, int length) {

     int i;
     int ii;

     ii = 0;

     for (i = 0; i < n_values; ++i) {
          if (mask & values[i])
               indexes[ii++] = i;

          if (ii >= length)
               return ii;
     }

     return ii;
}



/*******************************************************************************
 *
 ******************************************************************************/
long name_to_value(const char *name, long *values, const char **names, int n, const char *desc) {

     int i;

     for (i = 0; i < n; ++i) {
          if (strcmp(name, names[i]) == 0)
               return values[i];
     }

     fprintf(stderr, "ERROR: invalid %s name: %s\n", desc, name);

     return -1;
}



const char *value_to_name(long value, long *values, const char **names, int n, const char *desc) {

     int i;

     for (i = 0; i < n; ++i) {
          if (value == values[i])
               return names[i];
     }

     fprintf(stderr, "ERROR: invalid %s value: %ld\n", desc, value);

     return NULL;
}



char *mask_to_name_list(long mask, long *values, int n_values, const char **names, char *s, int length) {

     int i;
     int n;

     int flag;

     n = 0;

     flag = 0;
     for (i = 0; i < n_values; ++i) {
          if (mask & values[i]) {
               if (flag)
                    n += snprintf(s + n, length - n, ",");
               n += snprintf(s + n, length - n, "%s", names[i]);
               flag = 1;
          }
     }

     s[n] = '\0';

     return s;
}



long name_list_to_mask(char *s, long *values, const char **names, int n, const char *desc) {

     char *token;
     char *lasts;

     int temp;
     int mask;

     mask = 0;

     if ((token = strtok_r(s, ", \t", &lasts))) {
          do {
               if ((temp = name_to_value(token, values, names, n, desc)) < 0) {
                    fprintf(stderr, "ERROR: name_to_value()\n");
                    return -1;
               }
               mask |= temp;
          } while ((token = strtok_r(NULL, ", \t", &lasts)));
     }

     return mask;
}



/*******************************************************************************
 *
 ******************************************************************************/
int mask_to_value_list(long mask, long *values, int n_values, long *values2, int length) {

     int i;
     int ii;

     ii = 0;

     for (i = 0; i < n_values; ++i) {
          if (mask & values[i])
               values2[ii++] = values[i];

          if (ii >= length)
               return ii;
     }

     return ii;
}



long value_list_to_mask(long *values, int n_values) {

     int i;

     int mask;

     mask = 0;
     for (i = 0; i < n_values; ++i)
          mask |= values[i];

     return mask;
}
