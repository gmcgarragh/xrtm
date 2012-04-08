/******************************************************************************%
**
**    Copyright (C) 1998-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


void read_newline(FILE *fp) {

     while (fgetc(fp) != '\n') ;
}


int read_comment(FILE *fp, char comment) {

     int c;

     c = fgetc(fp);

     if (c == EOF)
          return 0;

     if (c == '\n')
          return 1;

     if (c != comment) {
          ungetc(c, fp);
          return 0;
     }

     do {
          c = fgetc(fp);
     } while (c != EOF && c != '\n');

     return 1;
}



void read_comments(FILE *fp, char comment) {

     while (read_comment(fp, comment)) ;
}
