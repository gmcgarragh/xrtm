/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


int gpass_newline(FILE *fp) {

     int c;

     while ((c = fgetc(fp)) != '\n')  {
          if (c == EOF)
               return 1;
     }

     return 0;
}



int gcount_lines(FILE *fp) {

     int i;

     long pos;

     pos = ftell(fp);

     for (i = 0; gpass_newline(fp) == 0; ++i) ;

     if (fseek(fp, pos, SEEK_SET)) {
          fprintf(stderr, "ERROR: Problem setting file position");
          return -1;
     }

     return i;
}



int gparse_comment(FILE *fp, char comment) {

     int c;

     c = fgetc(fp);

     if (c == EOF)
          return 1;

     if (c == '\n')
          return 0;

     if (c != comment) {
          ungetc(c, fp);
          return 1;
     }

     do {
          c = fgetc(fp);
     } while (c != EOF && c != '\n');

     return 0;
}



void gparse_comments(FILE *fp, char comment) {

     while (gparse_comment(fp, comment) == 0) ;
}



void trimm_string(char *s) {

     int i;
     int j;
     int k;

     int n;

     for (i = 0; s[i] == ' ' || s[i] == '\t'; ++i) ;

     n = strlen(s + i);

     if (n == 0)
          j = 0;
     else {
          for (j = n - 1; s[i + j] == ' ' || s[i + j] == '\t'; --j) ;
     }

     if (i > 0) {
          for (k = 0; k < j + 1; ++k)
               s[k] = s[i + k];
          s[k] = '\0';
     }

     s[j + 1] = '\0';
}



void gparse_init(gparse_data *d, char comment) {

     d->comment       = comment;
     d->count         = 0;
     d->newline_found = 1;
}



char gparse_string(gparse_data *d, FILE *fp, int n, char *s, const char *delims) {

     int i;

     int c;

     if (d->newline_found)
         d->count++;

     if (! d->newline_found)
          c = fgetc(fp);

     while (d->newline_found) {

     c = fgetc(fp);

     while (c != EOF && c != '\n' && c != d->comment && (c == ' ' || c == '\t'))
          c = fgetc(fp);

     if (c == EOF)
          return '\0';
     else
     if (c == '\n') {
          d->count++;
          continue;
     }
     else
     if (c == d->comment) {
          gpass_newline(fp);
          d->count++;
          continue;
     }
     else {
          d->newline_found = 0;
          break;
     }

     }

     while (c != EOF && c != '\n' && c != d->comment && strchr(delims, c) != NULL)
          c = fgetc(fp);

     if (c == EOF)
          return '\0';
     else
     if (c == '\n')
          return '\0';
     else
     if (c == d->comment)
          return '\0';

     i = 0;
     while (c != EOF && c != '\n' && c != d->comment && strchr(delims, c) == NULL && i < n - 1) {
          s[i++] = (char) c;
          c = fgetc(fp);
     }

     if (c == '\n')
          d->newline_found = 1;
     else
     if (c == d->comment) {
          gpass_newline(fp);
          d->newline_found = 1;
     }

     s[i] = '\0';

     return (char) c;
}



int gparse_newline(gparse_data *d, FILE *fp) {

     int c;

     if (d->newline_found)
         return 0;

     c = fgetc(fp);
     while (c != EOF && c != '\n' && c != d->comment && (c == ' ' || c == '\t'))
          c = fgetc(fp);

     if (c == EOF)
          return 1;
     else
     if (c == '\n') {
          d->newline_found = 1;
          return 0;
     }
     else
     if (c == d->comment) {
          gpass_newline(fp);
          d->newline_found = 1;
          return 0;
     }
     else
          return 1;
}



char gparse_char(gparse_data *d, FILE *fp, const char *delims) {

     char a;

     char temp[LN];

     if (gparse_string(d, fp, LN, temp, delims) == '\0')
          return CHAR_MAX;

     if (sscanf(temp, "%c", &a) < 1)
          return CHAR_MAX;

     return a;
}



char *gparse_quoted(gparse_data *d, FILE *fp, int n, char *s, const char *delims) {

     return NULL;
}



int gparse_int(gparse_data *d, FILE *fp, const char *delims) {

     int a;

     char temp[LN];

     if (gparse_string(d, fp, LN, temp, delims) == '\0')
          return INT_MAX;

     if (sscanf(temp, "%d", &a) < 1)
          return INT_MAX;

     return a;
}



double gparse_double(gparse_data *d, FILE *fp, const char *delims) {

     double a;

     char temp[LN];

     if (gparse_string(d, fp, LN, temp, delims) == '\0')
          return DBL_MAX;

     if (sscanf(temp, "%lf", &a) < 1)
          return DBL_MAX;

     return a;
}
