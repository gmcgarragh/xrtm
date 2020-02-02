/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


/*******************************************************************************
 *
 ******************************************************************************/
C_TYPE XCAT(parse_scaler_, FILE_NAME)(input_data *input) {

     int r;

     lex_type_data lex_type;

     r = yylex(&input->locus, &lex_type);

     if (r == LEX_TYPE)
          return lex_type.LEX_MEMBER;
#ifdef LEX_TYPE2
     else
     if (r == LEX_TYPE2)
          return lex_type.LEX_MEMBER2;
#endif
/*
     else
     if (r == LEX_IDENTIFIER)
          ;
*/
     else {
          input_error(&input->locus, "expected %s constant", FILE_STRING);
          exit(1);
     }
}



void *XCAT(FILE_NAME, _array_parse1)(input_data *input, int n, int *dimens, int i) {

     int j;

     int r;

     size_t size;

     void *array;

     lex_type_data lex_type;

     r = yylex(&input->locus, &lex_type);

     if (r != '{') {
          input_error(&input->locus, "expected '{'");
          exit(1);
     }

     if (i == n - 1)
          size = sizeof(C_TYPE);
     else
          size = sizeof(void *);

     array = malloc(dimens[i] * size);

     for (j = 0; j < dimens[i]; ++j) {
          if (i == n - 1)
               ((C_TYPE *) array)[j] = XCAT(parse_scaler_, FILE_NAME)(input);
          else
               ((void **)  array)[j] = XCAT(FILE_NAME, _array_parse1)(input, n, dimens, i + 1);
          if (j < dimens[i] - 1)
               parse_char(&input->locus, ',');
     }

     r = yylex(&input->locus, &lex_type);

     if (r != '}') {
          input_error(&input->locus, "expected '}'");
          exit(1);
     }

     return array;
}



void *XCAT(FILE_NAME, _array_parse2)(input_data *input, int n, int *dimens, int i) {

     int j;

     size_t size;

     void *array;

     if (i == n - 1)
          size = sizeof(C_TYPE);
     else
          size = sizeof(void *);

     array = malloc(dimens[i] * size);

     for (j = 0; j < dimens[i]; ++j) {
          if (i == n - 1) {
               ((C_TYPE *) array)[j] = XCAT(parse_scaler_, FILE_NAME)(input);
               if (j < dimens[i] - 1)
                    parse_char(&input->locus, ',');
          }
          else {
               ((void **)  array)[j] = XCAT(FILE_NAME, _array_parse2)(input, n, dimens, i + 1);
               if (j < dimens[i] - 1)
                    parse_char(&input->locus, ';');
          }
     }

     return array;
}



void *XCAT(parse_array_, FILE_NAME)(input_data *input, int n, int *dimens) {

     if (! input->use_array2)
          return XCAT(FILE_NAME, _array_parse1)(input, n, dimens, 0);
     else
          return XCAT(FILE_NAME, _array_parse2)(input, n, dimens, 0);
}



C_TYPE *XCAT(parse_array1_, FILE_NAME)(input_data *input, int n) {

     int dimens[1];

     dimens[0] = n;

     return (C_TYPE *) XCAT(parse_array_, FILE_NAME)(input, 1, dimens);
}



C_TYPE **XCAT(parse_array2_, FILE_NAME)(input_data *input, int m, int n) {

     int dimens[2];

     dimens[0] = m;
     dimens[1] = n;

     return (C_TYPE **) XCAT(parse_array_, FILE_NAME)(input, 2, dimens);
}



C_TYPE ***XCAT(parse_array3_, FILE_NAME)(input_data *input, int m, int n, int o) {

     int dimens[3];

     dimens[0] = m;
     dimens[1] = n;
     dimens[2] = o;

     return (C_TYPE ***) XCAT(parse_array_, FILE_NAME)(input, 3, dimens);
}



C_TYPE ****XCAT(parse_array4_, FILE_NAME)(input_data *input, int m, int n, int o, int p) {

     int dimens[4];

     dimens[0] = m;
     dimens[1] = n;
     dimens[2] = o;
     dimens[3] = p;

     return (C_TYPE ****) XCAT(parse_array_, FILE_NAME)(input, 4, dimens);
}



input_type_data XCAT(parse_, FILE_NAME)(input_data *input, int n, int *dimens) {

     input_type_data input_type;

     input_type.name  = NULL;
     input_type.flag  = 0;
     input_type.type  = FILE_TYPE;
     input_type.order = n;
     copy_dimens(input_type.dimens, dimens, n);

     if (input_type.order <= 0)
          input_type.d.FILE_MEMBER = XCAT(parse_scaler_, FILE_NAME)(input);
     else
          input_type.d.v           = XCAT(parse_array_,  FILE_NAME)(input, n, dimens);

     return input_type;
}



/*******************************************************************************
 *
 ******************************************************************************/
void XCAT(free_scaler_, FILE_NAME)(C_TYPE x) {
#ifdef FREE_SCALER
     free(x);
#endif
}



void XCAT(FILE_NAME, _array_free)(void *array, int n, int *dimens, int i) {

     int j;

#ifndef FREE_SCALER
     if (i < n - 1) {
#endif
          for (j = 0; j < dimens[i]; ++j) {
               if (i == n - 1)
                    XCAT(free_scaler_, FILE_NAME)(((C_TYPE *) array)[j]);
               else {
                    XCAT(FILE_NAME, _array_free)(((void **) array)[j], n, dimens, i + 1);
               }
          }
#ifndef FREE_SCALER
     }
#endif
     free(array);
}



void XCAT(free_array_, FILE_NAME)(void *array, int n, int *dimens) {

     XCAT(FILE_NAME, _array_free)(array, n, dimens, 0);
}



void XCAT(free_array1_, FILE_NAME)(C_TYPE *array, int n) {

     int dimens[1];

     dimens[0] = n;

     XCAT(free_array_, FILE_NAME)(array, 1, dimens);
}



void XCAT(free_array2_, FILE_NAME)(C_TYPE **array, int m, int n) {

     int dimens[2];

     dimens[0] = m;
     dimens[1] = n;

     XCAT(free_array_, FILE_NAME)(array, 2, dimens);
}



void XCAT(free_array3_, FILE_NAME)(C_TYPE ***array, int m, int n, int o) {

     int dimens[3];

     dimens[0] = m;
     dimens[1] = n;
     dimens[2] = o;

     XCAT(free_array_, FILE_NAME)(array, 3, dimens);
}



void XCAT(free_array4_, FILE_NAME)(C_TYPE ****array, int m, int n, int o, int p) {

     int dimens[4];

     dimens[0] = m;
     dimens[1] = n;
     dimens[2] = o;
     dimens[3] = o;

     XCAT(free_array_, FILE_NAME)(array, 4, dimens);
}



void XCAT(free_, FILE_NAME)(input_type_data *input_type) {

     if (input_type->order == 0)
          XCAT(free_scaler_, FILE_NAME)(input_type->d.FILE_MEMBER);
     else
          XCAT(free_array_, FILE_NAME)(input_type->d.v, input_type->order, input_type->dimens);
}



/*******************************************************************************
 *
 ******************************************************************************/
void XCAT(print_scaler_, FILE_NAME)(FILE *fp, C_TYPE x) {

     fprintf(fp, FORMAT, x);
}



void XCAT(FILE_NAME, _array_print1)(FILE *fp, void *array, int n, int *dimens, int i) {

     int j;

     fprintf(fp, "{");

     for (j = 0; j < dimens[i]; ++j) {
          if (i == n - 1)
               XCAT(print_scaler_, FILE_NAME)(fp, ((C_TYPE *) array)[j]);
          else {
               XCAT(FILE_NAME, _array_print1)(fp, ((void **) array)[j], n, dimens, i + 1);
          }

          if (j < dimens[i] - 1)
               fprintf(fp, ", ");

     }

     fprintf(fp, "}");
}



void XCAT(FILE_NAME, _array_print2)(FILE *fp, void *array, int n, int *dimens, int i) {

     int j;

     for (j = 0; j < dimens[i]; ++j) {
          if (i == n - 1) {
               XCAT(print_scaler_, FILE_NAME)(fp, ((C_TYPE *) array)[j]);
               if (j < dimens[i] - 1)
                     fprintf(fp, ", ");
          }
          else {
               XCAT(FILE_NAME, _array_print2)(fp, ((void **) array)[j], n, dimens, i + 1);
               if (j < dimens[i] - 1)
                    fprintf(fp, "; ");
          }
     }
}



void XCAT(print_array_, FILE_NAME)(FILE *fp, format_data *d, void *array, int n, int *dimens) {

     if (! d->use_array2)
          XCAT(FILE_NAME, _array_print1)(fp, array, n, dimens, 0);
     else
          XCAT(FILE_NAME, _array_print2)(fp, array, n, dimens, 0);
}



void XCAT(print_array1_, FILE_NAME)(FILE *fp, format_data *d, C_TYPE *array, int n) {

     int dimens[1];

     dimens[0] = n;

     XCAT(print_array_, FILE_NAME)(fp, d, array, 1, dimens);
}



void XCAT(print_array2_, FILE_NAME)(FILE *fp, format_data *d, C_TYPE **array, int m, int n) {

     int dimens[2];

     dimens[0] = m;
     dimens[1] = n;

     XCAT(print_array_, FILE_NAME)(fp, d, array, 2, dimens);
}



void XCAT(print_array3_, FILE_NAME)(FILE *fp, format_data *d, C_TYPE ***array, int m, int n, int o) {

     int dimens[3];

     dimens[0] = m;
     dimens[1] = n;
     dimens[2] = o;

     XCAT(print_array_, FILE_NAME)(fp, d, array, 3, dimens);
}



void XCAT(print_array4_, FILE_NAME)(FILE *fp, format_data *d, C_TYPE ****array, int m, int n, int o, int p) {

     int dimens[4];

     dimens[0] = m;
     dimens[1] = n;
     dimens[2] = o;
     dimens[3] = p;

     XCAT(print_array_, FILE_NAME)(fp, d, array, 4, dimens);
}



void XCAT(print_, FILE_NAME)(FILE *fp, format_data *d, input_type_data *input_type) {

     if (input_type->order == 0)
          XCAT(print_scaler_, FILE_NAME)(fp, input_type->d.FILE_MEMBER);
     else
          XCAT(print_array_, FILE_NAME)(fp, d, input_type->d.v, input_type->order, input_type->dimens);
}
