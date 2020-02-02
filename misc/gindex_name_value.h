/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef GINDEX_NAME_VALUE_H
#define GINDEX_NAME_VALUE_H

#ifdef __cplusplus
extern "C" {
#endif


#define STATIC_NO
#define STATIC_YES static


#define TEMPLATE_GINDEX_NAME_VALUE(NAME, NAME_STRING, N_VALUES, STATIC)		\
STATIC int XCAT(NAME, _n)() {							\
										\
     return N_VALUES;								\
}										\
										\
STATIC int XCAT(NAME, _max_name_length)() {					\
										\
     return max_name_length(XCAT(NAME, _names), N_VALUES);			\
}										\
										\
STATIC const char *XCAT(NAME, _index_to_name)(int index) {			\
										\
     return index_to_name(index, XCAT(NAME, _names), N_VALUES, NAME_STRING);	\
}										\
										\
STATIC int XCAT(NAME, _name_to_index)(const char *name) {			\
										\
     return name_to_index(name, XCAT(NAME, _names), N_VALUES, NAME_STRING);	\
}										\
										\
STATIC enum XCAT(NAME, _type) XCAT(NAME, _index_to_value)(int index) {		\
										\
     return (enum XCAT(NAME, _type)) index_to_value(index, XCAT(NAME, _types), N_VALUES, NAME_STRING);				\
}										\
										\
STATIC int XCAT(NAME, _value_to_index)(enum XCAT(NAME, _type) type) {		\
										\
     return value_to_index(type, XCAT(NAME, _types), N_VALUES, NAME_STRING);	\
}										\
										\
STATIC enum XCAT(NAME, _type) XCAT(NAME, _name_to_value)(const char *name) {	\
										\
     return (enum XCAT(NAME, _type)) name_to_value(name, XCAT(NAME, _types), XCAT(NAME, _names), N_VALUES, NAME_STRING);	\
}										\
										\
STATIC const char *XCAT(NAME, _value_to_name)(enum XCAT(NAME, _type) type) {	\
										\
     return value_to_name(type, XCAT(NAME, _types), XCAT(NAME, _names), N_VALUES, NAME_STRING);					\
}


#define GINDEX_NAME_VALUE_TEMPLATE(NAME, NAME_STRING, N_VALUES)			\
     TEMPLATE_GINDEX_NAME_VALUE(NAME, NAME_STRING, N_VALUES, STATIC_NO)


#define GINDEX_NAME_VALUE_TEMPLATE_STATIC(NAME, NAME_STRING, N_VALUES)		\
     TEMPLATE_GINDEX_NAME_VALUE(NAME, NAME_STRING, N_VALUES, STATIC_YES)


#define TEMPLATE_GINDEX_NAME_MASK(NAME, NAME_STRING, N_VALUES, STATIC)		\
STATIC int XCAT(NAME, _n)() {							\
										\
     return N_VALUES;								\
}										\
										\
STATIC int XCAT(NAME, _max_name_length)() {					\
										\
     return max_name_length(XCAT(NAME, _names), N_VALUES);			\
}										\
										\
STATIC const char *XCAT(NAME, _index_to_name)(int index) {			\
										\
     return index_to_name(index, XCAT(NAME, _names), N_VALUES, NAME_STRING);	\
}										\
										\
STATIC int XCAT(NAME, _name_to_index)(const char *name) {			\
										\
     return name_to_index(name, XCAT(NAME, _names), N_VALUES, NAME_STRING);	\
}										\
										\
STATIC enum XCAT(NAME, _mask) XCAT(NAME, _index_to_mask)(int index) {		\
										\
     return (enum XCAT(NAME, _mask)) index_to_value(index, XCAT(NAME, _masks), N_VALUES, NAME_STRING);				\
}										\
										\
STATIC int XCAT(NAME, _mask_to_index)(enum XCAT(NAME, _mask) mask) {		\
										\
     return value_to_index(mask, XCAT(NAME, _masks), N_VALUES, NAME_STRING);	\
}										\
										\
STATIC int XCAT(NAME, _mask_to_index_list)(enum XCAT(NAME, _mask) mask, int *indexes, int length) {				\
										\
     return mask_to_index_list(mask, XCAT(NAME, _masks), N_VALUES, indexes, length);						\
}										\
										\
STATIC enum XCAT(NAME, _mask) XCAT(NAME, _name_to_mask)(const char *name) {	\
										\
     return (enum XCAT(NAME, _mask)) name_to_value(name, XCAT(NAME, _masks), XCAT(NAME, _names), N_VALUES, NAME_STRING);	\
}										\
										\
STATIC const char *XCAT(NAME, _mask_to_name)(enum XCAT(NAME, _mask) mask) {	\
										\
     return value_to_name(mask, XCAT(NAME, _masks), XCAT(NAME, _names), N_VALUES, NAME_STRING);					\
}										\
										\
STATIC char *XCAT(NAME, _mask_to_name_list)(enum XCAT(NAME, _mask) mask, char *s, int length) {					\
										\
     return mask_to_name_list(mask, XCAT(NAME, _masks), N_VALUES, XCAT(NAME, _names), s, length);				\
}										\
										\
STATIC enum XCAT(NAME, _mask) XCAT(NAME, _name_list_to_mask)(char *s) {				\
										\
     return (enum XCAT(NAME, _mask)) name_list_to_mask(s, XCAT(NAME, _masks), XCAT(NAME, _names), N_VALUES, NAME_STRING);				\
}										\
										\
STATIC int XCAT(NAME, _mask_to_value_list)(enum XCAT(NAME, _mask) mask, long *values, int length) {				\
										\
     return mask_to_value_list(mask, XCAT(NAME, _masks), N_VALUES, values, length);						\
}										\
										\
STATIC enum XCAT(NAME, _mask) XCAT(NAME, _value_list_to_mask)(long *values, int n_values) {							\
										\
     return (enum XCAT(NAME, _mask)) value_list_to_mask(values, n_values);				\
}


#define GINDEX_NAME_MASK_TEMPLATE(NAME, NAME_STRING, N_VALUES)			\
     TEMPLATE_GINDEX_NAME_MASK(NAME, NAME_STRING, N_VALUES, STATIC_NO)


#define GINDEX_NAME_MASK_TEMPLATE_STATIC(NAME, NAME_STRING, N_VALUES)		\
     TEMPLATE_GINDEX_NAME_MASK(NAME, NAME_STRING, N_VALUES, STATIC_YES)


int max_name_length(const char **names, int n);
const char *index_to_name(int index, const char **names, int n, const char *desc);
int name_to_index(const char *name, const char **names, int n, const char *desc);
long index_to_value(int index, long *values, int n, const char *desc);
int value_to_index(long value, long *values, int n, const char *desc);
int mask_to_index_list(long mask, long *values, int n_values, int *indexes, int length);
long name_to_value(const char *name, long *values, const char **names, int n, const char *desc);
const char *value_to_name(long value, long *values, const char **names, int n, const char *desc);
char *mask_to_name_list(long mask, long *values, int n_values, const char **names, char *s, int length);
long name_list_to_mask(char *list, long *values, const char **names, int n, const char *desc);
int mask_to_value_list(long mask, long *values, int n_values, long *values2, int length);
long value_list_to_mask(long *values, int n_values);


#ifdef __cplusplus
}
#endif

#endif /* GINDEX_NAME_VALUE_H */
