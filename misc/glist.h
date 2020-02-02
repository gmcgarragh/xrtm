/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#ifndef GLIST_H
#define GLIST_H

#ifdef __cplusplus
extern "C" {
#endif


#ifndef __cplusplus

#define list_advance(elem) (elem = (void *) elem->next)

#define list_for_each(head, elem) 				\
     for (elem = (void *) (head)->next; elem != (void *) (head); elem = (void *) elem->next)

#define list_for_each_safe(head, elem, next) 			\
     for (elem = (void *) (head)->next, next = elem->next; elem != (void *) (head); elem = (void *) next, next = (void *) elem->next)


struct list_data {
     char *name;
     struct list_data *prev;
     struct list_data *next;
};


void list_init(void *v);
int list_is_empty(const void *v);
int list_count(const void *v);
int list_is_first_elem(const void *v1, const void *v2);
int list_is_last_elem(const void *v1, const void *v2);
void *list_last_elem(const void *v1);
void list_insert(void *v1, void *v2);
void *list_find(const void *v, const char *name);
void *list_append(void *v1, void *v2, int flag);
void list_free(void *v);

#endif


#ifdef __cplusplus
}
#endif

#endif /* GLIST_H */
