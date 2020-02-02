/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"

#include "glist.h"


#ifndef __cplusplus


/*******************************************************************************
 *
 ******************************************************************************/
void list_init(void *v) {

     struct list_data *head;

     head = (struct list_data *) v;

     head->next = head;
     head->prev = head;
}



int list_is_empty(const void *v) {

     struct list_data *head;

     head = (struct list_data *) v;

     return head->next == head;
}



int list_count(const void *v) {

     int count = 0;

     struct list_data *head;
     struct list_data *elem;

     head = (struct list_data *) v;

     list_for_each(head, elem)
          count++;

     return count;
}



int list_is_first_elem(const void *v1, const void *v2) {

     struct list_data *head;
     struct list_data *elem;

     head = (struct list_data *) v1;
     elem = (struct list_data *) v2;

     return elem == head->next;
}



int list_is_last_elem(const void *v1, const void *v2) {

     struct list_data *head;
     struct list_data *elem;

     head = (struct list_data *) v1;
     elem = (struct list_data *) v2;

     return elem->next == head;
}



void *list_last_elem(const void *v1) {

     struct list_data *head;

     head = (struct list_data *) v1;

     return head->prev;
}



void list_insert(void *v1, void *v2) {

     struct list_data *head;
     struct list_data *elem;

     head = (struct list_data *) v1;
     elem = (struct list_data *) v2;

     head->next->prev = elem;
     elem->next = head->next;
     elem->prev = head;
     head->next = elem;
}



void *list_find(const void *v, const char *name) {

     struct list_data *head;
     struct list_data *elem;

     head = (struct list_data *) v;

     list_for_each(head, elem) {
          if (strcmp(elem->name, name) == 0)
               return elem;
     }

     return NULL;
}



void *list_append(void *v1, void *v2, int flag) {

     struct list_data *head;
     struct list_data *elem;

     head = (struct list_data *) v1;
     elem = (struct list_data *) v2;

     if (flag) {
          if (list_find(head, elem->name))
               return NULL;
     }

     list_insert(head->prev, elem);

     return elem;
}



void list_free(void *v) {

     struct list_data *head;
     struct list_data *elem;
     struct list_data *next;

     head = (struct list_data *) v;

     list_for_each_safe(head, elem, next)
          free(elem);
}

#endif
