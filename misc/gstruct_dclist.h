/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


#define dclist                    XCAT(GSTRUCT_PREFIX,dclist)
#define dcelem                    XCAT(GSTRUCT_PREFIX,dcelem)
#define dclist_type               XCAT(GSTRUCT_PREFIX,dclist_type)
#define dclist_node               XCAT(GSTRUCT_PREFIX,dclist_node)

#define dclist_create             XCAT(GSTRUCT_PREFIX,dclist_create)
#define dclist_empty              XCAT(GSTRUCT_PREFIX,dclist_empty)
#define dclist_free               XCAT(GSTRUCT_PREFIX,dclist_free)
#define dclist_set_comp           XCAT(GSTRUCT_PREFIX,dclist_set_comp)
#define dclist_set_free           XCAT(GSTRUCT_PREFIX,dclist_set_free)
#define dclist_set_print          XCAT(GSTRUCT_PREFIX,dclist_set_print)
#define dclist_is_empty           XCAT(GSTRUCT_PREFIX,dclist_is_empty)
#define dclist_count              XCAT(GSTRUCT_PREFIX,dclist_count)
#define dclist_print              XCAT(GSTRUCT_PREFIX,dclist_print)
#define dclist_insert             XCAT(GSTRUCT_PREFIX,dclist_insert)
#define dclist_insert_end         XCAT(GSTRUCT_PREFIX,dclist_insert_end)
#define dclist_insert_front       XCAT(GSTRUCT_PREFIX,dclist_insert_front)
#define dclist_insert_sort_type   XCAT(GSTRUCT_PREFIX,dclist_insert_sort_type)
#define dclist_insert_sort_comp   XCAT(GSTRUCT_PREFIX,dclist_insert_sort_comp)
#define dclist_add                XCAT(GSTRUCT_PREFIX,dclist_add)
#define dclist_add_end            XCAT(GSTRUCT_PREFIX,dclist_add_end)
#define dclist_add_front          XCAT(GSTRUCT_PREFIX,dclist_add_front)
#define dclist_add_sort_type      XCAT(GSTRUCT_PREFIX,dclist_add_sort_type)
#define dclist_add_sort_comp      XCAT(GSTRUCT_PREFIX,dclist_add_sort_comp)
#define dclist_sort_type          XCAT(GSTRUCT_PREFIX,dclist_sort_type)
#define dclist_sort_comp          XCAT(GSTRUCT_PREFIX,dclist_sort_comp)
#define dclist_delete             XCAT(GSTRUCT_PREFIX,dclist_delete)
#define dclist_delete_type        XCAT(GSTRUCT_PREFIX,dclist_delete_type)
#define dclist_delete_comp        XCAT(GSTRUCT_PREFIX,dclist_delete_comp)
#define dclist_remove             XCAT(GSTRUCT_PREFIX,dclist_remove)
#define dclist_remove_type        XCAT(GSTRUCT_PREFIX,dclist_remove_type)
#define dclist_remove_comp        XCAT(GSTRUCT_PREFIX,dclist_remove_comp)
#define dclist_find_type          XCAT(GSTRUCT_PREFIX,dclist_find_type)
#define dclist_find_comp          XCAT(GSTRUCT_PREFIX,dclist_find_comp)
#define dclist_first              XCAT(GSTRUCT_PREFIX,dclist_first)
#define dclist_last               XCAT(GSTRUCT_PREFIX,dclist_last)
#define dcelem_retreat            XCAT(GSTRUCT_PREFIX,dcelem_retreat)
#define dcelem_advance            XCAT(GSTRUCT_PREFIX,dcelem_advance)
#define dcelem_value              XCAT(GSTRUCT_PREFIX,dcelem_value)
#define dcelem_pointer            XCAT(GSTRUCT_PREFIX,dcelem_pointer)


#ifdef GSTRUCT_POINTER_TYPE
#define POINTER(a)  (a)
#else
#define POINTER(a) &(a)
#endif


#if defined(DECLARATION) || defined(DEC_BEFORE)

typedef struct dclist_type dclist;
typedef struct dclist_node dcelem;

#endif


#if defined(DECLARATION) || defined(DEC_AFTER)

struct dclist_type {
     int (*comp) (const void *, const void *);
     void (*free_) (void *);
     int (*print) (FILE *, const void *);

     dcelem *first;
};

struct dclist_node {
     GSTRUCT_TYPE value;
     dcelem *prev;
     dcelem *next;
};


dclist *dclist_create(void);
void dclist_empty(dclist *);
void dclist_free(dclist *);
void dclist_set_comp(dclist *, int (*)(const void *, const void *));
void dclist_set_free(dclist *, void (*)(void *));
void dclist_set_print(dclist *, int (*)(FILE *, const void *));
int dclist_is_empty(dclist *);
int dclist_count(dclist *);
int dclist_print(dclist *, FILE *);
dcelem *dclist_insert(dclist *, dcelem *, GSTRUCT_TYPE);
dcelem *dclist_insert_end(dclist *, GSTRUCT_TYPE);
dcelem *dclist_insert_front(dclist *, GSTRUCT_TYPE);
dcelem *dclist_insert_sort_type(dclist *, GSTRUCT_TYPE);
dcelem *dclist_insert_sort_comp(dclist *, GSTRUCT_TYPE);
void dclist_add(dclist *, dcelem *, dcelem *);
void dclist_add_end(dclist *, dcelem *);
void dclist_add_front(dclist *, dcelem *);
void dclist_add_sort_type(dclist *, dcelem *);
void dclist_add_sort_comp(dclist *, dcelem *);
void dclist_sort_type(dclist *);
void dclist_sort_comp(dclist *);
void dclist_delete(dclist *, dcelem *);
void dclist_delete_type(dclist *, GSTRUCT_TYPE);
void dclist_delete_comp(dclist *, GSTRUCT_TYPE);
dcelem *dclist_remove(dclist *, dcelem *);
dcelem *dclist_remove_type(dclist *, GSTRUCT_TYPE);
dcelem *dclist_remove_comp(dclist *, GSTRUCT_TYPE);
dcelem *dclist_find_type(dclist *, GSTRUCT_TYPE);
dcelem *dclist_find_comp(dclist *, GSTRUCT_TYPE);
dcelem *dclist_first(dclist *);
dcelem *dclist_last(dclist *);
dcelem *dcelem_retreat(dcelem *);
dcelem *dcelem_advance(dcelem *);
GSTRUCT_TYPE dcelem_value(dcelem *);
GSTRUCT_TYPE *dcelem_pointer(dcelem *);

#endif


#ifdef IMPLEMENTATION

dclist *dclist_create(void) {

     dclist *list;

     list = (dclist *) malloc(sizeof(dclist));

     if (list == NULL)
          fprintf(stderr, "Memory allocation error: dclist_create()");

     list->comp   = NULL;
     list->free_  = NULL;
     list->print  = NULL;

     list->first  = NULL;

     return list;
}



void dclist_empty(dclist *list) {

     dcelem *pos;
     dcelem *tmp;

     if (list->first) {
          pos = list->first;

          do {
               tmp = pos->next;
               if (list->free_)
                    list->free_(POINTER(pos->value));
               free(pos);
               pos = tmp;
          } while (pos != list->first);

          list->first = NULL;
     }
}



void dclist_free(dclist *list) {

     dclist_empty(list);

     free(list);
}



void dclist_set_comp(dclist *list, int (*comp)(const void *, const void *)) {

     list->comp = comp;
}



void dclist_set_free(dclist *list, void (*free)(void *)) {

     list->free_ = free;
}



void dclist_set_print(dclist *list, int (*print)(FILE *, const void *)) {

     list->print = print;
}



int dclist_is_empty(dclist *list) {

     return list->first == NULL;
}



int dclist_count(dclist *list) {

     int count = 0;

     dcelem *pos;

     if (list->first) {
          pos = list->first;

          do {
               pos = pos->next;
               ++count;
          } while (pos != list->first);
     }

     return count;
}



int dclist_print(dclist *list, FILE *fp) {

     int count = 0;

     dcelem *pos;

     if (list->first) {
          pos = list->first;

          do {
               list->print(fp, POINTER(pos->value));
               pos = pos->next;

               ++count;
          } while (pos != list->first);
     }

     return count;
}



dcelem *dclist_insert(dclist *list, dcelem *pos, GSTRUCT_TYPE value) {

     dcelem *tmp;

     tmp = (dcelem *) malloc(sizeof(dcelem));

     if (tmp == NULL)
          fprintf(stderr, "Memory allocation error: dclist_insert()");

     tmp->value = value;

     if (list->first == NULL) {
          list->first = tmp;
          list->first->prev = tmp;
          list->first->next = tmp;
     }
     else {
          tmp->prev = pos;
          tmp->next = pos->next;
          tmp->prev->next = tmp;
          tmp->next->prev = tmp;
     }

     return tmp;
}



dcelem *dclist_insert_end(dclist *list, GSTRUCT_TYPE value) {

     return dclist_insert(list, dclist_last(list), value);
}



dcelem *dclist_insert_front(dclist *list, GSTRUCT_TYPE value) {

     dcelem *pos;
     dcelem *tmp;

     pos = dclist_last(list);

     tmp = dclist_insert(list, pos, value);

     list->first = pos->next;

     return tmp;
}



#ifdef GSTRUCT_BASIC_TYPE

dcelem *dclist_insert_sort_type(dclist *list, GSTRUCT_TYPE value) {

     dcelem *pos;

     /* First check for the case where the user is inserting presorted data.  In
        this case a doubly-linked circular list can efficiently insert after the
        last element. */
     pos = dclist_last(list);

     if (pos) {
          if (pos->value <= value) {
               return dclist_insert(list, pos, value);
          }
     }

     pos = list->first;

     if (pos) {
          if (value < pos->value) {
               return dclist_insert_front(list, value);
          }

          do {
               if (value < pos->next->value)
                    break;

               pos = pos->next;
          } while (pos->next != list->first);
     }

     return dclist_insert(list, pos, value);
}

#endif



dcelem *dclist_insert_sort_comp(dclist *list, GSTRUCT_TYPE value) {

     dcelem *pos;

     /* First check for the case where the user is inserting presorted data.  In
        this case a doubly-linked circular list can efficiently insert after the
        last element. */
     pos = dclist_last(list);

     if (pos) {
          if (list->comp(POINTER(pos->value), POINTER(value)) <= 0) {
               return dclist_insert(list, pos, value);
          }
     }

     pos = list->first;

     if (pos) {
          if (list->comp(POINTER(value), POINTER(pos->value)) < 0) {
               return dclist_insert_front(list, value);
          }

          do {
               if (list->comp(POINTER(value), POINTER(pos->next->value)) < 0)
                    break;

               pos = pos->next;
          } while (pos->next != list->first);
     }

     return dclist_insert(list, pos, value);
}



void dclist_add(dclist *list, dcelem *pos, dcelem *elem) {

     if (list->first == NULL) {
          list->first = elem;
          list->first->prev = elem;
          list->first->next = elem;
     }
     else {
          elem->prev = pos;
          elem->next = pos->next;
          elem->prev->next = elem;
          elem->next->prev = elem;
     }
}



void dclist_add_end(dclist *list, dcelem *elem) {

     dclist_add(list, dclist_last(list), elem);
}



void dclist_add_front(dclist *list, dcelem *elem) {

     dcelem *pos;

     pos = dclist_last(list);

     dclist_add(list, pos, elem);

     list->first = pos->next;
}



#ifdef GSTRUCT_BASIC_TYPE

void dclist_add_sort_type(dclist *list, dcelem *elem) {

     dcelem *pos;

     pos = list->first;

     if (pos) {
          if (elem->value < pos->value) {
               dclist_add_front(list, elem);
               return;
          }

          do {
               if (elem->value < pos->next->value)
                    break;

               pos = pos->next;
          } while (pos->next != list->first);
     }

     dclist_add(list, pos, elem);
}

#endif



void dclist_add_sort_comp(dclist *list, dcelem *elem) {

     dcelem *pos;

     pos = list->first;

     if (pos) {
          if (list->comp(POINTER(elem->value), POINTER(pos->value)) < 0) {
               dclist_add_front(list, elem);
               return;
          }

          do {
               if (list->comp(POINTER(elem->value), POINTER(pos->next->value)) < 0)
                    break;

               pos = pos->next;
          } while (pos->next != list->first);
     }

     dclist_add(list, pos, elem);
}



#ifdef GSTRUCT_BASIC_TYPE

void dclist_sort_type(dclist *list) {

     dcelem *pos;
     dcelem *tmp;

     if (list->first) {
          pos = list->first;

          while (pos->next != list->first) {
               if (pos->value > pos->next->value) {
                    tmp = dclist_remove(list, pos->next);
                    dclist_add_sort_type(list, tmp);
               }
               else
                    pos = pos->next;
          }
     }
}

#endif



void dclist_sort_comp(dclist *list) {

     dcelem *pos;
     dcelem *tmp;

     if (list->first) {
          pos = list->first;

          while (pos->next != list->first) {
               if (list->comp(POINTER(pos->value), POINTER(pos->next->value)) > 0) {
                    tmp = dclist_remove(list, pos->next);
                    dclist_add_sort_comp(list, tmp);
               }
               else
                    pos = pos->next;
          }
     }
}



void dclist_delete(dclist *list, dcelem *pos) {

     if (pos == pos->next) {
          list->first = NULL;
     }
     else {
          if (pos == list->first)
              list->first = pos->next;
          pos->prev->next = pos->next;
          pos->next->prev = pos->prev;
     }

     if (list->free_)
          list->free_(POINTER(pos->value));

     free(pos);
}



#ifdef GSTRUCT_BASIC_TYPE

void dclist_delete_type(dclist *list, GSTRUCT_TYPE value) {

     dcelem *pos;

     if ((pos = dclist_find_type(list, value)))
          dclist_delete(list, pos);
}

#endif



void dclist_delete_comp(dclist *list, GSTRUCT_TYPE value) {

     dcelem *pos;

     if ((pos = dclist_find_comp(list, value)))
          dclist_delete(list, pos);
}



dcelem *dclist_remove(dclist *list, dcelem *pos) {

     if (pos == pos->next) {
          list->first = NULL;
     }
     else {
          if (pos == list->first)
              list->first = pos->next;
          pos->prev->next = pos->next;
          pos->next->prev = pos->prev;
     }

     return pos;
}



#ifdef GSTRUCT_BASIC_TYPE

dcelem *dclist_remove_type(dclist *list, GSTRUCT_TYPE value) {

     dcelem *pos;

     if ((pos = dclist_find_type(list, value)))
          return dclist_remove(list, pos);
     else
          return NULL;
}

#endif



dcelem *dclist_remove_comp(dclist *list, GSTRUCT_TYPE value) {

     dcelem *pos;

     if ((pos = dclist_find_comp(list, value)))
          return dclist_remove(list, pos);
     else
          return NULL;
}



#ifdef GSTRUCT_BASIC_TYPE

dcelem *dclist_find_type(dclist *list, GSTRUCT_TYPE value) {

     int end;

     dcelem *pos;

     pos = list->first;

     end = 1;

     if (pos) {
          do {
               if (pos->value == value) {
                    end = 0;
                    break;
               }
               pos = pos->next;
          } while (pos != list->first);
     }

     if (end)
          return NULL;

     return pos;
}

#endif



dcelem *dclist_find_comp(dclist *list, GSTRUCT_TYPE value) {

     int end;

     dcelem *pos;

     pos = list->first;

     end = 1;

     if (pos) {
          do {
               if (list->comp(POINTER(value), POINTER(pos->value)) == 0) {
                    end = 0;
                    break;
               }
               pos = pos->next;
          } while (pos != list->first);
     }

     if (end)
          return NULL;

     return pos;
}



dcelem *dclist_first(dclist *list) {

     return list->first;
}



dcelem *dclist_last(dclist *list) {

     if (list->first == NULL)
          return NULL;

     return list->first->prev;
}



dcelem *dcelem_retreat(dcelem *pos) {

     return pos->prev;
}



dcelem *dcelem_advance(dcelem *pos) {

     return pos->next;
}



GSTRUCT_TYPE dcelem_value(dcelem *pos) {

     return pos->value;
}



GSTRUCT_TYPE *dcelem_pointer(dcelem *pos) {

     return &pos->value;
}

#endif


#undef dclist
#undef dcelem
#undef dclist_type
#undef dclist_node

#undef dclist_create
#undef dclist_empty
#undef dclist_free
#undef dclist_set_comp
#undef dclist_set_free
#undef dclist_set_print
#undef dclist_is_empty
#undef dclist_count
#undef dclist_print
#undef dclist_insert
#undef dclist_insert_end
#undef dclist_insert_front
#undef dclist_insert_sort_type
#undef dclist_insert_sort_comp
#undef dclist_add
#undef dclist_add_end
#undef dclist_add_front
#undef dclist_add_sort_type
#undef dclist_add_sort_comp
#undef dclist_sort_type
#undef dclist_sort_comp
#undef dclist_delete
#undef dclist_delete_type
#undef dclist_delete_comp
#undef dclist_remove
#undef dclist_remove_type
#undef dclist_remove_comp
#undef dclist_find_type
#undef dclist_find_comp
#undef dclist_first
#undef dclist_last
#undef dcelem_retreat
#undef dcelem_advance
#undef dcelem_value
#undef dcelem_pointer

#undef POINTER
