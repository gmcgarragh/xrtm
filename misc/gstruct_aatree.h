/*******************************************************************************
**
**    Copyright (C) 1998-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include "gutil.h"


#define aatree                   XCAT(GSTRUCT_PREFIX,aatree)
#define aatpos                   XCAT(GSTRUCT_PREFIX,aatpos)
#define aatree_type              XCAT(GSTRUCT_PREFIX,aatree_type)
#define aatree_node              XCAT(GSTRUCT_PREFIX,aatree_node)

#define aatree_null_node         XCAT(GSTRUCT_PREFIX,aatree_null_node)

#define aatree_create            XCAT(GSTRUCT_PREFIX,aatree_create)
#define aatree_empty             XCAT(GSTRUCT_PREFIX,aatree_empty)
#define aatree_free              XCAT(GSTRUCT_PREFIX,aatree_free)
#define aatree_set_comp          XCAT(GSTRUCT_PREFIX,aatree_set_comp)
#define aatree_set_free          XCAT(GSTRUCT_PREFIX,aatree_set_free)
#define aatree_set_print         XCAT(GSTRUCT_PREFIX,aatree_set_print)
#define aatree_is_empty          XCAT(GSTRUCT_PREFIX,aatree_is_empty)
#define aatree_count             XCAT(GSTRUCT_PREFIX,aatree_count)
#define aatree_print             XCAT(GSTRUCT_PREFIX,aatree_print)
#define aatree_insert_type       XCAT(GSTRUCT_PREFIX,aatree_insert_type)
#define aatree_insert_comp       XCAT(GSTRUCT_PREFIX,aatree_insert_comp)
#define aatree_delete_type       XCAT(GSTRUCT_PREFIX,aatree_delete_type)
#define aatree_delete_comp       XCAT(GSTRUCT_PREFIX,aatree_delete_comp)
#define aatree_find_type         XCAT(GSTRUCT_PREFIX,aatree_find_type)
#define aatree_find_comp         XCAT(GSTRUCT_PREFIX,aatree_find_comp)
#define aatree_find_min          XCAT(GSTRUCT_PREFIX,aatree_find_min)
#define aatree_find_max          XCAT(GSTRUCT_PREFIX,aatree_find_max)
#define aatpos_retreat_type      XCAT(GSTRUCT_PREFIX,aatpos_retreat_type)
#define aatpos_retreat_comp      XCAT(GSTRUCT_PREFIX,aatpos_retreat_comp)
#define aatpos_advance_type      XCAT(GSTRUCT_PREFIX,aatpos_advance_type)
#define aatpos_advance_comp      XCAT(GSTRUCT_PREFIX,aatpos_advance_comp)
#define aatpos_retrieve          XCAT(GSTRUCT_PREFIX,aatpos_retrieve)
#define aatpos_pointer           XCAT(GSTRUCT_PREFIX,aatpos_pointer)

#define aatpos_empty             XCAT(GSTRUCT_PREFIX,aatpos_empty)
#define aatpos_count             XCAT(GSTRUCT_PREFIX,aatpos_count)
#define aatpos_print             XCAT(GSTRUCT_PREFIX,aatpos_print)
#define aatree_rotate_with_left  XCAT(GSTRUCT_PREFIX,aatree_rotate_with_left)
#define aatree_rotate_with_right XCAT(GSTRUCT_PREFIX,aatree_rotate_with_right)
#define aatree_skew              XCAT(GSTRUCT_PREFIX,aatree_skew)
#define aatree_split             XCAT(GSTRUCT_PREFIX,aatree_split)
#define aatpos_insert_type       XCAT(GSTRUCT_PREFIX,aatpos_insert_type)
#define aatpos_insert_comp       XCAT(GSTRUCT_PREFIX,aatpos_insert_comp)
#define aatpos_rebalance         XCAT(GSTRUCT_PREFIX,aatpos_rebalance)
#define aatpos_delete_type       XCAT(GSTRUCT_PREFIX,aatpos_delete_type)
#define aatpos_delete_comp       XCAT(GSTRUCT_PREFIX,aatpos_delete_comp)
#define aatpos_find_type         XCAT(GSTRUCT_PREFIX,aatpos_find_type)
#define aatpos_find_comp         XCAT(GSTRUCT_PREFIX,aatpos_find_comp)
#define aatpos_find_min          XCAT(GSTRUCT_PREFIX,aatpos_find_min)
#define aatpos_find_max          XCAT(GSTRUCT_PREFIX,aatpos_find_max)
#define aatpos_retback_type      XCAT(GSTRUCT_PREFIX,aatpos_retback_type)
#define aatpos_retback_comp      XCAT(GSTRUCT_PREFIX,aatpos_retback_comp)
#define aatpos_advback_type      XCAT(GSTRUCT_PREFIX,aatpos_advback_type)
#define aatpos_advback_comp      XCAT(GSTRUCT_PREFIX,aatpos_advback_comp)


#ifdef GSTRUCT_POINTER_TYPE
#define POINTER(a)  (a)
#else
#define POINTER(a) &(a)
#endif


#if defined(DECLARATION) || defined(DEC_BEFORE)

typedef struct aatree_type aatree;
typedef struct aatree_node aatpos;

#endif


#if defined(DECLARATION) || defined(DEC_AFTER)

struct aatree_type {
     int (*comp) (const void *, const void *);
     void (*free_) (void *);
     int (*print) (FILE *, const void *);

     aatpos *base;

     aatpos *dele_node;
     aatpos *last_node;
};

struct aatree_node {
     GSTRUCT_TYPE value;
     int    level;
     aatpos *left;
     aatpos *right;
     aatpos *parent;
};


aatree *aatree_create(void);
void aatree_empty(aatree *);
void aatree_free(aatree *);
void aatree_set_comp(aatree *, int (*)(const void *, const void *));
void aatree_set_free(aatree *, void (*)(void *));
void aatree_set_print(aatree *, int (*)(FILE *, const void *));
int aatree_is_empty(aatree *);
int aatree_count(aatree *);
int aatree_print(aatree *, FILE *);
aatpos *aatree_insert_type(aatree *, GSTRUCT_TYPE);
aatpos *aatree_insert_comp(aatree *, GSTRUCT_TYPE);
int aatree_delete_type(aatree *, GSTRUCT_TYPE);
int aatree_delete_comp(aatree *, GSTRUCT_TYPE);
aatpos *aatree_find_type(aatree *, GSTRUCT_TYPE);
aatpos *aatree_find_comp(aatree *, GSTRUCT_TYPE);
aatpos *aatree_find_min(aatree *);
aatpos *aatree_find_max(aatree *);
aatpos *aatpos_retreat_type(aatpos *);
aatpos *aatpos_retreat_comp(aatpos *, int (*)(const void *, const void *));
aatpos *aatpos_advance_type(aatpos *);
aatpos *aatpos_advance_comp(aatpos *, int (*)(const void *, const void *));
GSTRUCT_TYPE aatpos_retrieve(aatpos *);
GSTRUCT_TYPE *aatpos_pointer(aatpos *);

#endif


#ifdef IMPLEMENTATION

static aatpos aatree_null_node = {
     GSTRUCT_NILL,
     0,
     &aatree_null_node,
     &aatree_null_node,
     &aatree_null_node
};


aatree *aatree_create(void) {

     aatree *tree;

     tree = (aatree *) malloc(sizeof(aatree));

     if (tree == NULL)
          fprintf(stderr, "Memory allocation error: aatree_create()");

     tree->comp  = NULL;
     tree->free_ = NULL;
     tree->print = NULL;

     tree->base  = &aatree_null_node;

     return tree;
}



static void aatpos_empty(aatree *tree, aatpos **tpos) {

     if (*tpos != &aatree_null_node) {
          aatpos_empty(tree, &(*tpos)->left);
          aatpos_empty(tree, &(*tpos)->right);

          if (tree->free_)
               tree->free_(POINTER((*tpos)->value));

          free(*tpos);

          *tpos = &aatree_null_node;
     }
}



void aatree_empty(aatree *tree) {

     aatpos_empty(tree, &tree->base);
}



void aatree_free(aatree *tree) {

     aatree_empty(tree);

     free(tree);
}



void aatree_set_comp(aatree *tree, int (*comp)(const void *, const void *)) {

     tree->comp = comp;
}



void aatree_set_free(aatree *tree, void (*free)(void *)) {

     tree->free_ = free;
}



void aatree_set_print(aatree *tree, int (*print)(FILE *, const void *)) {

     tree->print = print;
}



int aatree_is_empty(aatree *tree) {

     return tree->base == &aatree_null_node;
}



static int aatpos_count(aatpos *tpos) {

     int count;

     if (tpos == &aatree_null_node)
          return 0;

     count = 1;

     count += aatpos_count(tpos->left);
     count += aatpos_count(tpos->right);

     return count;
}



int aatree_count(aatree *tree) {

     return aatpos_count(tree->base);
}



static int aatpos_print(aatpos *tpos, FILE *fp, int (*func_print)(FILE *, const void *)) {

     int count;

     if (tpos == &aatree_null_node)
          return 0;

     count = 1;

     count += aatpos_print(tpos->left, fp, func_print);
     func_print(fp, POINTER(tpos->value));
     count += aatpos_print(tpos->right, fp, func_print);

     return count;
}



int aatree_print(aatree *tree, FILE *fp) {

     return aatpos_print(tree->base, fp, tree->print);
}



static aatpos *aatree_rotate_with_left(aatpos *K2) {

     aatpos *K1;

     K1        = K2->left;
     K2->left  = K1->right;
     K1->right =  K2;

     K2->left->parent = K2;
     K1->parent = K1->right->parent;
     K1->right->parent = K1;

     return K1;
}



static aatpos *aatree_rotate_with_right(aatpos *K1) {

     aatpos *K2;

     K2        = K1->right;
     K1->right = K2->left;
     K2->left  = K1;

     K1->right->parent = K1;
     K2->parent = K2->left->parent;
     K2->left->parent = K2;

     return K2;
}



static aatpos *aatree_skew(aatpos *pos) {

     if (pos->left->level == pos->level)
          pos = aatree_rotate_with_left(pos);

     return pos;
}



static aatpos *aatree_split(aatpos *pos) {

     if (pos->right->right->level == pos->level) {
          pos = aatree_rotate_with_right(pos);
          pos->level++;
     }

     return pos;
}



#ifdef GSTRUCT_BASIC_TYPE

static aatpos *aatpos_insert_type(aatpos **tpos_, aatpos *parent, GSTRUCT_TYPE value) {

     aatpos *tpos;
     aatpos *temp;

     tpos = *tpos_;

     if (tpos == &aatree_null_node) {
          tpos = (struct aatree_node *) malloc(sizeof(struct aatree_node));

          if (tpos == NULL)
               fprintf(stderr, "Memory allocation error: aatree_insert_type()");
          else {
               tpos->value  = value;
               tpos->level  = 1;
               tpos->left   = &aatree_null_node;
               tpos->right  = &aatree_null_node;
               tpos->parent = parent;
          }

          temp = tpos;
     }
     else if (value < tpos->value)
          temp = aatpos_insert_type(&tpos->left,  tpos, value);
     else if (value > tpos->value )
          temp = aatpos_insert_type(&tpos->right, tpos, value);
     else
          temp = NULL;

     tpos = aatree_skew (tpos);
     tpos = aatree_split(tpos);

     *tpos_ = tpos;

     return temp;
}



aatpos *aatree_insert_type(aatree *tree, GSTRUCT_TYPE value) {

     return aatpos_insert_type(&tree->base, &aatree_null_node, value);
}

#endif



static aatpos *aatpos_insert_comp(aatpos **tpos_, aatpos *parent, GSTRUCT_TYPE value,
                                  int (*func_comp)(const void *, const void *)) {

     int ret_val;

     aatpos *tpos;
     aatpos *temp;

     tpos = *tpos_;

     if (tpos == &aatree_null_node) {
          tpos = (struct aatree_node *) malloc(sizeof(struct aatree_node));

          if (tpos == NULL)
               fprintf(stderr, "Memory allocation error: aatree_insert_comp()");
          else {
               tpos->value  = value;
               tpos->level  = 1;
               tpos->left   = &aatree_null_node;
               tpos->right  = &aatree_null_node;
               tpos->parent = parent;
          }

          temp = tpos;
     }
     else {
          ret_val = func_comp(POINTER(value), POINTER(tpos->value));

          if (ret_val < 0)
               temp = aatpos_insert_comp(&tpos->left,  tpos, value, func_comp);
          else
          if (ret_val > 0)
               temp = aatpos_insert_comp(&tpos->right, tpos, value, func_comp);
          else
               temp = NULL;
     }

     tpos = aatree_skew (tpos);
     tpos = aatree_split(tpos);

     *tpos_ = tpos;

     return temp;
}



aatpos *aatree_insert_comp(aatree *tree, GSTRUCT_TYPE value) {

     return aatpos_insert_comp(&tree->base, &aatree_null_node, value, tree->comp);
}



static void aatpos_rebalance(aatpos **tpos_) {

     aatpos *tpos;

     tpos = *tpos_;

     if (tpos->left->level < tpos->level - 1 ||
          tpos->right->level < tpos->level - 1) {
          if (tpos->right->level > --tpos->level)
               tpos->right->level = tpos->level;
          tpos = aatree_skew(tpos);
          tpos->right = aatree_skew(tpos->right);
          tpos->right->right = aatree_skew(tpos->right->right);
          tpos = aatree_split(tpos);
          tpos->right = aatree_split(tpos->right);
     }

     *tpos_ = tpos;
}



#ifdef GSTRUCT_BASIC_TYPE

static int aatpos_delete_type(aatree *tree, aatpos **tpos_, GSTRUCT_TYPE value) {

     int ret;

     aatpos *tpos;

     tpos = *tpos_;

     if (tpos != &aatree_null_node) {
          tree->last_node = tpos;
          if (value < tpos->value) {
               ret = aatpos_delete_type(tree, &tpos->left,  value);
          }
          else {
               tree->dele_node = tpos;
               ret = aatpos_delete_type(tree, &tpos->right, value);
          }

          if (tpos == tree->last_node) {
               if (tree->dele_node != &aatree_null_node &&
                   value == tree->dele_node->value) {
                    if (tree->free_)
                         tree->free_(POINTER(tree->dele_node->value));
                    tree->dele_node->value = tpos->value;
                    tree->dele_node = &aatree_null_node;
                    tpos = tpos->right;
                    free(tree->last_node);
                    *tpos_ = tpos;
                    return 1;
               }
               else {
                    *tpos_ = tpos;
                    return 0;
               }
          }
          else {
               aatpos_rebalance(&tpos);
          }

          *tpos_ = tpos;
          return ret;
     }

     *tpos_ = tpos;

     return 0;
}



int aatree_delete_type(aatree *tree, GSTRUCT_TYPE value) {

     return aatpos_delete_type(tree, &tree->base, value);
}

#endif



static int aatpos_delete_comp(aatree *tree, aatpos **tpos_, GSTRUCT_TYPE value) {

     int ret;

     aatpos *tpos;

     tpos = *tpos_;

     if (tpos != &aatree_null_node) {
          tree->last_node = tpos;
          if (tree->comp(POINTER(value), POINTER(tpos->value)) < 0) {
               ret = aatpos_delete_comp(tree, &tpos->left,  value);
          }
          else {
               tree->dele_node = tpos;
               ret = aatpos_delete_comp(tree, &tpos->right, value);
          }

          if (tpos == tree->last_node) {
               if (tree->dele_node != &aatree_null_node &&
                   tree->comp(POINTER(value), POINTER(tree->dele_node->value)) == 0) {
                    if (tree->free_)
                         tree->free_(POINTER(tree->dele_node->value));
                    tree->dele_node->value = tpos->value;
                    tree->dele_node = &aatree_null_node;
                    tpos = tpos->right;
                    free(tree->last_node);
                    *tpos_ = tpos;
                    return 1;
               }
               else {
                    *tpos_ = tpos;
                    return 0;
               }
          }
          else {
               aatpos_rebalance(&tpos);
          }

          *tpos_ = tpos;
          return ret;
     }

     *tpos_ = tpos;

     return 0;
}



int aatree_delete_comp(aatree *tree, GSTRUCT_TYPE value) {

     return aatpos_delete_comp(tree, &tree->base, value);
}



#ifdef GSTRUCT_BASIC_TYPE

static aatpos *aatpos_find_type(aatpos *tpos, GSTRUCT_TYPE value) {

     if (tpos == &aatree_null_node)
          return NULL;

     if (value < tpos->value)
          return aatpos_find_type(tpos->left,  value);
     else
     if (value > tpos->value)
          return aatpos_find_type(tpos->right, value);
     else
          return tpos;
}



aatpos *aatree_find_type(aatree *tree, GSTRUCT_TYPE value) {

     return aatpos_find_type(tree->base, value);
}

#endif



static aatpos *aatpos_find_comp(aatpos *tpos, GSTRUCT_TYPE value,
                                int (*func_comp)(const void *, const void *)) {

     int ret_val;

     if (tpos == &aatree_null_node)
          return NULL;

     ret_val = func_comp(POINTER(value), POINTER(tpos->value));

     if (ret_val < 0)
          return aatpos_find_comp(tpos->left,  value, func_comp);
     if (ret_val > 0)
          return aatpos_find_comp(tpos->right, value, func_comp);

     return tpos;
}



aatpos *aatree_find_comp(aatree *tree, GSTRUCT_TYPE value) {

     return aatpos_find_comp(tree->base, value, tree->comp);
}



static aatpos *aatpos_find_min(aatpos *tpos) {

     if (tpos == &aatree_null_node)
          return NULL;
     if (tpos->left == &aatree_null_node)
          return tpos;

     return aatpos_find_min(tpos->left);
}



aatpos *aatree_find_min(aatree *tree) {

     return aatpos_find_min(tree->base);
}



static aatpos *aatpos_find_max(aatpos *tpos) {

     if (tpos == &aatree_null_node)
          return NULL;
     if (tpos->right == &aatree_null_node)
          return tpos;

     return aatpos_find_max(tpos->right);
}



aatpos *aatree_find_max(aatree *tree) {

     return aatpos_find_max(tree->base);
}



#ifdef GSTRUCT_BASIC_TYPE

static aatpos *aatpos_retback_type(aatpos *tpos) {

     if (tpos->parent != &aatree_null_node) {
          if (tpos->value > tpos->parent->value)
               return tpos->parent;

          return aatpos_retback_type(tpos->parent);
     }

     return NULL;
}



aatpos *aatpos_retreat_type(aatpos *tpos) {

     if (tpos == &aatree_null_node)
          return NULL;
     if (tpos->left != &aatree_null_node)
          return aatpos_find_max(tpos->left);

     return aatpos_retback_type(tpos);
}

#endif



static aatpos *aatpos_retback_comp(aatpos *tpos,
                                   int (*func_comp)(const void *, const void *)) {

     if (tpos->parent != &aatree_null_node) {
          if (func_comp(POINTER(tpos->value), POINTER(tpos->parent->value)) > 0)
               return tpos->parent;

          return aatpos_retback_comp(tpos->parent, func_comp);
     }

     return NULL;
}



aatpos *aatpos_retreat_comp(aatpos *tpos,
                            int (*func_comp)(const void *, const void *)) {

     if (tpos == &aatree_null_node)
          return NULL;
     if (tpos->left != &aatree_null_node)
          return aatpos_find_max(tpos->left);

     return aatpos_retback_comp(tpos, func_comp);
}



#ifdef GSTRUCT_BASIC_TYPE

static aatpos *aatpos_advback_type(aatpos *tpos) {

     if (tpos->parent != &aatree_null_node) {
          if (tpos->value < tpos->parent->value)
               return tpos->parent;

          return aatpos_advback_type(tpos->parent);
     }

     return NULL;
}



aatpos *aatpos_advance_type(aatpos *tpos) {

     if (tpos == &aatree_null_node)
          return NULL;
     if (tpos->right != &aatree_null_node)
          return aatpos_find_min(tpos->right);

     return aatpos_advback_type(tpos);
}

#endif



static aatpos *aatpos_advback_comp(aatpos *tpos,
                                   int (*func_comp)(const void *, const void *)) {

     if (tpos->parent != &aatree_null_node) {
          if (func_comp(POINTER(tpos->value), POINTER(tpos->parent->value)) < 0)
               return tpos->parent;

          return aatpos_advback_comp(tpos->parent, func_comp);
     }

     return NULL;
}



aatpos *aatpos_advance_comp(aatpos *tpos,
                            int (*func_comp)(const void *, const void *)) {

     if (tpos == &aatree_null_node)
          return NULL;
     if (tpos->right != &aatree_null_node)
          return aatpos_find_min(tpos->right);

     return aatpos_advback_comp(tpos, func_comp);
}



GSTRUCT_TYPE aatpos_retrieve(aatpos *tpos) {

     return tpos->value;
}



GSTRUCT_TYPE *aatpos_pointer(aatpos *tpos) {

     return &tpos->value;
}

#endif


#undef aatree
#undef aatpos
#undef aatree_type
#undef aatree_node

#undef aatree_null_node

#undef aatree_create
#undef aatree_empty
#undef aatree_free
#undef aatree_set_comp
#undef aatree_set_free
#undef aatree_set_print
#undef aatree_is_empty
#undef aatree_count
#undef aatree_print
#undef aatree_insert_type
#undef aatree_insert_comp
#undef aatree_delete_type
#undef aatree_delete_comp
#undef aatree_find_type
#undef aatree_find_comp
#undef aatree_find_min
#undef aatree_find_max
#undef aatpos_retreat_type
#undef aatpos_retreat_comp
#undef aatpos_advance_type
#undef aatpos_advance_comp
#undef aatpos_retrieve
#undef aatpos_pointer

#undef aatpos_empty
#undef aatpos_count
#undef aatpos_print
#undef aatree_rotate_with_left
#undef aatree_rotate_with_right
#undef aatree_skew
#undef aatree_split
#undef aatpos_insert_type
#undef aatpos_insert_comp
#undef aatpos_rebalance
#undef aatpos_delete_type
#undef aatpos_delete_comp
#undef aatpos_find_type
#undef aatpos_find_comp
#undef aatpos_find_min
#undef aatpos_find_max
#undef aatpos_retback_type
#undef aatpos_retback_comp
#undef aatpos_advback_type
#undef aatpos_advback_comp

#undef POINTER
