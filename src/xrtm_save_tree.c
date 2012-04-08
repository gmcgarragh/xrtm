/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>
#include <gstruct.h>

#include "xrtm.h"
#include "xrtm_save_tree.h"


static ulong djb2(ulong hash, const uchar *str) {
/*
     ulong hash = 5381;
*/
     int c;

     while ((c = *str++))
         hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

     return hash;
}



static int save_tree_comp(save_node_data *d1, save_node_data *d2) {

     if (d1->hash < d2->hash)
          return -1;
     if (d1->hash > d2->hash)
          return 1;

     return 0;
}



static void save_node_free(save_node_data *d) {

     forward_save_data *d2 = (forward_save_data *) d->d;

     if (! d->is_proxy) {
          if (d2->free)
               d2->free(d2);

          free(d2);
     }

     free(d);
}



int save_tree_init(save_tree_data *d) {

     d->temp = (char *) malloc(1024);

     d->i = 0;

     d->hash = (ulong *) malloc(128 * sizeof(ulong));
     d->hash[0] = 5381;

     d->t = gaatree_create();
     gaatree_set_comp(d->t, (int (*)(const void *, const void *)) save_tree_comp);
     gaatree_set_free(d->t, (void (*)(void *)) save_node_free);

     return 0;
}



void save_tree_free(save_tree_data *d) {

     free(d->temp);

     free(d->hash);

     gaatree_free(d->t);
}



void save_tree_encode_s(save_tree_data *d, const char *s) {

     d->hash[d->i + 1] = djb2(d->hash[d->i], (uchar *) s);
     d->i++;
}



void save_tree_encode_i(save_tree_data *d, int i) {
/*
     sprintf(d->temp, "%d", i);
*/
     int ii;

     d->temp[0] = 48 + i % 10;
     for (ii = 1; i /= 10; ++ii)
          d->temp[ii] = 48 + i % 10;
     d->temp[ii] = '\0';

     save_tree_encode_s(d, d->temp);
}



void save_tree_encode_i_j(save_tree_data *d, int i, int j) {

     sprintf(d->temp, "%d_%d", i, j);

     save_tree_encode_s(d, d->temp);
}



void save_tree_decode_s(save_tree_data *d, const char *s) {

     d->i--;
}



void save_tree_decode_i(save_tree_data *d) {

     d->i--;
}



void save_tree_decode_i_j(save_tree_data *d) {

     d->i--;
}



void save_tree_recode_s(save_tree_data *d, const char *s) {

     save_tree_decode_s(d, s);

     save_tree_encode_s(d, s);
}



void save_tree_recode_i(save_tree_data *d, int i, int flag) {

     if (! flag)
          save_tree_decode_i(d);

     save_tree_encode_i(d, i);
}



void save_tree_recode_i_j(save_tree_data *d, int i, int j, int flag) {

     if (! flag)
          save_tree_decode_i(d);

     save_tree_encode_i_j(d, i, j);
}



static int save_tree_data_retrieve(gaatree *t, ulong hash, size_t size, void **v) {

     int flag = 0;

     gaatpos *pos;

     save_node_data node;
     save_node_data *node2;

     node.hash = hash;

     pos = gaatree_find_comp(t, &node);

     if (! pos) {
          flag = 1;

          node2 = (save_node_data *) malloc(sizeof(save_node_data));
          node2->is_proxy = 0;
          node2->hash     = node.hash;
          node2->d        = malloc(size);

          pos = gaatree_insert_comp(t, node2);
     }

     *v = ((save_node_data *) gaatpos_retrieve(pos))->d;

     return flag;
}



int __save_tree_retrieve_data(save_tree_data *d, size_t size, void **v) {

     return save_tree_data_retrieve(d->t, d->hash[d->i], size, v);
}



int __save_tree_retrieve_proxy(save_tree_data *d, const char *s, size_t size, void **v) {

     int flag = 0;

     gaatpos *pos;

     save_node_data node;
     save_node_data *node2;

     node.hash = d->hash[d->i];

     pos = gaatree_find_comp(d->t, &node);

     if (! pos) {
          flag = save_tree_data_retrieve(d->t, djb2(5381, (uchar *) s), size, v);

          node2 = (save_node_data *) malloc(sizeof(save_node_data));
          node2->is_proxy = 1;
          node2->hash     = node.hash;
          node2->d        = *v;

          pos = gaatree_insert_comp(d->t, node2);
     }

     *v = ((save_node_data *) gaatpos_retrieve(pos))->d;

     return flag;
}
