/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/

#include <gutil.h>

#include <setjmp.h>

#include "xrtm_stacks.h"
#include "xrtm_support.h"
#include "xrtm_utility.h"


typedef struct grid_data_ {
     int active;
     int varies;
     struct stack_data_ *s;
} grid_data;



/*******************************************************************************
 *
 ******************************************************************************/
void print_stacks(stack_data *stacks, int n_stacks) {

     int i;

     for (i = 0; i < n_stacks; ++i) {
          printf("%2d, %2d, %2d, %2d, %2d", i,
                 stacks[i].i1, stacks[i].i2,
                 stacks[i].stack, stacks[i].count);

          if (stacks[i].i1 != stacks[i].i2)
               printf(", %2d, %2d, %2d, %2d",
                      stacks[i].p1->i1, stacks[i].p1->i2,
                      stacks[i].p2->i1, stacks[i].p2->i2);

          printf("\n");
     }

     printf("\n");
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static void print_grid(grid_data **grid, int n_layers, int n_derivs) {

     int i;
     int j;

     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_derivs; ++j) {
               if (  grid[i][j].varies)
                    printf(" varied  ");
               else
               if (! grid[i][j].active)
                    printf("combined ");
               else
                    printf("(%02d, %02d) ",
                           grid[i][j].s->i1, grid[i][j].s->i2);
          }
          printf("\n");
     }
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
static int next_layer(grid_data **grid, int i, int j, int n_layers) {

     for ( ; i < n_layers; ++i) {
          if (grid[i][j].varies)
               return -1;
          if (grid[i][j].active)
               return  i;
     }

     return -1;
}



/*******************************************************************************
 *
 ******************************************************************************/
static int max_set_of_pairs1(grid_data **grid, int n_layers, int n_derivs,
                             stack_data **s1, stack_data **s2) {

     int i;
     int j;

     int count;
     int max_c;

     stack_data *cur_s1;
     stack_data *cur_s2;

     stack_data *max_s1;
     stack_data *max_s2;

     count = 0;
     max_c = 0;
     for (j = 0; j < n_derivs; ++j) {
/*
          if (grid[0][j].varies || ! grid[0][j].active)
               goto L1;
*/
          if (  grid[0][j].varies)
               goto L1;

          if (! grid[0][j].active)
               continue;

          i = next_layer(grid, 1, j, n_layers);

          if (i < 0 || grid[i][j].varies)
               goto L1;

          if (count == 0) {
               cur_s1 = grid[0][j].s;
               cur_s2 = grid[i][j].s;
          }

          if (grid[0][j].s != cur_s1 || grid[i][j].s != cur_s2)
               goto L1;

          count++;

          continue;

L1:       if (count > max_c) {
               max_c  = count;
               max_s1 = cur_s1;
               max_s2 = cur_s2;
          }

          count = 0;
     }

     if (count > max_c) {
          max_c  = count;
          max_s1 = cur_s1;
          max_s2 = cur_s2;
     }

     if (max_c != 0) {
          *s1 = max_s1;
          *s2 = max_s2;
     }

     return max_c;
}



/*******************************************************************************
 *
 ******************************************************************************/
/*
static int max_set_of_pairs2(grid_data **grid, int n_layers, int n_derivs,
                             stack_data **s1, stack_data **s2) {

     int i;
     int j;
     int jj;

     int nn;

     int max_c;
     int j_max;

     int *counts;

     stack_data ***index;

     counts = alloc_array1_i(n_derivs);

     init_array1_i(counts, n_derivs, 0);

     index  = (stack_data ***) alloc_array2(2, n_derivs, sizeof(stack_data **));

     nn = 0;
     for (j = 0; j < n_derivs; ++j) {
          if (grid[0][j].varies || ! grid[0][j].active)
               continue;

          i = next_layer(grid, 1, j, n_layers);

          if (i < 0 || grid[i][j].varies)
               continue;

          for (jj = 0; jj < nn; ++jj) {
               if (grid[0][j].s == index[0][jj] && grid[i][j].s == index[1][jj])
                    counts[jj]++;
          }

          if (jj == nn) {
               counts[jj]++;
               index[0][jj] = grid[0][j].s;
               index[1][jj] = grid[i][j].s;
               nn++;
          }
     }

     max_c =  0;
     j_max = -1;
     for (j = 0; j < n_derivs; ++j) {
          if (counts[j] > max_c) {
               max_c = counts[j];
               j_max = j;
          }
     }

     if (max_c != 0) {
          *s1 = index[0][j_max];
          *s2 = index[1][j_max];
     }

     free_array1_i(counts);
     free_array2(void **) index);

     return max_c;
}
*/


/*******************************************************************************
 *
 ******************************************************************************/
static void combine_grid_layers(grid_data **grid, stack_data *stacks,
                    int n_layers, int n_derivs, stack_data *s1, stack_data *s2) {

     int i;
     int jj;

     s1->count++;
     s2->count++;

     stacks[0].i1    = s1->i1;
     stacks[0].i2    = s2->i2;
     stacks[0].stack = 0;
     stacks[0].count = 0;
     stacks[0].p1    = s1;
     stacks[0].p2    = s2;

     for (jj = 0; jj < n_derivs; ++jj) {
          if (grid[0][jj].active) {
               i = next_layer(grid, 1, jj, n_layers);
               if (i >= 0) {
                    if ((! grid[0][jj].varies && grid[0][jj].s == s1) &&
                        (! grid[i][jj].varies && grid[i][jj].s == s2)) {
                         grid[0][jj].s = &stacks[0];

                         grid[i][jj].active = 0;
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static int check_for_duplicates(stack_data *stacks, int n_stacks, int n_layers) {

     int i;

     int count = 0;

     uchar **grid;

     grid = alloc_array2_uc(n_layers, n_layers);

     init_array2_uc(grid, n_layers, n_layers, 0);

     for (i = 0; i < n_stacks; ++i) {

          grid[stacks[i].i1][stacks[i].i2]++;

          if (grid[stacks[i].i1][stacks[i].i2] > 1) {
               printf("Duplicate found (%2d, %2d), count = %d\n",
                      stacks[i].i1, stacks[i].i2, grid[stacks[i].i1][stacks[i].i2]);
               count++;
          }
     }

     free_array2_uc(grid);

     return count;
}



/*******************************************************************************
 *
 ******************************************************************************/
int build_stack_chain(int n_layers, int n_derivs, uchar **derivs, stack_data *stacks) {

     int i;
     int j;

     int flag;

     int i_stack;

     int n_stacks;

     int count;
     int max_c;
     int max_i;

     stack_data *s1;
     stack_data *s2;

     stack_data *max_s1;
     stack_data *max_s2;

     grid_data **grid;

     grid = (grid_data **) alloc_array2(n_layers, n_derivs, sizeof(grid_data));

     for (i = 0; i < n_layers; ++i) {
          stacks[i].i1    = i;
          stacks[i].i2    = i;

          stacks[i].stack = 0;

          stacks[i].count = 0;

          stacks[i].p1    = NULL;
          stacks[i].p2    = NULL;
     }

     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_derivs; ++j) {
               flag = derivs[i][j];

               if (flag) {
                    grid[i][j].active = 0;
                    grid[i][j].varies = 1;
               }
               else {
                    grid[i][j].active = 1;
                    grid[i][j].varies = 0;
                    grid[i][j].s = &stacks[i];
               }
          }
     }
/*
     print_grid(grid, n_layers, n_derivs);
     printf("\n");
     print_stacks(stacks, n_layers);
     printf("\n");
*/
     i_stack = n_layers;

     while(1) {
          max_c = 0;
          for (i = 0; i < n_layers - 1; ++i) {

               count = max_set_of_pairs1(&grid[i], n_layers - i,
                                         n_derivs, &s1, &s2);
/*
               count = max_set_of_pairs2(&grid[i], n_layers - i,
                                         n_derivs, &s1, &s2);
*/
               if (count > max_c) {
                    max_c  = count;
                    max_i  = i;
                    max_s1 = s1;
                    max_s2 = s2;
               }
          }

          if (max_c == 0)
               break;

          combine_grid_layers(&grid[max_i], &stacks[i_stack],
                              n_layers - max_i, n_derivs, max_s1, max_s2);

          i_stack++;
/*
          printf("%04d\n", i_stack);
          print_grid(grid, n_layers, n_derivs);
          printf("\n");
          print_stacks(stacks, i_stack);
          printf("\n");
*/
     }

     n_stacks = i_stack;

     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_derivs; ++j) {
               if (grid[i][j].active)
                    grid[i][j].s->stack = 1;
          }
     }

     free_array2((void **) grid);

     check_for_duplicates(stacks, n_stacks, n_layers);
/*
     print_stacks(stacks, n_stacks);
     printf("\n");
*/
     return n_stacks;
}



/*******************************************************************************
 *
 ******************************************************************************/
int stack_chain_alloc(int n_four, int n_quad, int n_derivs, int n_layers, int n_stacks0, stack_data *stacks0, int n_stacks, stack_data *stacks, uchar **derivs_layers, uchar **derivs_s, uchar **derivs_d) {

     uchar *derivs_s2;

     int i;
     int ii;
     int j;
     int k;

     int r;

     jmp_buf env;

     if (! (derivs_s2 = flags_alloc(n_derivs))) {
          fprintf(stderr, "ERROR: flags_alloc()");
          return -1;
     }

     if ((r = setjmp(env)) == 0) {

     for (i = 0; i < n_stacks; ++i) {
          stacks[i].allocated = 1;

          if (stacks0 && (ii = stack_chain_find(n_stacks0, stacks0, stacks[i].i1, stacks[i].i2)) >= 0) {
               for (j = 0; j < n_four; ++j)
                    stacks[i].new_[j] = 0;

               stacks[i].R_p = stacks0[ii].R_p;
               stacks[i].T_p = stacks0[ii].T_p;
               stacks[i].R_m = stacks0[ii].R_m;
               stacks[i].T_m = stacks0[ii].T_m;

               stacks[i].S_p = stacks0[ii].S_p;
               stacks[i].S_m = stacks0[ii].S_m;

               stacks[i].Sl_p = stacks0[ii].Sl_p;
               stacks[i].Sl_m = stacks0[ii].Sl_m;

               if (i < n_layers) {
                    stacks[i].R_p_l = stacks0[ii].R_p_l;
                    stacks[i].T_p_l = stacks0[ii].T_p_l;
                    stacks[i].R_m_l = stacks0[ii].R_m_l;
                    stacks[i].T_m_l = stacks0[ii].T_m_l;

                    stacks[i].S_p_l = stacks0[ii].S_p_l;
                    stacks[i].S_m_l = stacks0[ii].S_m_l;

                    stacks[i].Sl_p_l = stacks0[ii].Sl_p_l;
                    stacks[i].Sl_m_l = stacks0[ii].Sl_m_l;
               }
               else {
                    stacks[i].S_p_l = stacks0[ii].S_p_l;
                    stacks[i].S_m_l = stacks0[ii].S_m_l;

                    stacks[i].Sl_p_l = stacks0[ii].Sl_p_l;
                    stacks[i].Sl_m_l = stacks0[ii].Sl_m_l;
               }

               stacks[i].d = stacks0[ii].d;

               stacks0[ii].allocated = 0;
          }
          else {
               for (j = 0; j < n_four; ++j)
                    stacks[i].new_[j] = 1;

               stacks[i].R_p = ALLOC_ARRAY3_D(n_four, n_quad, n_quad);
               stacks[i].T_p = ALLOC_ARRAY3_D(n_four, n_quad, n_quad);
               stacks[i].R_m = ALLOC_ARRAY3_D(n_four, n_quad, n_quad);
               stacks[i].T_m = ALLOC_ARRAY3_D(n_four, n_quad, n_quad);

               stacks[i].S_p = ALLOC_ARRAY2_D(n_four, n_quad);
               stacks[i].S_m = ALLOC_ARRAY2_D(n_four, n_quad);

               stacks[i].Sl_p = ALLOC_ARRAY2_D(n_four, n_quad);
               stacks[i].Sl_m = ALLOC_ARRAY2_D(n_four, n_quad);

               if (i < n_layers) {
/*
                    if (flags_or(derivs_layers[i], n_derivs)) {
*/
                         stacks[i].R_p_l = ALLOC_ARRAY2(n_four, n_derivs, double **);
                         stacks[i].T_p_l = ALLOC_ARRAY2(n_four, n_derivs, double **);
                         stacks[i].R_m_l = ALLOC_ARRAY2(n_four, n_derivs, double **);
                         stacks[i].T_m_l = ALLOC_ARRAY2(n_four, n_derivs, double **);

                         for (j = 0; j < n_four; ++j) {
                              for (k = 0; k < n_derivs; ++k) {
/*
if (derivs_layers[i][k]) {
*/
                                   stacks[i].R_p_l[j][k] = ALLOC_ARRAY2_D(n_quad, n_quad);
                                   stacks[i].T_p_l[j][k] = ALLOC_ARRAY2_D(n_quad, n_quad);
                                   stacks[i].R_m_l[j][k] = ALLOC_ARRAY2_D(n_quad, n_quad);
                                   stacks[i].T_m_l[j][k] = ALLOC_ARRAY2_D(n_quad, n_quad);
/*
}
*/
                              }
                         }
/*
                    }

                    if (flags_or(derivs_d[i], n_derivs)) {
*/
                         stacks[i].S_p_l = ALLOC_ARRAY2(n_four, n_derivs, double *);
                         stacks[i].S_m_l = ALLOC_ARRAY2(n_four, n_derivs, double *);

                         stacks[i].Sl_p_l = ALLOC_ARRAY2(n_four, n_derivs, double *);
                         stacks[i].Sl_m_l = ALLOC_ARRAY2(n_four, n_derivs, double *);

                         for (j = 0; j < n_four; ++j) {
                              for (k = 0; k < n_derivs; ++k) {
/*
if (derivs_d[i][k]) {
*/
                                   stacks[i].S_p_l[j][k] = ALLOC_ARRAY1_D(n_quad);
                                   stacks[i].S_m_l[j][k] = ALLOC_ARRAY1_D(n_quad);
/*
}
*/
                                   stacks[i].Sl_p_l[j][k] = ALLOC_ARRAY1_D(n_quad);
                                   stacks[i].Sl_m_l[j][k] = ALLOC_ARRAY1_D(n_quad);
                              }
                         }
/*
                    }
*/
               }
               else {
                    derivs_merge_s(n_derivs, derivs_s[stacks[i].i1], derivs_s[stacks[i].i2], derivs_s2);
/*
                    if (flags_or(derivs_s2, n_derivs)) {
*/
                         stacks[i].S_p_l = ALLOC_ARRAY2(n_four, n_derivs, double *);
                         stacks[i].S_m_l = ALLOC_ARRAY2(n_four, n_derivs, double *);

                         stacks[i].Sl_p_l = ALLOC_ARRAY2(n_four, n_derivs, double *);
                         stacks[i].Sl_m_l = ALLOC_ARRAY2(n_four, n_derivs, double *);

                         for (j = 0; j < n_four; ++j) {
                              for (k = 0; k < n_derivs; ++k) {
/*
if (derivs_s2[k]) {
*/
                                   stacks[i].S_p_l[j][k] = ALLOC_ARRAY1_D(n_quad);
                                   stacks[i].S_m_l[j][k] = ALLOC_ARRAY1_D(n_quad);
/*
}
*/
                                   stacks[i].Sl_p_l[j][k] = ALLOC_ARRAY1_D(n_quad);
                                   stacks[i].Sl_m_l[j][k] = ALLOC_ARRAY1_D(n_quad);
                              }
                         }
/*
                    }
*/
                    stacks[i].d = ALLOC_ARRAY1(n_four, layer_add_aux_data);

                    for (j = 0; j < n_four; ++j) {
                         if (layer_add_aux_alloc(&stacks[i].d[j], n_quad)) {
                              fprintf(stderr, "ERROR: layer_add_aux_alloc()");
                              return -1;
                         }
                    }
               }
          }
     }

     }
     else if (r < 0)
          return -1;

     flags_free(derivs_s2);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int stack_chain_free(int n_four, int n_quad, int n_derivs, int n_layers, int n_stacks, stack_data *stacks, uchar **derivs_layers, uchar **derivs_s, uchar **derivs_d) {

     uchar *derivs_s2;

     int i;
     int j;
     int k;

     if (! stacks)
          return 0;

     if (! (derivs_s2 = flags_alloc(n_derivs))) {
          fprintf(stderr, "ERROR: flags_alloc()");
          return -1;
     }

     for (i = 0; i < n_stacks; ++i) {
          if (stacks[i].allocated) {
               free_array3_d(stacks[i].R_p);
               free_array3_d(stacks[i].T_p);
               free_array3_d(stacks[i].R_m);
               free_array3_d(stacks[i].T_m);

               free_array2_d(stacks[i].S_p);
               free_array2_d(stacks[i].S_m);

               free_array2_d(stacks[i].Sl_p);
               free_array2_d(stacks[i].Sl_m);

               if (i < n_layers) {
/*
                    if (flags_or(derivs_layers[i], n_derivs)) {
*/
                         for (j = 0; j < n_four; ++j) {
                              for (k = 0; k < n_derivs; ++k) {
/*
if (derivs_layers[i][k]) {
*/
                                   free_array2_d(stacks[i].R_p_l[j][k]);
                                   free_array2_d(stacks[i].T_p_l[j][k]);
                                   free_array2_d(stacks[i].R_m_l[j][k]);
                                   free_array2_d(stacks[i].T_m_l[j][k]);
/*
}
*/
                              }
                         }

                         free_array2((void **) stacks[i].R_p_l);
                         free_array2((void **) stacks[i].T_p_l);
                         free_array2((void **) stacks[i].R_m_l);
                         free_array2((void **) stacks[i].T_m_l);
/*
                    }

                    if (flags_or(derivs_d[i], n_derivs)) {
*/
                         for (j = 0; j < n_four; ++j) {
                              for (k = 0; k < n_derivs; ++k) {
/*
if (derivs_d[i][k]) {
*/
                                   free_array1_d(stacks[i].S_p_l[j][k]);
                                   free_array1_d(stacks[i].S_m_l[j][k]);
/*
}
*/
                                   free_array1_d(stacks[i].Sl_p_l[j][k]);
                                   free_array1_d(stacks[i].Sl_m_l[j][k]);
                              }
                         }

                         free_array2((void **) stacks[i].S_p_l);
                         free_array2((void **) stacks[i].S_m_l);

                         free_array2((void **) stacks[i].Sl_p_l);
                         free_array2((void **) stacks[i].Sl_m_l);
/*
                    }
*/
               }
               else {
                    derivs_merge_s(n_derivs, derivs_s[stacks[i].i1], derivs_s[stacks[i].i2], derivs_s2);
/*
                    if (flags_or(derivs_s2, n_derivs)) {
*/
                         for (j = 0; j < n_four; ++j) {
                              for (k = 0; k < n_derivs; ++k) {
/*
if (derivs_s2[k]) {
*/
                                   free_array1_d(stacks[i].S_p_l[j][k]);
                                   free_array1_d(stacks[i].S_m_l[j][k]);
/*
}
*/
                                   free_array1_d(stacks[i].Sl_p_l[j][k]);
                                   free_array1_d(stacks[i].Sl_m_l[j][k]);
                              }
                         }

                         free_array2((void **) stacks[i].S_p_l);
                         free_array2((void **) stacks[i].S_m_l);

                         free_array2((void **) stacks[i].Sl_p_l);
                         free_array2((void **) stacks[i].Sl_m_l);
/*
                    }
*/
                    for (j = 0; j < n_four; ++j)
                         layer_add_aux_free(&stacks[i].d[j]);

                    free_array1(stacks[i].d);
               }
          }
     }

     flags_free(derivs_s2);

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
int stack_chain_find(int n_stacks, stack_data *stacks, int i1, int i2) {

    int i;

    for (i = 0; i < n_stacks; ++i) {
         if (stacks[i].i1 == i1 && stacks[i].i2 == i2) {
              return i;
         }
    }

    return -1;
}



int stack_chain_find2(int n_stacks, stack_data *stacks, int i11, int i12, int i21, int i22) {

    int i;

    for (i = 0; i < n_stacks; ++i) {
         if (stacks[i].p1 && stacks[i].p2) {
              if (stacks[i].p1->i1 == i11 && stacks[i].p1->i2 == i12 &&
                  stacks[i].p2->i1 == i21 && stacks[i].p2->i2 == i22) {
                   return i;
              }
         }
    }

    return -1;
}
