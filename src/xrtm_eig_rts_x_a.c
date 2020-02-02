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
void FORWARD_SAVE_CALC_GLOBAL_R_AND_T_FREE(FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA *d);

int FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC(FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA *d, int n_quad) {

     d->free = (void (*)(void *)) FORWARD_SAVE_CALC_GLOBAL_R_AND_T_FREE;

     d->X_m_ip = alloc_array1_i(n_quad);
     d->c_ip   = alloc_array1_i(n_quad);

     d->X_m_lu = XCAT(alloc_array2_, TYPE_POSTFIX)(n_quad, n_quad);

     d->lambda = XCAT(alloc_array1_, TYPE_POSTFIX)(n_quad);
     d->a      = XCAT(alloc_array2_, TYPE_POSTFIX)(n_quad, n_quad);
     d->b      = XCAT(alloc_array2_, TYPE_POSTFIX)(n_quad, n_quad);
     d->c_lu   = XCAT(alloc_array2_, TYPE_POSTFIX)(n_quad, n_quad);
     d->d      = XCAT(alloc_array2_, TYPE_POSTFIX)(n_quad, n_quad);
     d->e      = XCAT(alloc_array2_, TYPE_POSTFIX)(n_quad, n_quad);

     return 0;
}



void FORWARD_SAVE_CALC_GLOBAL_R_AND_T_FREE(FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA *d) {

     free_array1_i(d->X_m_ip);
     free_array1_i(d->c_ip);

     XCAT(free_array2_, TYPE_POSTFIX)(d->X_m_lu);

     XCAT(free_array1_, TYPE_POSTFIX)(d->lambda);
     XCAT(free_array2_, TYPE_POSTFIX)(d->a);
     XCAT(free_array2_, TYPE_POSTFIX)(d->b);
     XCAT(free_array2_, TYPE_POSTFIX)(d->c_lu);
     XCAT(free_array2_, TYPE_POSTFIX)(d->d);
     XCAT(free_array2_, TYPE_POSTFIX)(d->e);
}



/*******************************************************************************
 *
 ******************************************************************************/
void CALC_GLOBAL_R_AND_T_A(int n_quad,
                           double ltau, double *ltau_a,
                           TYPE  *nu, TYPE  **X_p, TYPE  **X_m,
                           double  **R, double  **T,
                           TYPE *nu_a, TYPE **X_p_a, TYPE **X_m_a,
                           double **R_a, double **T_a,
                           int symmetric, save_tree_data save_tree, work_data work) {

     int i;
#ifndef REAL
     int j;

     TYPE **w1;
#endif
     TYPE **T_a2;
     TYPE **R_a2;

     TYPE **t1_a;

     TYPE **lambda_a;

     TYPE **f_a;
     TYPE **g_a;
     TYPE **h_a;
     TYPE **p_a;
     TYPE **q_a;

     FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, SAVE_TREE_STRING);

     save_tree_retrieve_data(&save_tree, FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef REAL
     w1       = get_work1(&work, WORK_XXX);
#endif
     t1_a     = get_work1(&work, WORK_XXX);
     lambda_a = get_work1(&work, WORK_XXX);
     f_a      = get_work1(&work, WORK_XXX);
     g_a      = get_work1(&work, WORK_XXX);
     h_a      = get_work1(&work, WORK_XXX);
     p_a      = get_work1(&work, WORK_XXX);
     q_a      = get_work1(&work, WORK_XXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifdef REAL
     T_a2 = T_a;
#else
     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               w1[i][j] = T_a[i][j];
          }
     }

     T_a2 = w1;
#endif
     xmat_copy(q_a, T_a2, n_quad, n_quad);

     xmat_gxgxmx(0, T_a2, 1, save->b, -1., p_a, 0., n_quad, n_quad, n_quad);

     xmat_gxgxmx(1, save->d, 0, T_a2, -1., g_a, 0., n_quad, n_quad, n_quad);
#ifdef REAL
     R_a2 = R_a;
#else
     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               w1[i][j] = R_a[i][j];
          }
     }

     R_a2 = w1;
#endif
     xmat_gxgxmx(0, R_a2, 1, save->b, 1., q_a, 1., n_quad, n_quad, n_quad);

     xmat_gxgxmx(1, save->e, 0, R_a2, 1., g_a, 1., n_quad, n_quad, n_quad);

     xmat_sub(p_a, R_a2, p_a, n_quad, n_quad);

     xmat_copy(t1_a, q_a, n_quad, n_quad);
     xmat_getrs2('t', save->c_lu, t1_a, n_quad, n_quad, save->c_ip);

     xmat_gxdxmx(0, t1_a, save->lambda, 1., X_m_a, 1., n_quad, n_quad);

     xmat_gxgxmx(1, X_m, 0, t1_a, 1., lambda_a, 0., n_quad, n_quad, n_quad);

     xmat_gxgxmx(1, save->e, 0, t1_a, -1., h_a, 0., n_quad, n_quad, n_quad);

     xmat_copy(t1_a, p_a, n_quad, n_quad);
     xmat_getrs2('t', save->c_lu, t1_a, n_quad, n_quad, save->c_ip);

     xmat_add(X_p_a, t1_a, X_p_a, n_quad, n_quad);

     xmat_gxgxmx(1, save->d, 0, t1_a, -1., h_a, 1., n_quad, n_quad, n_quad);

     xmat_add(X_m_a, h_a, X_m_a, n_quad, n_quad);

     xmat_gxgxmx(0, h_a, 1, save->a, -1., g_a, 1., n_quad, n_quad, n_quad);

     xmat_gxgxmx(1, save->b, 0, h_a, -1., f_a, 0., n_quad, n_quad, n_quad);

     xmat_copy(t1_a, g_a, n_quad, n_quad);
     xmat_getrs2('t', save->X_m_lu, t1_a, n_quad, n_quad, save->X_m_ip);

     xmat_add(f_a, t1_a, f_a, n_quad, n_quad);

     xmat_gxgxmx(1, save->b, 0, t1_a, -1., X_m_a, 1., n_quad, n_quad, n_quad);

     xmat_gxdxmx(0, f_a, save->lambda, 1., X_p_a, 1., n_quad, n_quad);

     xmat_gxgxmx(1, X_p, 0, f_a, 1., lambda_a, 1., n_quad, n_quad, n_quad);

     xmat_mul_diag(lambda_a, save->lambda, t1_a, n_quad, n_quad);

     for (i = 0; i < n_quad; ++i) {
          nu_a[i] -= t1_a[i][i] * ltau;
          *ltau_a += -XREAL(nu[i] * t1_a[i][i]);
     }

     dmat_zero(R_a, n_quad, n_quad);
     dmat_zero(T_a, n_quad, n_quad);
}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_eig_rts_x_a2.c"
#endif
