/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


/*******************************************************************************
 *
 ******************************************************************************/
static void CALC_GLOBAL_R_AND_T(int n_quad, int n_derivs,
                                double ltau, double *ltau_l,
                                TYPE  *nu, TYPE  **X_p, TYPE  **X_m,
                                double  **R, double  **T,
                                TYPE **nu_l, TYPE ***X_p_l, TYPE ***X_m_l,
                                double ***R_l, double ***T_l,
                                int symmetric, uchar *derivs,
                                save_tree_data save_tree, work_data work) {

     int i;
     int j;

     int *i1;
     int *i2;

     TYPE *lambda;

     TYPE **w1;

     TYPE **X_m_lu;

     TYPE **a;
     TYPE **b;
     TYPE **c_lu;
     TYPE **d;
     TYPE **e;

     FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, SAVE_TREE_STRING);

          if (save_tree_retrieve_data(&save_tree, FORWARD_SAVE_CALC_GLOBAL_R_AND_T_DATA, &save))
               FORWARD_SAVE_CALC_GLOBAL_R_AND_T_ALLOC(save, n_quad);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1     = get_work1(&work, WORK_IX);
     i2     = get_work1(&work, WORK_IX);

     lambda = get_work1(&work, WORK_XX);

     w1     = get_work1(&work, WORK_XXX);

     X_m_lu = get_work1(&work, WORK_XXX);

     a      = get_work1(&work, WORK_XXX);
     b      = get_work1(&work, WORK_XXX);
     c_lu   = get_work1(&work, WORK_XXX);
     d      = get_work1(&work, WORK_XXX);
     e      = get_work1(&work, WORK_XXX);


     /*-------------------------------------------------------------------------
      * lambda
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_quad; ++i)
          lambda[i] = XEXP(-nu[i] * ltau);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     xmat_mul_diag(X_p, lambda, a, n_quad, n_quad);


     xmat_trans(X_m, X_m_lu, n_quad, n_quad);
     xmat_getrf(X_m_lu, n_quad, n_quad, i1);

     xmat_copy(b, a, n_quad, n_quad);
     xmat_getrs2('n', X_m_lu, b, n_quad, n_quad, i1);


     xmat_mul(b, a, n_quad, n_quad, n_quad, w1);
     xmat_sub(X_m, w1, w1, n_quad, n_quad);

     xmat_trans(w1, c_lu, n_quad, n_quad);
     xmat_getrf(c_lu, n_quad, n_quad, i2);


     xmat_copy(d, X_p, n_quad, n_quad);
     xmat_getrs2('n', c_lu, d, n_quad, n_quad, i2);


     xmat_mul_diag(X_m, lambda, e, n_quad, n_quad);
     xmat_getrs2('n', c_lu, e, n_quad, n_quad, i2);
#ifdef REAL
if (! symmetric) {
     xmat_mul(e, b, n_quad, n_quad, n_quad, w1);
     xmat_sub(w1, d, R, n_quad, n_quad);
}
else {
     xsym_mul(e, b, n_quad, n_quad, w1, 0);
     xsym_sub(w1, d, R, n_quad, 1);
}
if (! symmetric) {
     xmat_mul(d, b, n_quad, n_quad, n_quad, w1);
     xmat_sub(e, w1, T, n_quad, n_quad);
}
else {
     xsym_mul(d, b, n_quad, n_quad, w1, 0);
     xsym_sub(e, w1, T, n_quad, 1);
}
#else
     xmat_mul(e, b, n_quad, n_quad, n_quad, w1);
     xmat_sub(w1, d, w1, n_quad, n_quad);

     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               R[i][j] = XREAL(w1[i][j]);
          }
     }

     xmat_mul(d, b, n_quad, n_quad, n_quad, w1);
     xmat_sub(e, w1, w1, n_quad, n_quad);

     for (i = 0; i < n_quad; ++i) {
          for (j = 0; j < n_quad; ++j) {
               T[i][j] = XREAL(w1[i][j]);
          }
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          copy_array1_i(save->X_m_ip, i1, n_quad);
          copy_array1_i(save->c_ip,   i2, n_quad);

          xmat_copy((TYPE **) save->X_m_lu, X_m_lu, n_quad, n_quad);

          xvec_copy((TYPE *)  save->lambda, lambda, n_quad);
          xmat_copy((TYPE **) save->a,      a,      n_quad, n_quad);
          xmat_copy((TYPE **) save->b,      b,      n_quad, n_quad);
          xmat_copy((TYPE **) save->c_lu,   c_lu,   n_quad, n_quad);
          xmat_copy((TYPE **) save->d,      d,      n_quad, n_quad);
          xmat_copy((TYPE **) save->e,      e,      n_quad, n_quad);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef USE_AD_FOR_TL_CALC_CALC_GLOBAL_R_AND_T
     if (flags_or(derivs, n_derivs)) {
#ifndef REAL
          int k;

          TYPE **w2;
#endif
          TYPE *lambda_l;

          TYPE **f;
          TYPE **g;
          TYPE **h;
          TYPE **p;
          TYPE **q;

          lambda_l = get_work1(&work, WORK_XX);
#ifndef REAL
          w2       = get_work1(&work, WORK_XXX);
#endif
          f        = get_work1(&work, WORK_XXX);
          g        = get_work1(&work, WORK_XXX);
          h        = get_work1(&work, WORK_XXX);
          p        = get_work1(&work, WORK_XXX);
          q        = get_work1(&work, WORK_XXX);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (i = 0; i < n_derivs; ++i) {

               if (! derivs[i])
                    continue;


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               for (j = 0; j < n_quad; ++j)
                    lambda_l[j] = -(nu_l[i][j] * ltau + nu[j] * ltau_l[i]) * lambda[j];


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               xmat_mul_diag(X_p_l[i], lambda, f, n_quad, n_quad);

               xmat_mul_diag(X_p, lambda_l, w1, n_quad, n_quad);
               xmat_add(f, w1, f, n_quad, n_quad);
/*
               xmat_gxdxmx(0, X_p, lambda_l, 1., f, 1., n_quad, n_quad);
*/
               xmat_mul(b, X_m_l[i], n_quad, n_quad, n_quad, w1);
               xmat_sub(f, w1, g, n_quad, n_quad);
               xmat_getrs2('n', X_m_lu, g, n_quad, n_quad, i1);

               xmat_mul(g, a, n_quad, n_quad, n_quad, w1);
               xmat_sub(X_m_l[i], w1, h, n_quad, n_quad);

               xmat_mul(b, f, n_quad, n_quad, n_quad, w1);
               xmat_sub(h, w1, h, n_quad, n_quad);
/*
               xmat_gxgxmx(0, b, 0, f, -1., h, 1., n_quad, n_quad, n_quad);
*/
               xmat_mul(d, h, n_quad, n_quad, n_quad, w1);
               xmat_sub(X_p_l[i], w1, p, n_quad, n_quad);
               xmat_getrs2('n', c_lu, p, n_quad, n_quad, i2);

               xmat_mul_diag(X_m_l[i], lambda, q, n_quad, n_quad);
               xmat_mul_diag(X_m, lambda_l, w1, n_quad, n_quad);
               xmat_add(q, w1, q, n_quad, n_quad);
/*
               xmat_gxdxmx(0, X_m_l[i], lambda, 1., q, 0., n_quad, n_quad);
               xmat_gxdxmx(0, X_m, lambda_l, 1., q, 1., n_quad, n_quad);
*/
               xmat_mul(e, h, n_quad, n_quad, n_quad, w1);
               xmat_sub(q, w1, q, n_quad, n_quad);
/*
               xmat_gxgxmx(0, e, 0, h, -1., q, 1., n_quad, n_quad, n_quad);
*/
               xmat_getrs2('n', c_lu, q, n_quad, n_quad, i2);
#ifdef REAL
if (! symmetric) {
               xmat_mul(q, b, n_quad, n_quad, n_quad, R_l[i]);

               xmat_mul(e, g, n_quad, n_quad, n_quad, w1);
               xmat_add(R_l[i], w1, R_l[i], n_quad, n_quad);
/*
               xmat_gxgxmx(0, e, 0, g, 1., R_l[i], 1., n_quad, n_quad, n_quad);
*/
               xmat_sub(R_l[i], p, R_l[i], n_quad, n_quad);
}
else {
               xsym_mul(q, b, n_quad, n_quad, R_l[i], 0);
               xsym_mul(e, g, n_quad, n_quad, w1, 0);
               xsym_add(R_l[i], w1, R_l[i], n_quad, 0);
               xsym_sub(R_l[i], p, R_l[i], n_quad, 1);
}
if (! symmetric) {
               xmat_mul(p, b, n_quad, n_quad, n_quad, w1);
               xmat_sub(q, w1, T_l[i], n_quad, n_quad);

               xmat_mul(d, g, n_quad, n_quad, n_quad, w1);
               xmat_sub(T_l[i], w1, T_l[i], n_quad, n_quad);
/*
               xmat_gxgxmx(0, d, 0, g, -1., T_l[i], 1., n_quad, n_quad, n_quad);
*/
}
else {
               xsym_mul(p, b, n_quad, n_quad, w1, 0);
               xsym_sub(q, w1, T_l[i], n_quad, 0);
               xsym_mul(d, g, n_quad, n_quad, w1, 0);
               xsym_sub(T_l[i], w1, T_l[i], n_quad, 1);
}
#else
               xmat_mul(q, b, n_quad, n_quad, n_quad, w1);

               xmat_mul(e, g, n_quad, n_quad, n_quad, w2);
               xmat_add(w1, w2, w1, n_quad, n_quad);
/*
               xmat_gxgxmx(0, e, 0, g, 1., w1, 1., n_quad, n_quad, n_quad);
*/
               xmat_sub(w1, p, w1, n_quad, n_quad);

               for (j = 0; j < n_quad; ++j) {
                    for (k = 0; k < n_quad; ++k) {
                         R_l[i][j][k] = XREAL(w1[j][k]);
                    }
               }

               xmat_mul(p, b, n_quad, n_quad, n_quad, w1);
               xmat_sub(q, w1, w1, n_quad, n_quad);

               xmat_mul(d, g, n_quad, n_quad, n_quad, w2);
               xmat_sub(w1, w2, w1, n_quad, n_quad);
/*
               xmat_gxgxmx(0, d, 0, g, -1., w1, 1., n_quad, n_quad, n_quad);
*/
               for (j = 0; j < n_quad; ++j) {
                    for (k = 0; k < n_quad; ++k) {
                         T_l[i][j][k] = XREAL(w1[j][k]);
                    }
               }
#endif
          }
     }
#else
     CALC_GLOBAL_R_AND_T_TL_WITH_AD(n_quad, n_derivs, ltau, ltau_l, nu, X_p, X_m,
                                    R, T, nu_l, X_p_l, X_m_l, R_l, T_l, symmetric, derivs, save_tree, work);
#endif
}
