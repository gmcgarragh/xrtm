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
void CALC_GLOBAL_R_AND_T_TL_WITH_AD(int n_quad, int n_derivs,
                                    double ltau, double *ltau_l,
                                    TYPE  *nu, TYPE  **X_p, TYPE  **X_m,
                                    double  **R, double  **T,
                                    TYPE **nu_l, TYPE ***X_p_l, TYPE ***X_m_l,
                                    double ***R_l, double ***T_l,
                                    int symmetric, uchar *derivs,
                                    save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double ltau_a;

     double **R_a;
     double **T_a;

     double ***a;

     TYPE *nu_a;
     TYPE **X_p_a;
     TYPE **X_m_a;


     if (! flags_or(derivs, n_derivs))
          return;


     save_tree_decode_s(&save_tree, SAVE_TREE_STRING);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R_a   = get_work1(&work, WORK_DXX);
     T_a   = get_work1(&work, WORK_DXX);

     nu_a  = get_work1(&work, WORK_XX);
     X_p_a = get_work1(&work, WORK_XXX);
     X_m_a = get_work1(&work, WORK_XXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R_a, n_quad, n_quad);
     dmat_zero(T_a, n_quad, n_quad);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * n_quad * n_quad;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad * n_quad) {
               i_input = 0;
               ii = i;
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }
          else {
               i_input = 1;
               ii = i - n_quad * n_quad;
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          ltau_a = 0.;

          xvec_zero(nu_a, n_quad);
          xmat_zero(X_p_a, n_quad, n_quad);
          xmat_zero(X_m_a, n_quad, n_quad);


          if (i_input == 0)
               R_a[ii_input][jj_input] = 1.;
          else
               T_a[ii_input][jj_input] = 1.;


          CALC_GLOBAL_R_AND_T_A(n_quad, ltau, &ltau_a, nu, X_p, X_m, R, T,
                                nu_a, X_p_a, X_m_a, R_a, T_a, 0, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = R_l;
               else
                    a = T_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    a[j][ii_input][jj_input] = 0.;

                    a[j][ii_input][jj_input] += ltau_l[j] * ltau_a;

                    for (k = 0; k < n_quad; ++k) {
                         a[j][ii_input][jj_input] += XREAL(nu_l[j][k] * nu_a[k]);

                         for (l = 0; l < n_quad; ++l) {
                              a[j][ii_input][jj_input] += XREAL(X_p_l[j][k][l] * X_p_a[k][l]);
                              a[j][ii_input][jj_input] += XREAL(X_m_l[j][k][l] * X_m_a[k][l]);
                         }
                    }
               }
          }
     }
}
