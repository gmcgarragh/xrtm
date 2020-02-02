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
void SOLVE_BVP_TL_WITH_AD(int n_quad, int n_stokes, int n_derivs, int n_layers,
            double *ltau, double **ltau_l,
            double **Rs_qq, double ***Rs_qq_l,
            double *atran, double **atran_l,
            TYPE **nu, TYPE ***X_p, TYPE ***X_m,
            double **F_p, double **F_m,
            double **F0_p, double **F0_m, double **F1_p, double **F1_m,
            TYPE ***nu_l, TYPE ****X_p_l, TYPE ****X_m_l,
            double ***F_p_l, double ***F_m_l,
            double ***F0_p_l, double ***F0_m_l, double ***F1_p_l, double ***F1_m_l,
            TYPE *B, TYPE **B_l,
            double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
            int surface, int thermal, uchar **derivs_layers, uchar **derivs_beam_down,
            save_tree_data save_tree, work_data work) {

     uchar *derivs_layers_union;
     uchar *derivs_beam_down_union;

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;

     int n_quad_v;
     int n_comp;

     int i_input;
     int n_inputs;

     int ii_input;

     double *ltau_a;

     double **Rs_qq_a;

     double *atran_a;

     double **F_p_a;
     double **F_m_a;

     double *I1_m_a;
     double *In_p_a;

     TYPE **a;

     TYPE *B_a;

     TYPE **nu_a;
     TYPE ***X_p_a;
     TYPE ***X_m_a;


     if (! flags_or2(derivs_beam_down, n_layers, n_derivs))
          return;


     n_quad_v  = n_quad * n_stokes;

     n_comp = 2 * n_quad_v * n_layers;


     save_tree_decode_s(&save_tree, SOLVE_BVP_SAVE_TREE_STRING);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     derivs_layers_union = get_work_uc1(&work, n_layers + 1);
     derivs_beam_down_union = get_work_uc1(&work, n_layers + 1);

     derivs_union_logical_or2(n_layers + 1, n_derivs, derivs_layers, derivs_layers_union);
     derivs_union_bitwise_or2(n_layers + 1, n_derivs, derivs_beam_down, derivs_beam_down_union);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     B_a     = get_work_x1(&work, n_comp);

     ltau_a  = get_work_d1(&work, n_layers);
     Rs_qq_a = get_work_d2(&work, n_quad_v, n_quad_v);
     atran_a = get_work_d1(&work, n_layers);
     nu_a    = get_work_x2(&work, n_layers, n_quad_v);
     X_p_a   = get_work_x3(&work, n_layers, n_quad_v, n_quad_v);
     X_m_a   = get_work_x3(&work, n_layers, n_quad_v, n_quad_v);
     F_p_a   = get_work_d2(&work, n_layers, n_quad_v);
     F_m_a   = get_work_d2(&work, n_layers, n_quad_v);
     I1_m_a  = get_work_d1(&work, n_quad_v);
     In_p_a  = get_work_d1(&work, n_quad_v);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     XCAT(init_array1_, TYPE_POSTFIX)(B_a, n_comp, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = n_comp;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_comp) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array1_d(ltau_a,  n_layers, 0.);
          init_array2_d(Rs_qq_a, n_quad_v, n_quad_v, 0.);
          init_array1_d(atran_a, n_layers, 0.);
          init_array2_d(F_p_a,   n_layers, n_quad_v, 0.);
          init_array2_d(F_m_a,   n_layers, n_quad_v, 0.);
          init_array1_d(I1_m_a,  n_quad_v, 0.);
          init_array1_d(In_p_a,  n_quad_v, 0.);

          XCAT(init_array2_, TYPE_POSTFIX)(nu_a,  n_layers, n_quad_v, 0.);
          XCAT(init_array3_, TYPE_POSTFIX)(X_p_a, n_layers, n_quad_v, n_quad_v, 0.);
          XCAT(init_array3_, TYPE_POSTFIX)(X_m_a, n_layers, n_quad_v, n_quad_v, 0.);


          if (i_input == 0)
               B_a[ii_input] = 1.;


          SOLVE_BVP_A(n_quad, n_stokes, n_layers, ltau, ltau_a, Rs_qq, Rs_qq_a, atran, atran_a, nu, X_p, X_m, F_p, F_m, F0_p, F0_m, F1_p, F1_m, nu_a, X_p_a, X_m_a, F_p_a, F_m_a, NULL, NULL, NULL, NULL, B, B_a, I1_m, I1_m_a, In_p, In_p_a, surface, thermal, derivs_layers_union, derivs_beam_down_union, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = B_l;

               for (k = 0; k < n_derivs; ++k)
                    a[k][ii_input] = 0.;

               for (j = 0; j < n_layers; ++j) {
                    for (k = 0; k < n_derivs; ++k) {
                         a[k][ii_input] += ltau_l [j][k] * ltau_a [j];
                         a[k][ii_input] += atran_l[j][k] * atran_a[j];

                         for (l = 0; l < n_quad_v; ++l) {
                              if (derivs_layers[j][k]) {
                                   a[k][ii_input] += nu_l[j][k][l] * nu_a[j][l];

                                   for (m = 0; m < n_quad_v; ++m) {
                                        a[k][ii_input] += X_p_l[j][k][l][m] * X_p_a[j][l][m];
                                        a[k][ii_input] += X_m_l[j][k][l][m] * X_m_a[j][l][m];
                                   }
                              }

                              if (derivs_beam_down[j][k]) {
                                   a[k][ii_input] += F_p_l[j][k][l] * F_p_a[j][l];
                                   a[k][ii_input] += F_m_l[j][k][l] * F_m_a[j][l];
                              }
                         }
                    }
               }

               for (k = 0; k < n_derivs; ++k) {
                    for (l = 0; l < n_quad_v; ++l) {
                         a[k][ii_input] += I1_m_l[k][l] * I1_m_a[l];
                         a[k][ii_input] += In_p_l[k][l] * In_p_a[l];

                         for (m = 0; m < n_quad_v; ++m) {
                              if (surface && derivs_layers[j][k])
                                   a[k][ii_input] += Rs_qq_l[k][l][m] * Rs_qq_a[l][m];
                         }
                    }
               }
          }
     }
}
