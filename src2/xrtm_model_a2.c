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
void update_opt_props_tl_with_ad(xrtm_data *d, int i_layer, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int n_coef;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double *a;
     double ***b;


     if (! (d->options & XRTM_OPTION_CALC_DERIVS))
          return;

     n_coef = MIN(d->n_coef_layer[i_layer], d->n_coef2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     d->ltau_a [i_layer] = 0.;
     d->omega_a[i_layer] = 0.;
     init_array2_d(d->coef_a[i_layer], d->n_elem, n_coef, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 + d->n_elem * n_coef;

     for (i = 0; i < n_inputs; ++i) {
          if (i < 1)
               i_input = 0;
          else
          if (i < 2)
               i_input = 1;
          else {
               i_input = 2;
               ii = i - 2;
               ii_input = ii / n_coef;
               jj_input = ii % n_coef;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_opt_props0_a(d, i_layer, n_coef, work);


          if (i_input == 0)
               d->ltau_a[i_layer]  = 1.;
          else
          if (i_input == 1)
               d->omega_a[i_layer] = 1.;
          else
               d->coef_a[i_layer][ii_input][jj_input] = 1.;


          update_opt_props_a(d, i_layer, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = d->ltau_l [i_layer];
               else
                    a = d->omega_l[i_layer];

               for (j = 0; j < d->n_derivs; ++j) {
                    a[j]  = 0.;

                    a[j] += d->omega0_l[i_layer][j] * d->omega0_a[i_layer];
                    a[j] += d->ltau0_l [i_layer][j] * d->ltau0_a [i_layer];

                    a[j] += d->coef0_l[i_layer][j][0][d->n_coef2] * d->coef0_a[i_layer][0][d->n_coef2];

                    for (k = 0; k < d->n_elem; ++k) {
                         for (l = 0; l < n_coef; ++l) {
                              a[j] += d->coef0_l[i_layer][j][k][l] * d->coef0_a[i_layer][k][l];
                         }
                    }
               }
          }


          if (i_input == 2) {
               if (i_input == 2)
                    b = d->coef_l[i_layer];

               for (j = 0; j < d->n_derivs; ++j) {
                    b[j][ii_input][jj_input]  = 0.;

                    b[j][ii_input][jj_input] += d->omega0_l[i_layer][j] * d->omega0_a[i_layer];
                    b[j][ii_input][jj_input] += d->ltau0_l [i_layer][j] * d->ltau0_a [i_layer];

                    b[j][ii_input][jj_input] += d->coef0_l[i_layer][j][0][d->n_coef2] * d->coef0_a[i_layer][0][d->n_coef2];

                    for (k = 0; k < d->n_elem; ++k) {
                         for (l = 0; l < n_coef; ++l) {
                              b[j][ii_input][jj_input] += d->coef0_l[i_layer][j][k][l] * d->coef0_a[i_layer][k][l];
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void get_local_r_t_u_w_tl_with_ad
     (xrtm_data *d, int i_four, int i_layer,
      double **r_p, double **t_p, double **r_m, double **t_m,
      double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
      int flag, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     int n_coef;

     double **r_p_a;
     double **t_p_a;
     double **r_m_a;
     double **t_m_a;

     double ***a;


     if (! flags_or(d->derivs.layers[i_layer], d->n_derivs))
          return;


     save_tree_decode_s(&save_tree, "get_local_r_t_u_w");


     n_coef = MIN(d->n_coef_layer[i_layer], d->n_coef2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     r_p_a = get_work1(&work, WORK_DXX);
     t_p_a = get_work1(&work, WORK_DXX);

     if (d->options & XRTM_OPTION_VECTOR) {
          r_m_a = get_work1(&work, WORK_DXX);
          t_m_a = get_work1(&work, WORK_DXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(r_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(t_p_a, d->n_quad_v_x, d->n_quad_v_x);

     if (d->options & XRTM_OPTION_VECTOR) {
          dmat_zero(r_m_a, d->n_quad_v_x, d->n_quad_v_x);
          dmat_zero(t_m_a, d->n_quad_v_x, d->n_quad_v_x);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (! (d->options & XRTM_OPTION_VECTOR))
          n_inputs = 2 * d->n_quad_v_x * d->n_quad_v_x;
     else
          n_inputs = 4 * d->n_quad_v_x * d->n_quad_v_x;

     for (i = 0; i < n_inputs; ++i) {
          if (i < d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 0;
               ii = i;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 2 * d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 1;
               ii = i - d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 3 * d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 2;
               ii = i - 2 * d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else {
               i_input = 3;
               ii = i - 3 * d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          d->omega_a[i_layer] = 0.;

          for (j = 0; j < d->n_elem; ++j) {
               for (k = i_four; k < n_coef; ++k) {
                    d->coef_a[i_layer][j][k] = 0.;
               }
          }


          if (i_input == 0)
               r_p_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 1)
               t_p_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 2)
               r_m_a[ii_input][jj_input] = 1.;
          else
          t_m_a[ii_input][jj_input] = 1.;


          get_local_r_t_u_w_a(d, i_four, i_layer, r_p, t_p, r_m, t_m, r_p_a, t_p_a, r_m_a, t_m_a, flag, save_tree, &work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1 || i_input == 2 || i_input == 3) {
               if (i_input == 0)
                    a = r_p_l;
               else
               if (i_input == 1)
                    a = t_p_l;
               else
               if (i_input == 2)
                    a = r_m_l;
               else
                    a = t_m_l;

               for (j = 0; j < d->n_derivs; ++j) {
                    if (! d->derivs.layers[i_layer][j])
                         continue;

                    a[j][ii_input][jj_input]  = 0.;

                    a[j][ii_input][jj_input] += d->omega_l[i_layer][j] * d->omega_a[i_layer];

                    for (k = 0; k < d->n_elem; ++k) {
                         for (l = i_four; l < n_coef; ++l) {
                              a[j][ii_input][jj_input] += d->coef_l[i_layer][j][k][l] * d->coef_a[i_layer][k][l];
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void get_layer_R_T_S_U_W_V_tl_with_ad
     (xrtm_data *d, int i_four, int i_layer, int solver,
      double **R_p, double **T_p, double **R_m, double **T_m,
      double *S_p, double *S_m,
      double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
      double **S_p_l, double **S_m_l,
      int flag, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double **R_p_a;
     double **T_p_a;
     double **R_m_a;
     double **T_m_a;
     double *S_p_a;
     double *S_m_a;

     double **a;
     double ***b;


     if (! flags_or(d->derivs.sources[i_layer], d->n_derivs))
          return;


     save_tree_decode_s(&save_tree, "get_layer_R_T_S_U_W_V");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R_p_a = get_work1(&work, WORK_DXX);
     T_p_a = get_work1(&work, WORK_DXX);
     R_m_a = get_work1(&work, WORK_DXX);
     T_m_a = get_work1(&work, WORK_DXX);
     S_p_a = get_work1(&work, WORK_DX);
     S_m_a = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(R_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dvec_zero(S_p_a, d->n_quad_v_x);
     dvec_zero(S_m_a, d->n_quad_v_x);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 4 * d->n_quad_v_x * d->n_quad_v_x + 2 * d->n_quad_v_x;

     for (i = 0; i < n_inputs; ++i) {
          if (i < d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 0;
               ii = i;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 2 * d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 1;
               ii = i - d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 3 * d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 2;
               ii = i - 2 * d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 4 * d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 3;
               ii = i - 3 * d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 4 * d->n_quad_v_x * d->n_quad_v_x + d->n_quad_v_x) {
               i_input = 4;
               ii = i - 4 * d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii;
          }
          else {
               i_input = 5;
               ii = i - (4 * d->n_quad_v_x * d->n_quad_v_x + d->n_quad_v_x);
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_opt_props0_up_a(d, i_layer, d->n_coef, work);

          if (d->options & XRTM_OPTION_DELTA_M)
               init_opt_props_up_a(d, i_layer, d->n_coef, work);

          init_beam_params_up_a(d, i_layer, work);


          if (i_input == 0)
               R_p_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 1)
               T_p_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 2)
               R_m_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 3)
               T_m_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 4)
               S_p_a[ii_input] = 1.;
          else
               S_m_a[ii_input] = 1;


          get_layer_R_T_S_U_W_V_a(d, i_four, i_layer, solver, R_p, T_p, R_m, T_m, S_p, S_m, R_p_a, T_p_a, R_m_a, T_m_a, S_p_a, S_m_a, flag, save_tree, &work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 4 || i_input == 5) {
               if (i_input == 4)
                    a = S_p_l;
               else
                    a = S_m_l;

               for (j = 0; j < d->n_derivs; ++j) {
                    if (d->derivs.sources[i_layer][j])
                         a[j][ii_input] = 0.;
               }

               for (j = 0; j <= i_layer; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         if (! d->derivs.sources[i_layer][k])
                              continue;

                         a[k][ii_input] += d->omega0_l[j][k] * d->omega0_a[j];
                         a[k][ii_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                         if (d->options & XRTM_OPTION_DELTA_M)
                              a[k][ii_input] += d->coef0_l[j][k][0][d->n_coef2] * d->coef0_a[j][0][d->n_coef2];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = i_four; m < d->n_coef2; ++m) {
                                   a[k][ii_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }
          }


          if (i_input == 0 || i_input == 1 || i_input == 2 || i_input == 3) {
               if (i_input == 0)
                    b = R_p_l;
               else
               if (i_input == 1)
                    b = T_p_l;
               else
               if (i_input == 2)
                    b = R_m_l;
               else
                    b = T_m_l;

               for (j = 0; j < d->n_derivs; ++j) {
                    if (d->derivs.layers[i_layer][j])
                         b[j][ii_input][jj_input] = 0.;
               }

               for (j = 0; j <= i_layer; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         if (! d->derivs.layers[i_layer][k])
                              continue;

                         b[k][ii_input][jj_input] += d->omega0_l[j][k] * d->omega0_a[j];
                         b[k][ii_input][jj_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                         if (d->options & XRTM_OPTION_DELTA_M)
                              b[k][ii_input][jj_input] += d->coef0_l[j][k][0][d->n_coef2] * d->coef0_a[j][0][d->n_coef2];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = i_four; m < d->n_coef2; ++m) {
                                   b[k][ii_input][jj_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void get_total_TOA_R_S_U_V_tl_with_ad
     (xrtm_data *d, int i_four, int solver,
      double **Rs_qq, double *S_s, double ***Rs_qq_l, double **S_s_l,
      double **R_m, double *S_p, double ***R_m_l, double **S_p_l,
      int surface, int flag, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double **Rs_qq_a;
     double *S_s_a;

     double **R_m_a;
     double *S_p_a;

     double **a;
     double ***b;


     if (! flags_or(d->derivs.sources[d->n_layers], d->n_derivs))
          return;


     save_tree_decode_s(&save_tree, "get_total_TOA_R_S_U_V");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     Rs_qq_a = get_work1(&work, WORK_DXX);
     S_s_a   = get_work1(&work, WORK_DX);

     R_m_a   = get_work1(&work, WORK_DXX);
     S_p_a   = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dvec_zero(S_p_a, d->n_quad_v_x);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = d->n_quad_v_x * d->n_quad_v_x + d->n_quad_v_x;

     for (i = 0; i < n_inputs; ++i) {
          if (i < d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 0;
               ii = i;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < d->n_quad_v_x * d->n_quad_v_x + d->n_quad_v_x) {
               i_input = 1;
               ii = i - d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_opt_props0_all_a(d, d->n_coef2, work);

          if (d->options & XRTM_OPTION_DELTA_M)
               init_opt_props_all_a(d, d->n_coef2, work);

          init_beam_params_all_a(d, work);

          dmat_zero(Rs_qq_a, d->n_quad_v_x, d->n_quad_v_x);
          dvec_zero(S_s_a,   d->n_quad_v_x);


          if (i_input == 0)
               R_m_a[ii_input][jj_input] = 1.;
          else
               S_p_a[ii_input] = 1.;


          get_total_TOA_R_S_U_V_a(d, i_four, solver, Rs_qq, S_s, Rs_qq_a, S_s_a, R_m, S_p, R_m_a, S_p_a, surface, flag, save_tree, &work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 1) {
               if (i_input == 1)
                    a = S_p_l;

               for (j = 0; j < d->n_derivs; ++j) {
                    if (d->derivs.sources[d->n_layers][j])
                         a[j][ii_input] = 0.;
               }

               for (j = 0; j < d->n_layers; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         if (! d->derivs.sources[d->n_layers][k])
                              continue;

                         a[k][ii_input] += d->omega0_l[j][k] * d->omega0_a[j];
                         a[k][ii_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                         if (d->options & XRTM_OPTION_DELTA_M)
                              a[k][ii_input] += d->coef0_l[j][k][0][d->n_coef2] * d->coef0_a[j][0][d->n_coef2];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = i_four; m < d->n_coef2; ++m) {
                                   a[k][ii_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }

               if (surface) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         for (l = 0; l < d->n_quad_v_x; ++l) {
                              if (d->derivs.sources[d->n_layers][k])
                                   a[k][ii_input] += S_s_l[k][l] * S_s_a[l];

                              for (m = 0; m < d->n_quad_v_x; ++m) {
                                   if (d->derivs.layers[d->n_layers][k])
                                        a[k][ii_input] += Rs_qq_l[k][l][m] * Rs_qq_a[l][m];
                              }
                         }
                    }
               }
          }


          if (i_input == 0) {
               if (i_input == 0)
                    b = R_m_l;

               for (j = 0; j < d->n_derivs; ++j) {
                    if (d->derivs.layers[d->n_layers][j])
                         b[j][ii_input][jj_input] = 0.;
               }

               for (j = 0; j < d->n_layers; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         if (! d->derivs.sources[d->n_layers][k])
                              continue;

                         b[k][ii_input][jj_input] += d->omega0_l[j][k] * d->omega0_a[j];
                         b[k][ii_input][jj_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                         if (d->options & XRTM_OPTION_DELTA_M)
                              b[k][ii_input][jj_input] += d->coef0_l[j][k][0][d->n_coef2] * d->coef0_a[j][0][d->n_coef2];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = i_four; m < d->n_coef2; ++m) {
                                   b[k][ii_input][jj_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }

               if (surface) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         for (l = 0; l < d->n_quad_v_x; ++l) {
                              if (d->derivs.sources[d->n_layers][k])
                                   b[k][ii_input][jj_input] += S_s_l[k][l] * S_s_a[l];

                              for (m = 0; m < d->n_quad_v_x; ++m) {
                                   if (d->derivs.layers[d->n_layers][k])
                                        b[k][ii_input][jj_input] += Rs_qq_l[k][l][m] * Rs_qq_a[l][m];
                              }
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void get_total_R_T_S_U_W_V_tl_with_ad
     (xrtm_data *d, int i_four, int solver,
      double **R_p, double **T_p, double **R_m, double **T_m,
      double *S_p, double *S_m,
      double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
      double **S_p_l, double **S_m_l,
      int flag, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double **R_p_a;
     double **T_p_a;
     double **R_m_a;
     double **T_m_a;
     double *S_p_a;
     double *S_m_a;

     double **a;
     double ***b;


     if (! flags_or(d->derivs.sources[d->n_layers], d->n_derivs))
          return;


     save_tree_decode_s(&save_tree, "get_total_R_T_S_U_W_V");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R_p_a = get_work1(&work, WORK_DXX);
     T_p_a = get_work1(&work, WORK_DXX);
     R_m_a = get_work1(&work, WORK_DXX);
     T_m_a = get_work1(&work, WORK_DXX);
     S_p_a = get_work1(&work, WORK_DX);
     S_m_a = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T_p_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(R_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dmat_zero(T_m_a, d->n_quad_v_x, d->n_quad_v_x);
     dvec_zero(S_p_a, d->n_quad_v_x);
     dvec_zero(S_m_a, d->n_quad_v_x);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 4 * d->n_quad_v_x * d->n_quad_v_x + 2 * d->n_quad_v_x;

     for (i = 0; i < n_inputs; ++i) {
          if (i < d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 0;
               ii = i;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 2 * d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 1;
               ii = i - d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 3 * d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 2;
               ii = i - 2 * d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 4 * d->n_quad_v_x * d->n_quad_v_x) {
               i_input = 3;
               ii = i - 3 * d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii / d->n_quad_v_x;
               jj_input = ii % d->n_quad_v_x;
          }
          else
          if (i < 4 * d->n_quad_v_x * d->n_quad_v_x + d->n_quad_v_x) {
               i_input = 4;
               ii = i - 4 * d->n_quad_v_x * d->n_quad_v_x;
               ii_input = ii;
          }
          else {
               i_input = 5;
               ii = i - (4 * d->n_quad_v_x * d->n_quad_v_x + d->n_quad_v_x);
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_opt_props0_all_a(d, d->n_coef2, work);

          if (d->options & XRTM_OPTION_DELTA_M)
               init_opt_props_all_a(d, d->n_coef2, work);

          init_beam_params_all_a(d, work);


          if (i_input == 0)
               R_p_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 1)
               T_p_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 2)
               R_m_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 3)
               T_m_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 4)
               S_p_a[ii_input] = 1.;
          else
               S_m_a[ii_input] = 1.;


          get_total_R_T_S_U_W_V_a(d, i_four, solver, R_p, T_p, R_m, T_m, S_p, S_m, R_p_a, T_p_a, R_m_a, T_m_a, S_p_a, S_m_a, flag, save_tree, &work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 4 || i_input == 5) {
               if (i_input == 4)
                    a = S_p_l;
               else
                    a = S_m_l;

               for (j = 0; j < d->n_derivs; ++j) {
                    if (d->derivs.sources[d->n_layers - 1][j])
                         a[j][ii_input] = 0.;
               }

               for (j = 0; j < d->n_layers; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         if (! d->derivs.sources[d->n_layers - 1][k])
                              continue;

                         a[k][ii_input] += d->omega0_l[j][k] * d->omega0_a[j];
                         a[k][ii_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                         if (d->options & XRTM_OPTION_DELTA_M)
                              a[k][ii_input] += d->coef0_l[j][k][0][d->n_coef2] * d->coef0_a[j][0][d->n_coef2];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = i_four; m < d->n_coef2; ++m) {
                                   a[k][ii_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }
          }


          if (i_input == 0 || i_input == 1 || i_input == 2 || i_input == 3) {
               if (i_input == 0)
                    b = R_p_l;
               else
               if (i_input == 1)
                    b = T_p_l;
               else
               if (i_input == 2)
                    b = R_m_l;
               else
                    b = T_m_l;

               for (j = 0; j < d->n_derivs; ++j) {
                    if (d->derivs.sources[d->n_layers - 1][j])
                         b[j][ii_input][jj_input] = 0.;
               }

               for (j = 0; j < d->n_layers; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         if (! d->derivs.sources[d->n_layers - 1][k])
                              continue;

                         b[k][ii_input][jj_input] += d->omega0_l[j][k] * d->omega0_a[j];
                         b[k][ii_input][jj_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                         if (d->options & XRTM_OPTION_DELTA_M)
                              b[k][ii_input][jj_input] += d->coef0_l[j][k][0][d->n_coef2] * d->coef0_a[j][0][d->n_coef2];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = i_four; m < d->n_coef2; ++m) {
                                   b[k][ii_input][jj_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void update_diff_bound_input_tl_with_ad(xrtm_data *d, int i_four, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;

     int i_input;
     int n_inputs;

     int ii_input;

     double *In_p_a;
     double *I1_m_a;

     double **a;


     save_tree_decode_s(&save_tree, "update_diff_bound_input");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     In_p_a = get_work1(&work, WORK_DD);
     I1_m_a = get_work1(&work, WORK_DD);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(In_p_a, d->n_quad_v + d->n_umus_v, 0.);
     init_array1_d(I1_m_a, d->n_quad_v + d->n_umus_v, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * (d->n_quad_v + d->n_umus_v);

     for (i = 0; i < n_inputs; ++i) {
          if (i < d->n_quad_v + d->n_umus_v) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }
          else {
               i_input = 1;
               ii = i - (d->n_quad_v + d->n_umus_v);
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_opt_props_all_a(d, d->n_coef2, work);

          init_beam_params_all_a(d, work);

          for (j = 0; j < d->n_kernels; ++j)
               d->kernel_ampfac_a[j] = 0.;


          if (i_input == 0)
               In_p_a[ii_input] = 1.;
          else
               I1_m_a[ii_input] = 1.;


          update_diff_bound_input_a(d, i_four, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = d->In_p_l[i_four];
               else
                    a = d->I1_m_l[i_four];

               for (j = 0; j < d->n_derivs; ++j)
                    a[j][ii_input] = 0.;

               for (j = 0; j < d->n_layers; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         a[k][ii_input] += d->ltau_l [j][k] * d->ltau_a [j];
                         a[k][ii_input] += d->omega_l[j][k] * d->omega_a[j];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = i_four; m < d->n_coef2; ++m) {
                                   a[k][ii_input] += d->coef_l[j][k][l][m] * d->coef_a[j][l][m];
                              }
                         }
                    }
               }

               for (j = 0; j < d->n_kernels; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         a[k][ii_input] += d->kernel_ampfac_l[j][k] * d->kernel_ampfac_a[j];
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void get_single_tl_with_ad(xrtm_data *d, int n_phis, double **phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;
/*
     int n_quad_x;
*/
     int n_mus2;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;
     int ll_input;
     int kk_input;

     double ****I_p_a;
     double ****I_m_a;

     double *****a;


     save_tree_decode_s(&save_tree, "get_single");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     n_quad_x = d->n_quad + d->n_umus;
*/
     if (d->n_umus == 0)
          n_mus2   = d->n_quad;
     else
          n_mus2   = d->n_umus;
/*
     else
          n_mus2   = d->n_quad + d->n_umus;
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_p_a = get_work_d4(&work, d->n_ulevels, n_mus2, n_phis, d->n_stokes);
     I_m_a = get_work_d4(&work, d->n_ulevels, n_mus2, n_phis, d->n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array4_d(I_p_a, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);
     init_array4_d(I_m_a, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * (d->n_ulevels * n_mus2 * n_phis * d->n_stokes);

     for (i = 0; i < n_inputs; ++i) {
          if (i < d->n_ulevels * n_mus2 * n_phis * d->n_stokes) {
               i_input = 0;
               ii = i;
               ii_input = ii / (n_mus2 * n_phis * d->n_stokes);
               jj_input = ii % (n_mus2 * n_phis * d->n_stokes) / (n_phis * d->n_stokes);
               kk_input = ii % (n_mus2 * n_phis * d->n_stokes) % (n_phis * d->n_stokes) / d->n_stokes;
               ll_input = ii % (n_mus2 * n_phis * d->n_stokes) % (n_phis * d->n_stokes) % d->n_stokes;
          }
          else {
               i_input = 1;
               ii = i - (d->n_ulevels * n_mus2 * n_phis * d->n_stokes);
               ii_input = ii / (n_mus2 * n_phis * d->n_stokes);
               jj_input = ii % (n_mus2 * n_phis * d->n_stokes) / (n_phis * d->n_stokes);
               kk_input = ii % (n_mus2 * n_phis * d->n_stokes) % (n_phis * d->n_stokes) / d->n_stokes;
               ll_input = ii % (n_mus2 * n_phis * d->n_stokes) % (n_phis * d->n_stokes) % d->n_stokes;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_opt_props0_all_a(d, d->n_coef, work);

          if (d->options & XRTM_OPTION_DELTA_M)
               init_opt_props_all_a(d, d->n_coef2, work);

          init_beam_params_all_a(d, work);

          init_diff_bound_input_a(d, 0, work);

          for (j = 0; j < d->n_kernels; ++j)
               d->kernel_ampfac_a[j] = 0.;


          if (i_input == 0)
               I_p_a[ii_input][jj_input][kk_input][ll_input] = 1.;
          else
               I_m_a[ii_input][jj_input][kk_input][ll_input] = 1.;


          get_single_a(d, n_phis, phis, I_p, I_m, I_p_a, I_m_a, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = I_p_l;
               else
                    a = I_m_l;

               for (j = 0; j < d->n_derivs; ++j)
                    a[ii_input][j][jj_input][kk_input][ll_input] = 0.;

               for (j = 0; j < d->n_layers; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         a[ii_input][k][jj_input][kk_input][ll_input] += d->omega0_l[j][k] * d->omega0_a[j];
                         a[ii_input][k][jj_input][kk_input][ll_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = 0; m < d->n_coef; ++m) {
                                   a[ii_input][k][jj_input][kk_input][ll_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }

               for (j = 0; j < d->n_kernels; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         a[ii_input][k][jj_input][kk_input][ll_input] += d->kernel_ampfac_l[j][k] * d->kernel_ampfac_a[j];
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void fourier_get_add_both_ways_tl_with_ad(xrtm_data *d, int i_four, int solver, double **I_p, double **I_m, double ***I_p_l, double ***I_m_l, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;

     int i_input;
     int n_inputs;

     int ii_input;

     double **I_p_a;
     double **I_m_a;

     double **a;


     save_tree_decode_s(&save_tree, "fourier_get_add_both_ways");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_p_a = get_work_d2(&work, d->n_ulevels, d->n_quad_v_x);
     I_m_a = get_work_d2(&work, d->n_ulevels, d->n_quad_v_x);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_p_a, d->n_ulevels, d->n_quad_v_x, 0.);
     init_array2_d(I_m_a, d->n_ulevels, d->n_quad_v_x, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 4 * d->n_quad_v_x;

     for (i = 0; i < n_inputs; ++i) {
          if (i < d->n_quad_v_x) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }
          else
          if (i < 2 * d->n_quad_v_x) {
               i_input = 1;
               ii = i - d->n_quad_v_x;
               ii_input = ii;
          }
          else
          if (i < 3 * d->n_quad_v_x) {
               i_input = 2;
               ii = i - 2 * d->n_quad_v_x;
               ii_input = ii;
          }
          else {
               i_input = 3;
               ii = i - 3 * d->n_quad_v_x;
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_opt_props0_all_a(d, d->n_coef2, work);

          if (d->options & XRTM_OPTION_DELTA_M)
               init_opt_props_all_a(d, d->n_coef2, work);

          init_beam_params_all_a(d, work);

          init_diff_bound_input_a(d, i_four, work);

          for (j = 0; j < d->n_kernels; ++j)
               d->kernel_ampfac_a[j] = 0.;


          if (i_input == 0)
               I_p_a[0][ii_input] = 1.;
          else
          if (i_input == 1)
               I_p_a[1][ii_input] = 1.;
          else
          if (i_input == 2)
               I_m_a[0][ii_input] = 1.;
          else
               I_m_a[1][ii_input] = 1.;


          fourier_get_add_both_ways_a(d, i_four, solver, I_p, I_m, I_p_a, I_m_a, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1 || i_input == 2 || i_input == 3) {
               if (i_input == 0)
                    a = I_p_l[0];
               else
               if (i_input == 1)
                    a = I_p_l[1];
               else
               if (i_input == 2)
                    a = I_m_l[0];
               else
                    a = I_m_l[1];

               for (j = 0; j < d->n_derivs; ++j)
                    a[j][ii_input] = 0.;

               for (j = 0; j < d->n_layers; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         a[k][ii_input] += d->omega0_l[j][k] * d->omega0_a[j];
                         a[k][ii_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                         if (d->options & XRTM_OPTION_DELTA_M)
                              a[k][ii_input] += d->coef0_l[j][k][0][d->n_coef2] * d->coef0_a[j][0][d->n_coef2];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = i_four; m < d->n_coef2; ++m) {
                                   a[k][ii_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }

               for (j = 0; j < d->n_kernels; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         a[k][ii_input] += d->kernel_ampfac_l[j][k] * d->kernel_ampfac_a[j];
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void get_solution_internal_tl_with_ad(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double *flux_div, double **flux_div_l, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;
     int m;

     int n_mus2;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;
     int kk_input;
     int ll_input;

     double ****I_p_a;
     double ****I_m_a;

     double *****a;


     save_tree_decode_s(&save_tree, "get_solution_internal");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (d->n_umus == 0)
          n_mus2   = d->n_quad;
     else
          n_mus2   = d->n_umus;
/*
     else
          n_mus2   = d->n_quad + d->n_umus;
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     I_p_a = get_work_d4(&work, d->n_ulevels, n_mus2, n_phis, d->n_stokes);
     I_m_a = get_work_d4(&work, d->n_ulevels, n_mus2, n_phis, d->n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array4_d(I_p_a, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);
     init_array4_d(I_m_a, d->n_ulevels, n_mus2, n_phis, d->n_stokes, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * d->n_ulevels * n_mus2 * n_phis * d->n_stokes;

     for (i = 0; i < n_inputs; ++i) {
          if (i < d->n_ulevels * n_mus2 * n_phis * d->n_stokes) {
               i_input = 0;
               ii = i;
               ii_input = ii / (n_mus2 * n_phis * d->n_stokes);
               jj_input = ii % (n_mus2 * n_phis * d->n_stokes) / (n_phis * d->n_stokes);
               kk_input = ii % (n_mus2 * n_phis * d->n_stokes) % (n_phis * d->n_stokes) / d->n_stokes;
               ll_input = ii % (n_mus2 * n_phis * d->n_stokes) % (n_phis * d->n_stokes) % d->n_stokes;
          }
          else {
               i_input = 1;
               ii = i - d->n_ulevels * n_mus2 * n_phis * d->n_stokes;
               ii_input = ii / (n_mus2 * n_phis * d->n_stokes);
               jj_input = ii % (n_mus2 * n_phis * d->n_stokes) / (n_phis * d->n_stokes);
               kk_input = ii % (n_mus2 * n_phis * d->n_stokes) % (n_phis * d->n_stokes) / d->n_stokes;
               ll_input = ii % (n_mus2 * n_phis * d->n_stokes) % (n_phis * d->n_stokes) % d->n_stokes;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_opt_props0_all_a(d, d->n_coef, work);

          if (d->options & XRTM_OPTION_DELTA_M)
               init_opt_props_all_a(d, d->n_coef2, work);

          init_beam_params_all_a(d, work);

          init_diff_bound_input_all_a(d, d->n_four, work);

          for (j = 0; j < d->n_kernels; ++j)
               d->kernel_ampfac_a[j] = 0.;


          if (i_input == 0)
               I_p_a[ii_input][jj_input][kk_input][ll_input] = 1.;
          else
               I_m_a[ii_input][jj_input][kk_input][ll_input] = 1.;


          get_solution_internal_a(d, solver, solutions, n_phis, phis, I_p, I_m, I_p_a, I_m_a, mean_p, mean_m, NULL, NULL, flux_p, flux_m, NULL, NULL, NULL, NULL, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = I_p_l;
               else
                    a = I_m_l;

               for (j = 0; j < d->n_derivs; ++j)
                    a[ii_input][j][jj_input][kk_input][ll_input] = 0.;

               for (j = 0; j < d->n_kernels; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         a[ii_input][k][jj_input][kk_input][ll_input] += d->kernel_ampfac_l[j][k] * d->kernel_ampfac_a[j];
                    }
               }

               for (j = 0; j < d->n_layers; ++j) {
                    for (k = 0; k < d->n_derivs; ++k) {
                         a[ii_input][k][jj_input][kk_input][ll_input] += d->omega0_l[j][k] * d->omega0_a[j];
                         a[ii_input][k][jj_input][kk_input][ll_input] += d->ltau0_l [j][k] * d->ltau0_a [j];

                         for (l = 0; l < d->n_elem; ++l) {
                              for (m = 0; m < d->n_coef; ++m) {
                                   a[ii_input][k][jj_input][kk_input][ll_input] += d->coef0_l[j][k][l][m] * d->coef0_a[j][l][m];
                              }
                         }
                    }
               }
          }
     }
}
