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
void rtm_eig_bvp_tl_with_ad(int i_four, int n_quad, int n_stokes, int n_derivs, int n_layers,
       double qf, double *qx_v, double *qw_v, double F_0, double mu_0,
       int *ulevels, double *utaus, int n_ulevels,
       double *umus, int n_umus, double *planck,
       double *omega, double **omega_l, double *ltau, double **ltau_l,
       double **Rs_qq, double ***Rs_qq_l,
       double *Rs_u0, double **Rs_u0_l, double **Rs_uq, double ***Rs_uq_l,
       double *btran, double **btran_l,
       double *as_0, double **as_0_l, double *atran, double **atran_l,
       double **P_q0_mm, double **P_q0_pm, double **P_u0_mm, double **P_u0_pm,
       double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm,
       double ***r_p, double ***t_p, double ***r_m, double ***t_m,
       double ***P_q0_mm_l, double ***P_q0_pm_l, double ***P_u0_mm_l, double ***P_u0_pm_l,
       double ****P_uq_pp_l, double ****P_uq_mp_l, double ****P_uq_mm_l, double ****P_uq_pm_l,
       double ****r_p_l, double ****t_p_l, double ****r_m_l, double ****t_m_l,
       double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
       double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
       int sfi, int surface, int thermal, int upwelling,
       int downwelling, int utau_output, int vector,
       int eigen_solver_real, int eigen_solver_complex,
       uchar **derivs_layers, uchar **derivs_beam_down, save_tree_data save_tree, work_data work) {

     uchar *derivs_layers_union;
     uchar *derivs_beam_down_union;

     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_layer;

     int n_quad_x;
     int n_quad_v;
     int n_quad_v_x;
     int n_umus_v;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double *omega_a;
     double *ltau_a;
     double *as_0_a;
     double *btran_a;
     double *atran_a;

     double **Rs_qq_a;
     double *Rs_u0_a;
     double **Rs_uq_a;

     double **P_q0_mm_a;
     double **P_q0_pm_a;
     double **P_u0_mm_a;
     double **P_u0_pm_a;

     double ***P_uq_pp_a;
     double ***P_uq_mp_a;
     double ***P_uq_mm_a;
     double ***P_uq_pm_a;

     double ***r_p_a;
     double ***t_p_a;
     double ***r_m_a;
     double ***t_m_a;

     double *I1_m_a;
     double *In_p_a;

     double **I_p_a;
     double **I_m_a;

     double ***a;

/*
     if (! flags_or(derivs_layers, n_derivs) && ! flags_or(derivs_beam_down, n_derivs))
          return;
*/

     calc_quad_and_umus(n_quad, n_umus, n_stokes, &n_quad_x, NULL, &n_quad_v, &n_quad_v_x, NULL, &n_umus_v, sfi);


     save_tree_decode_s(&save_tree, "rtm_eig_bvp");


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
     omega_a   = get_work_d1(&work, n_layers);
     ltau_a    = get_work_d1(&work, n_layers);
     Rs_qq_a   = get_work_d2(&work, n_quad_v_x, n_quad_v_x);
     Rs_u0_a   = get_work_d1(&work, n_umus_v);
     Rs_uq_a   = get_work_d2(&work, n_umus_v, n_quad_v_x);
     btran_a   = get_work_d1(&work, n_layers + 1);
     as_0_a    = get_work_d1(&work, n_layers);
     atran_a   = get_work_d1(&work, n_layers);
     P_q0_mm_a = get_work_d2(&work, n_layers, n_quad_v_x);
     P_q0_pm_a = get_work_d2(&work, n_layers, n_quad_v_x);
     P_u0_mm_a = get_work_d2(&work, n_layers, n_umus_v);
     P_u0_pm_a = get_work_d2(&work, n_layers, n_umus_v);
     P_uq_pp_a = get_work_d3(&work, n_layers, n_umus_v, n_quad_v_x);
     P_uq_mp_a = get_work_d3(&work, n_layers, n_umus_v, n_quad_v_x);
     P_uq_mm_a = get_work_d3(&work, n_layers, n_umus_v, n_quad_v_x);
     P_uq_pm_a = get_work_d3(&work, n_layers, n_umus_v, n_quad_v_x);
     r_p_a     = get_work_d3(&work, n_layers, n_quad_v_x, n_quad_v_x);
     t_p_a     = get_work_d3(&work, n_layers, n_quad_v_x, n_quad_v_x);
     r_m_a     = get_work_d3(&work, n_layers, n_quad_v_x, n_quad_v_x);
     t_m_a     = get_work_d3(&work, n_layers, n_quad_v_x, n_quad_v_x);
     I1_m_a    = get_work_d1(&work, n_quad_v_x);
     In_p_a    = get_work_d1(&work, n_quad_v_x);

     I_p_a     = get_work_d2(&work, n_ulevels, n_quad_v_x);
     I_m_a     = get_work_d2(&work, n_ulevels, n_quad_v_x);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_p_a, n_ulevels, n_quad_v_x, 0.);
     init_array2_d(I_m_a, n_ulevels, n_quad_v_x, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * n_ulevels * n_quad_v_x;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_ulevels * n_quad_v_x) {
               i_input = 0;
               ii = i;
               ii_input = ii / n_quad_v_x;
               jj_input = ii % n_quad_v_x;
          }
          else {
               i_input = 1;
               ii = i - n_ulevels * n_quad_v_x;
               ii_input = ii / n_quad_v_x;
               jj_input = ii % n_quad_v_x;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array1_d(omega_a,   n_layers, 0.);
          init_array1_d(ltau_a,    n_layers, 0.);
          init_array2_d(Rs_qq_a,   n_quad_v_x, n_quad_v_x, 0.);
          init_array1_d(Rs_u0_a,   n_umus_v, 0.);
          init_array2_d(Rs_uq_a,   n_umus_v, n_quad_v_x, 0.);
          init_array1_d(btran_a,   n_layers + 1, 0.);
          init_array1_d(as_0_a ,   n_layers, 0.);
          init_array1_d(atran_a,   n_layers, 0.);

          init_array2_d(P_q0_mm_a, n_layers, n_quad_v_x, 0.);
          init_array2_d(P_q0_pm_a, n_layers, n_quad_v_x, 0.);
          init_array2_d(P_u0_mm_a, n_layers, n_umus_v, 0.);
          init_array2_d(P_u0_pm_a, n_layers, n_umus_v, 0.);
          init_array3_d(P_uq_pp_a, n_layers, n_umus_v, n_quad_v_x, 0.);
          init_array3_d(P_uq_mp_a, n_layers, n_umus_v, n_quad_v_x, 0.);
          init_array3_d(P_uq_mm_a, n_layers, n_umus_v, n_quad_v_x, 0.);
          init_array3_d(P_uq_pm_a, n_layers, n_umus_v, n_quad_v_x, 0.);
          init_array3_d(r_p_a,     n_layers, n_quad_v_x, n_quad_v_x, 0.);
          init_array3_d(t_p_a,     n_layers, n_quad_v_x, n_quad_v_x, 0.);
          init_array3_d(r_m_a,     n_layers, n_quad_v_x, n_quad_v_x, 0.);
          init_array3_d(t_m_a,     n_layers, n_quad_v_x, n_quad_v_x, 0.);
          init_array1_d(I1_m_a,    n_quad_v_x, 0.);
          init_array1_d(In_p_a,    n_quad_v_x, 0.);


          if (i_input == 0)
               I_p_a[ii_input][jj_input] = 1.;
          else
               I_m_a[ii_input][jj_input] = 1.;


          rtm_eig_bvp_a(i_four, n_quad, n_stokes, n_layers, qf, qx_v, qw_v, F_0, mu_0, ulevels, utaus, n_ulevels, umus, n_umus, planck, omega, omega_a, ltau, ltau_a, Rs_qq, Rs_qq_a, Rs_u0, Rs_u0_a, Rs_uq, Rs_uq_a, btran, btran_a, as_0, as_0_a, atran, atran_a, P_q0_mm, P_q0_pm, P_u0_mm, P_u0_pm, P_uq_pp, P_uq_mp, P_uq_mm, P_uq_pm, r_p, t_p, r_m, t_m, P_q0_mm_a, P_q0_pm_a, P_u0_mm_a, P_u0_pm_a, P_uq_pp_a, P_uq_mp_a, P_uq_mm_a, P_uq_pm_a, r_p_a, t_p_a, r_m_a, t_m_a, I1_m, I1_m_a, In_p, In_p_a, I_p, I_m, I_p_a, I_m_a, sfi, surface, thermal, upwelling, downwelling, utau_output, vector, derivs_layers_union, derivs_beam_down_union, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = I_p_l;
               else
                    a = I_m_l;

               for (j = 0; j < n_derivs; ++j) {
                    a[ii_input][j][jj_input] = 0.;

                    for (i_layer = 0; i_layer < n_layers; ++i_layer) {
                         a[ii_input][j][jj_input] += omega_l[i_layer][j] * omega_a[i_layer];
                         a[ii_input][j][jj_input] += ltau_l [i_layer][j] * ltau_a [i_layer];
                    }

                    if (surface) {
                         if (derivs_layers[n_layers][j]) {
                              for (k = 0; k < n_quad_v_x; ++k) {
                                   for (l = 0; l < n_quad_v_x; ++l) {
                                        a[ii_input][j][jj_input] += Rs_qq_l[j][k][l] * Rs_qq_a[k][l];
                                   }
                              }
                              if (sfi && n_umus > 0) {
                                   for (k = 0; k < n_umus_v; ++k) {
                                        a[ii_input][j][jj_input] += Rs_u0_l[j][k] * Rs_u0_a[k];

                                        for (l = 0; l < n_quad_v_x; ++l) {
                                             a[ii_input][j][jj_input] += Rs_uq_l[j][k][l] * Rs_uq_a[k][l];
                                        }
                                   }
                              }
                         }
                    }

                    for (i_layer = 0; i_layer < n_layers; ++i_layer) {
                         a[ii_input][j][jj_input] += btran_l[i_layer][j] * btran_a[i_layer];
                         a[ii_input][j][jj_input] += as_0_l [i_layer][j] * as_0_a [i_layer];
                         a[ii_input][j][jj_input] += atran_l[i_layer][j] * atran_a[i_layer];
                    }

                    a[ii_input][j][jj_input] += btran_l[i_layer][j] * btran_a[i_layer];

                    for (i_layer = 0; i_layer < n_layers; ++i_layer) {
                         if (derivs_layers[i_layer][j]) {
                              for (k = 0; k < n_quad_v_x; ++k) {
                                   a[ii_input][j][jj_input] += P_q0_mm_l[i_layer][j][k] * P_q0_mm_a[i_layer][k];
                                   a[ii_input][j][jj_input] += P_q0_pm_l[i_layer][j][k] * P_q0_pm_a[i_layer][k];
                              }

                              if (sfi && n_umus > 0) {
                                   for (k = 0; k < n_umus_v; ++k) {
                                        a[ii_input][j][jj_input] += P_u0_mm_l[i_layer][j][k] * P_u0_mm_a[i_layer][k];
                                        a[ii_input][j][jj_input] += P_u0_pm_l[i_layer][j][k] * P_u0_pm_a[i_layer][k];

                                        for (l = 0; l < n_quad_v_x; ++l) {
                                             a[ii_input][j][jj_input] += P_uq_pp_l[i_layer][j][k][l] * P_uq_pp_a[i_layer][k][l];
                                             a[ii_input][j][jj_input] += P_uq_mp_l[i_layer][j][k][l] * P_uq_mp_a[i_layer][k][l];
                                             a[ii_input][j][jj_input] += P_uq_mm_l[i_layer][j][k][l] * P_uq_mm_a[i_layer][k][l];
                                             a[ii_input][j][jj_input] += P_uq_pm_l[i_layer][j][k][l] * P_uq_pm_a[i_layer][k][l];
                                        }
                                   }
                              }

                              for (k = 0; k < n_quad_v_x; ++k) {
                                   for (l = 0; l < n_quad_v_x; ++l) {
                                        a[ii_input][j][jj_input] += r_p_l[i_layer][j][k][l] * r_p_a[i_layer][k][l];
                                        a[ii_input][j][jj_input] += t_p_l[i_layer][j][k][l] * t_p_a[i_layer][k][l];
/*
if (vector) {
                                        a[ii_input][j][jj_input] += r_m_l[i_layer][j][k][l] * r_m_a[i_layer][k][l];
                                        a[ii_input][j][jj_input] += t_m_l[i_layer][j][k][l] * t_m_a[i_layer][k][l];
}
*/
                                   }
                              }
                         }
                    }

                    for (k = 0; k < n_quad_v_x; ++k) {
                         a[ii_input][j][jj_input] += I1_m_l[j][k] * I1_m_a[k];
                         a[ii_input][j][jj_input] += In_p_l[j][k] * In_p_a[k];
                    }
               }
          }
     }
}
