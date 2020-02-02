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
void build_source_vectors_1n_tl_with_ad(int n_quad, int n_stokes, int n_derivs, double *qx_v,
                                        double F_0, double omega, double *omega_l,
                                        double as_0, double *as_0_l,
                                        double  *P_q0_mm, double  *P_q0_pm,
                                        double  **tpr, double  **tmr, double  **gamma,
                                        double  *F_p, double  *F_m,
                                        double **P_q0_mm_l, double **P_q0_pm_l,
                                        double ***tpr_l, double ***tmr_l, double ***gamma_l,
                                        double **F_p_l, double **F_m_l,
                                        uchar *derivs_layers, uchar *derivs_beam_down,
                                        save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int n_quad_v;

     int i_input;
     int n_inputs;

     int ii_input;

     double omega_a;
     double as_0_a;

     double *P_q0_mm_a;
     double *P_q0_pm_a;
     double **tpr_a;
     double **tmr_a;
     double **gamma_a;

     double *F_p_a;
     double *F_m_a;

     double **a;


     if (! flags_or(derivs_layers, n_derivs) && ! flags_or(derivs_beam_down, n_derivs))
          return;


     save_tree_decode_s(&save_tree, "build_source_vectors_1n");


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     P_q0_mm_a = get_work1(&work, WORK_DX);
     P_q0_pm_a = get_work1(&work, WORK_DX);
     tpr_a     = get_work1(&work, WORK_DXX);
     tmr_a     = get_work1(&work, WORK_DXX);
     gamma_a   = get_work1(&work, WORK_DXX);

     F_p_a     = get_work1(&work, WORK_DX);
     F_m_a     = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(F_p_a, n_quad_v);
     dvec_zero(F_m_a, n_quad_v);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * n_quad_v;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad_v) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }
          else {
               i_input = 1;
               ii = i - n_quad_v;
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          omega_a = 0.;
          as_0_a  = 0.;

          dvec_zero(P_q0_mm_a, n_quad_v);
          dvec_zero(P_q0_pm_a, n_quad_v);
          dmat_zero(tpr_a, n_quad_v, n_quad_v);
          dmat_zero(tmr_a, n_quad_v, n_quad_v);
          dmat_zero(gamma_a, n_quad_v, n_quad_v);


          if (i_input == 0)
               F_p_a[ii_input] = 1.;
          else
               F_m_a[ii_input] = 1.;


          build_source_vectors_1n_a(n_quad, n_stokes, qx_v, F_0, omega, &omega_a, as_0, &as_0_a, P_q0_mm, P_q0_pm, tpr, tmr, gamma, F_p, F_m, P_q0_mm_a, P_q0_pm_a, tpr_a, tmr_a, gamma_a, F_p_a, F_m_a, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = F_p_l;
               else
                    a = F_m_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs_layers[j] && ! derivs_beam_down[j])
                         continue;

                    a[j][ii_input] = 0.;

                    if (derivs_layers[j]) {
                         a[j][ii_input] += omega_l[j] * omega_a;

                         for (k = 0; k < n_quad_v; ++k) {
                              a[j][ii_input] += P_q0_mm_l[j][k] * P_q0_mm_a[k];
                              a[j][ii_input] += P_q0_pm_l[j][k] * P_q0_pm_a[k];

                              for (l = 0; l < n_quad_v; ++l) {
                                   a[j][ii_input] += tpr_l[j][k][l] * tpr_a[k][l];
                                   a[j][ii_input] += tmr_l[j][k][l] * tmr_a[k][l];
                                   a[j][ii_input] += gamma_l[j][k][l] * gamma_a[k][l];
                              }
                         }
                    }

                    if (derivs_beam_down[j])
                         a[j][ii_input] += as_0_l [j] * as_0_a;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_global_source_tl_with_ad2(int n_quad_v, int n_derivs, double atran, double *atran_l,
                                double **R_p, double **T_p, double **R_m, double **T_m,
                                double *F_p, double *F_m, double *S_p, double *S_m,
                                double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                                double **F_p_l, double **F_m_l, double **S_p_l, double **S_m_l,
                                uchar *derivs_layers, uchar *derivs_beam_down, work_data work) {
     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;

     double atran_a;

     double **R_p_a;
     double **T_p_a;
     double **R_m_a;
     double **T_m_a;
     double *F_p_a;
     double *F_m_a;

     double *S_p_a;
     double *S_m_a;

     double **a;


     if (! flags_or(derivs_layers, n_derivs) && ! flags_or(derivs_beam_down, n_derivs))
          return;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R_p_a = get_work1(&work, WORK_DXX);
     T_p_a = get_work1(&work, WORK_DXX);
     R_m_a = get_work1(&work, WORK_DXX);
     T_m_a = get_work1(&work, WORK_DXX);

     F_p_a = get_work1(&work, WORK_DX);
     F_m_a = get_work1(&work, WORK_DX);

     S_p_a = get_work1(&work, WORK_DX);
     S_m_a = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(S_p_a, n_quad_v);
     dvec_zero(S_m_a, n_quad_v);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * n_quad_v;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad_v) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }
          else {
               i_input = 1;
               ii = i - n_quad_v;
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          atran_a = 0.;

          dmat_zero(R_p_a, n_quad_v, n_quad_v);
          dmat_zero(T_p_a, n_quad_v, n_quad_v);
          dmat_zero(R_m_a, n_quad_v, n_quad_v);
          dmat_zero(T_m_a, n_quad_v, n_quad_v);
          dvec_zero(F_p_a, n_quad_v);
          dvec_zero(F_m_a, n_quad_v);


          if (i_input == 0)
               S_p_a[ii_input] = 1.;
          else
               S_m_a[ii_input] = 1.;


          build_global_source_a2(n_quad_v, atran, &atran_a, R_p, T_p, R_m, T_m, F_p, F_m, S_p, S_m, R_p_a, T_p_a, R_m_a, T_m_a, F_p_a, F_m_a, S_p_a, S_m_a, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = S_p_l;
               else
                    a = S_m_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs_layers[j] && ! derivs_beam_down[j])
                         continue;

                    a[j][ii_input] = 0.;

                    if (derivs_layers[j]) {
                         for (k = 0; k < n_quad_v; ++k) {
                              for (l = 0; l < n_quad_v; ++l) {
                                   a[j][ii_input] += R_p_l[j][k][l] * R_p_a[k][l];
                                   a[j][ii_input] += T_p_l[j][k][l] * T_p_a[k][l];
                                   a[j][ii_input] += R_m_l[j][k][l] * R_m_a[k][l];
                                   a[j][ii_input] += T_m_l[j][k][l] * T_m_a[k][l];
                              }
                         }
                    }

                    if (derivs_beam_down[j]) {
                         for (k = 0; k < n_quad_v; ++k) {
                              a[j][ii_input] += F_p_l[j][k] * F_p_a[k];
                              a[j][ii_input] += F_m_l[j][k] * F_m_a[k];
                         }
                    }

                    if (derivs_beam_down[j])
                         a[j][ii_input] += atran_l[j] * atran_a;
               }
          }
     }
}
