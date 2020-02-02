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
void rtm_eig_rts_tl_with_ad(int n_quad, int n_stokes, int n_derivs, double F_0,
                            double *qx_v, double *qw_v, double planck0, double planck1,
                            double omega, double *omega_l, double ltau, double *ltau_l,
                            double as_0, double *as_0_l, double atran, double *atran_l,
                            double  *P_0p, double  *P_0m,
                            double  **r_p, double  **t_p, double  **r_m, double  **t_m,
                            double  **R_p, double  **T_p, double  **R_m, double  **T_m,
                            double  *S_p, double  *S_m,
                            double **P_0p_l, double **P_0m_l,
                            double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
                            double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                            double **S_p_l, double **S_m_l,
                            int symmetric, int thermal, int vector,
                            int eigen_solver_real, int eigen_solver_complex,
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
     int jj_input;

     double omega_a;
     double ltau_a;
     double as_0_a;
     double atran_a;

     double *P_0p_a;
     double *P_0m_a;

     double **r_p_a;
     double **t_p_a;
     double **r_m_a;
     double **t_m_a;

     double **R_p_a;
     double **T_p_a;
     double **R_m_a;
     double **T_m_a;

     double *S_p_a;
     double *S_m_a;

     double **a;
     double ***b;


     if (! flags_or(derivs_layers, n_derivs) && ! flags_or(derivs_beam_down, n_derivs))
          return;


     save_tree_decode_s(&save_tree, "rtm_eig_rts");


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     P_0p_a = get_work1(&work, WORK_DX);
     P_0m_a = get_work1(&work, WORK_DX);

     r_p_a  = get_work1(&work, WORK_DXX);
     t_p_a  = get_work1(&work, WORK_DXX);
     r_m_a  = get_work1(&work, WORK_DXX);
     t_m_a  = get_work1(&work, WORK_DXX);

     R_p_a  = get_work1(&work, WORK_DXX);
     T_p_a  = get_work1(&work, WORK_DXX);
     R_m_a  = get_work1(&work, WORK_DXX);
     T_m_a  = get_work1(&work, WORK_DXX);

     S_p_a  = get_work1(&work, WORK_DX);
     S_m_a  = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R_p_a, n_quad_v, n_quad_v);
     dmat_zero(T_p_a, n_quad_v, n_quad_v);
     dmat_zero(R_m_a, n_quad_v, n_quad_v);
     dmat_zero(T_m_a, n_quad_v, n_quad_v);
     dvec_zero(S_p_a, n_quad_v);
     dvec_zero(S_m_a, n_quad_v);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 4 * n_quad_v * n_quad_v + 2 * n_quad_v;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad_v * n_quad_v) {
               i_input = 0;
               ii = i;
               ii_input = ii / n_quad_v;
               jj_input = ii % n_quad_v;
          }
          else
          if (i < 2 * n_quad_v * n_quad_v) {
               i_input = 1;
               ii = i - n_quad_v * n_quad_v;
               ii_input = ii / n_quad_v;
               jj_input = ii % n_quad_v;
          }
          else
          if (i < 3 * n_quad_v * n_quad_v) {
               i_input = 2;
               ii = i - 2 * n_quad_v * n_quad_v;
               ii_input = ii / n_quad_v;
               jj_input = ii % n_quad_v;
          }
          else
          if (i < 4 * n_quad_v * n_quad_v) {
               i_input = 3;
               ii = i - 3 * n_quad_v * n_quad_v;
               ii_input = ii / n_quad_v;
               jj_input = ii % n_quad_v;
          }
          else
          if (i < 4 * n_quad_v * n_quad_v + n_quad_v) {
               i_input = 4;
               ii = i - 4 * n_quad_v * n_quad_v;
               ii_input = ii;
          }
          else {
               i_input = 5;
               ii = i - (4 * n_quad_v * n_quad_v + n_quad_v);
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          omega_a = 0.;
          ltau_a  = 0.;
          as_0_a  = 0.;
          atran_a = 0.;

          dvec_zero(P_0p_a, n_quad_v);
          dvec_zero(P_0m_a, n_quad_v);

          dmat_zero(r_p_a, n_quad_v, n_quad_v);
          dmat_zero(t_p_a, n_quad_v, n_quad_v);
          dmat_zero(r_m_a, n_quad_v, n_quad_v);
          dmat_zero(t_m_a, n_quad_v, n_quad_v);


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


          rtm_eig_rts_a(n_quad, n_stokes, F_0,
                        qx_v, qw_v, planck0, planck1,
                        omega, &omega_a, ltau, &ltau_a,
                        as_0, &as_0_a, atran, &atran_a,
                        P_0p, P_0m,
                        r_p, t_p, r_m, t_m,
                        R_p, T_p, R_m, T_m,
                        S_p, S_m,
                        P_0p_a,  P_0m_a,
                        r_p_a, t_p_a, r_m_a, t_m_a,
                        R_p_a, T_p_a, R_m_a, T_m_a,
                        S_p_a, S_m_a,
                        symmetric, thermal, vector,
                        eigen_solver_real, eigen_solver_complex,
                        derivs_layers[0], derivs_beam_down[0], save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 4 || i_input == 5) {
               if (i_input == 4)
                    a = S_p_l;
               else
                    a = S_m_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs_layers[j] && ! derivs_beam_down[j])
                         continue;

                    a[j][ii_input] = 0.;

                    if (derivs_layers[j]) {
                         a[j][ii_input] += omega_l[j] * omega_a;
                         a[j][ii_input] += ltau_l [j] * ltau_a;

                         for (k = 0; k < n_quad_v; ++k) {
                              a[j][ii_input] += P_0p_l[j][k] * P_0p_a[k];
                              a[j][ii_input] += P_0m_l[j][k] * P_0m_a[k];

                              for (l = 0; l < n_quad_v; ++l) {
                                   a[j][ii_input] += r_p_l[j][k][l] * r_p_a[k][l];
                                   a[j][ii_input] += t_p_l[j][k][l] * t_p_a[k][l];
                                   a[j][ii_input] += r_m_l[j][k][l] * r_m_a[k][l];
                                   a[j][ii_input] += t_m_l[j][k][l] * t_m_a[k][l];
                              }
                         }
                    }

                    if (derivs_beam_down[j]) {
                         a[j][ii_input] += as_0_l [j] * as_0_a;
                         a[j][ii_input] += atran_l[j] * atran_a;
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

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs_layers[j])
                         continue;

                    b[j][ii_input][jj_input] = 0.;

                    b[j][ii_input][jj_input] += omega_l[j] * omega_a;
                    b[j][ii_input][jj_input] += ltau_l [j] * ltau_a;
                    b[j][ii_input][jj_input] += as_0_l [j] * as_0_a;
                    b[j][ii_input][jj_input] += atran_l[j] * atran_a;

                    for (k = 0; k < n_quad_v; ++k) {
                         b[j][ii_input][jj_input] += P_0p_l[j][k] * P_0p_a[k];
                         b[j][ii_input][jj_input] += P_0m_l[j][k] * P_0m_a[k];

                         for (l = 0; l < n_quad_v; ++l) {
                              b[j][ii_input][jj_input] += r_p_l[j][k][l] * r_p_a[k][l];
                              b[j][ii_input][jj_input] += t_p_l[j][k][l] * t_p_a[k][l];
                              b[j][ii_input][jj_input] += r_m_l[j][k][l] * r_m_a[k][l];
                              b[j][ii_input][jj_input] += t_m_l[j][k][l] * t_m_a[k][l];
                         }
                    }
               }
          }
     }
}
