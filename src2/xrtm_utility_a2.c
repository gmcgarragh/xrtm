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
void build_local_r_and_t_tl_with_ad(int i_four, int n_quad, int n_derivs,
                                    double *qx_v, double *qw_v, double omega, double *omega_l,
                                    double  **P_p, double  **P_m, double  **r, double  **t,
                                    double ***P_p_l, double ***P_m_l, double ***r_l, double ***t_l,
                                    uchar *derivs, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double omega_a;

     double **P_p_a;
     double **P_m_a;

     double **r_a;
     double **t_a;

     double ***a;


     if (! flags_or(derivs, n_derivs))
          return;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     P_p_a = get_work1(&work, WORK_DXX);
     P_m_a = get_work1(&work, WORK_DXX);

     r_a   = get_work1(&work, WORK_DXX);
     t_a   = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(r_a, n_quad, n_quad);
     dmat_zero(t_a, n_quad, n_quad);


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
          omega_a = 0.;
          dmat_zero(P_p_a, n_quad, n_quad);
          dmat_zero(P_m_a, n_quad, n_quad);


          if (i_input == 0)
               r_a[ii_input][jj_input] = 1.;
          else
               t_a[ii_input][jj_input] = 1.;


          build_local_r_and_t_a(i_four, n_quad, qx_v, qw_v, omega, &omega_a, P_p, P_m, r, t, P_p_a, P_m_a, r_a, t_a, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = r_l;
               else
                    a = t_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    a[j][ii_input][jj_input] = 0.;

                    a[j][ii_input][jj_input] += omega_l[j] * omega_a;

                    for (k = 0; k < n_quad; ++k) {
                         for (l = 0; l < n_quad; ++l) {
                              a[j][ii_input][jj_input] += P_p_l[j][k][l] * P_p_a[k][l];
                              a[j][ii_input][jj_input] += P_m_l[j][k][l] * P_m_a[k][l];
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_txr_tl_with_ad(int n_quad, int n_stokes, int n_derivs,
                          double **r_p, double **t_p, double **tpr, double **tmr,
                          double ***r_p_l, double ***t_p_l, double ***tpr_l, double ***tmr_l,
                          uchar *derivs, work_data work) {

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

     double **r_p_a;
     double **t_p_a;

     double **tpr_a;
     double **tmr_a;

     double ***a;


     if (! flags_or(derivs, n_derivs))
          return;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     r_p_a = get_work1(&work, WORK_DXX);
     t_p_a = get_work1(&work, WORK_DXX);

     tpr_a = get_work1(&work, WORK_DXX);
     tmr_a = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(tpr_a, n_quad_v, n_quad_v);
     dmat_zero(tmr_a, n_quad_v, n_quad_v);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * n_quad_v * n_quad_v;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad_v * n_quad_v) {
               i_input = 0;
               ii = i;
               ii_input = ii / n_quad_v;
               jj_input = ii % n_quad_v;
          }
          else {
               i_input = 1;
               ii = i - n_quad_v * n_quad_v;
               ii_input = ii / n_quad_v;
               jj_input = ii % n_quad_v;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          dmat_zero(r_p_a, n_quad_v, n_quad_v);
          dmat_zero(t_p_a, n_quad_v, n_quad_v);


          if (i_input == 0)
               tpr_a[ii_input][jj_input] = 1.;
          else
               tmr_a[ii_input][jj_input] = 1.;


          build_txr_a(n_quad, n_stokes, r_p, t_p, tpr, tmr, r_p_a, t_p_a, tpr_a, tmr_a, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = tpr_l;
               else
                    a = tmr_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    a[j][ii_input][jj_input] = 0.;

                    for (k = 0; k < n_quad_v; ++k) {
                         for (l = 0; l < n_quad_v; ++l) {
                              a[j][ii_input][jj_input] += r_p_l[j][k][l] * r_p_a[k][l];
                              a[j][ii_input][jj_input] += t_p_l[j][k][l] * t_p_a[k][l];
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_gamma_tl_with_ad(int n_quad, int n_derivs,
                            double **tpr, double **tmr, double **gamma,
                            double ***tpr_l, double ***tmr_l, double ***gamma_l,
                            uchar *derivs, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double **tpr_a;
     double **tmr_a;

     double **gamma_a;

     double ***a;


     if (! flags_or(derivs, n_derivs))
          return;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     tpr_a = get_work1(&work, WORK_DXX);
     tmr_a = get_work1(&work, WORK_DXX);

     gamma_a = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(gamma_a, n_quad, n_quad);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = n_quad * n_quad;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad * n_quad) {
               i_input = 0;
               ii = i;
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          dmat_zero(tpr_a, n_quad, n_quad);
          dmat_zero(tmr_a, n_quad, n_quad);


          if (i_input == 0)
               gamma_a[ii_input][jj_input] = 1.;


          build_gamma_a(n_quad, tpr, tmr, gamma, tpr_a, tmr_a, gamma_a, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0) {
               if (i_input == 0)
                    a = gamma_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    a[j][ii_input][jj_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         for (l = 0; l < n_quad; ++l) {
                              a[j][ii_input][jj_input] += tpr_l[j][k][l] * tpr_a[k][l];
                              a[j][ii_input][jj_input] += tmr_l[j][k][l] * tmr_a[k][l];
                         }
                    }
               }
          }
     }
}
