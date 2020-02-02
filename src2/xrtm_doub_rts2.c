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
void rtm_doub_rts2(int n_quad, int n_stokes, int n_derivs, double F_0,
                   double *qx_v, double *qwv, double planck0, double planck1,
                   double omega, double *omega_l, double ltau, double *ltau_l,
                   double btran, double *btran_l, double as_0, double *as_0_l,
                   double *P_q0_mm, double *P_q0_pm,
                   double **r_p, double **t_p, double **r_m, double **t_m,
                   double **R_p, double **T_p, double **R_m, double **T_m,
                   double *S_p, double *S_m,
                   double **P_q0_mm_l, double **P_q0_pm_l,
                   double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l,
                   double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l,
                   double **S_p_l, double **S_m_l,
                   double d_tau, int thermal, uchar *derivs_layers, uchar *derivs_beam_down,
                   save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int jj;

     int n_quad_v;

     int n_doub;

     double a;
     double b;
     double c;
     double d;

     double d_tran;

     double lin_fac;

     double *v1;
     double *v2;

     double *d_tau_l;

     double *d_path_l;
     double *d_tran_l;

     double *Sl_p;
     double *Sl_m;

     double **w1;

     double **Sl_p_l;
     double **Sl_m_l;


     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1       = get_work1(&work, WORK_DX);

     if (thermal) {
          v2   = get_work1(&work, WORK_DX);

          Sl_p = get_work1(&work, WORK_DX);
          Sl_m = get_work1(&work, WORK_DX);
     }

     d_tau_l  = get_work1(&work, WORK_DDERIVS);

     d_path_l = get_work1(&work, WORK_DDERIVS);
     d_tran_l = get_work1(&work, WORK_DDERIVS);

     if (flags_or(derivs_layers, n_derivs)) {
          if (thermal) {
               Sl_p_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_layers);
               Sl_m_l = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_layers);
          }

          w1 = get_work1(&work, WORK_DXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     calc_doub_params(n_derivs, ltau, ltau_l, &n_doub, &d_tau, d_tau_l, derivs_layers);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_scale(-d_tau, r_p, R_p, n_quad_v, n_quad_v);
     dmat_scale(-d_tau, r_m, R_m, n_quad_v, n_quad_v);

     dmat_scale(-d_tau, t_p, T_p, n_quad_v, n_quad_v);
     dmat_i_sub(T_p, T_p, n_quad_v);
     dmat_scale(-d_tau, t_m, T_m, n_quad_v, n_quad_v);
     dmat_i_sub(T_m, T_m, n_quad_v);

     if (flags_or(derivs_layers, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_layers[i])
                    continue;

               dmat_scale(-d_tau, r_p_l[i], R_p_l[i], n_quad_v, n_quad_v);
               dmat_scale(-d_tau_l[i], r_p, w1 , n_quad_v, n_quad_v);
               dmat_add(R_p_l[i], w1, R_p_l[i], n_quad_v, n_quad_v);

               dmat_scale(-d_tau, r_m_l[i], R_m_l[i], n_quad_v, n_quad_v);
               dmat_scale(-d_tau_l[i], r_m, w1 , n_quad_v, n_quad_v);
               dmat_add(R_m_l[i], w1, R_m_l[i], n_quad_v, n_quad_v);

               dmat_scale( d_tau, t_p_l[i], T_p_l[i], n_quad_v, n_quad_v);
               dmat_scale( d_tau_l[i], t_p, w1 , n_quad_v, n_quad_v);
               dmat_add(T_p_l[i], w1, T_p_l[i], n_quad_v, n_quad_v);

               dmat_scale( d_tau, t_m_l[i], T_m_l[i], n_quad_v, n_quad_v);
               dmat_scale( d_tau_l[i], t_m, w1 , n_quad_v, n_quad_v);
               dmat_add(T_m_l[i], w1, T_m_l[i], n_quad_v, n_quad_v);
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (F_0 == 0.)
          d_tran = 0.;
     else {
          d_tran = exp(-d_tau * as_0);

          if (flags_or(derivs_layers, n_derivs)) {
               for (i = 0; i < n_derivs; ++i) {
                    if (! derivs_layers[i])
                         continue;

                    d_path_l[i] = -d_tau_l[i] * as_0 - d_tau * as_0_l[i];

                    d_tran_l[i] = d_path_l[i] * d_tran;
               }
          }

          a = F_0 / (4. * PI);

          b = a * omega * d_tran * d_tau;

          dm_v_dinv_mul(qx_v, P_q0_pm, S_p, n_quad_v);
          dvec_scale(b, S_p, S_p, n_quad_v);

          dm_v_dinv_mul(qx_v, P_q0_mm, S_m, n_quad_v);
          dvec_scale(b, S_m, S_m, n_quad_v);
     }

     if (! thermal)
          lin_fac = 0.;
     else {
          if (planck0 == 0.)
               lin_fac = 0.;
          else
               lin_fac = (planck1 / planck0 - 1.) / pow(2., n_doub);
          dvec_zero(v1, n_quad_v);

          ii = 0;
          c = (1. - omega) * planck0 * d_tau;
          for (i = 0; i < n_quad; ++i) {
               v1[ii] = c;
               ii += n_stokes;
          }
          dm_v_dinv_mul(qx_v, v1, v1, n_quad_v);
          dvec_copy(Sl_p, v1, n_quad_v);
          dvec_copy(Sl_m, v1, n_quad_v);
     }


     if (flags_or(derivs_layers, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_layers[i])
                    continue;

               if (F_0 != 0.) {
                    c = a * (omega_l[i] * d_tran      * d_tau +
                             omega      * d_tran_l[i] * d_tau +
                             omega      * d_tran      * d_tau_l[i]);
/*
                    c = a * (omega_l[i] * btran      * d_tran      * d_tau +
                             omega      * btran_l[i] * d_tran      * d_tau +
                             omega      * btran      * d_tran_l[i] * d_tau +
                             omega      * btran      * d_tran      * d_tau_l[i]);
*/
                    dvec_scale(c, P_q0_pm, S_p_l[i], n_quad_v);
                    dvec_scale(b, P_q0_pm_l[i], v1, n_quad_v);
                    dvec_add(S_p_l[i], v1, S_p_l[i], n_quad_v);
                    dm_v_dinv_mul(qx_v, S_p_l[i], S_p_l[i], n_quad_v);

                    dvec_scale(c, P_q0_mm, S_m_l[i], n_quad_v);
                    dvec_scale(b, P_q0_mm_l[i], v1, n_quad_v);
                    dvec_add(S_m_l[i], v1, S_m_l[i], n_quad_v);
                    dm_v_dinv_mul(qx_v, S_m_l[i], S_m_l[i], n_quad_v);
               }

               if (thermal) {
                    dvec_zero(v1, n_quad_v);

                    jj = 0;
                    d = -omega_l[i] * d_tau + (1. - omega) * d_tau_l[i];
                    for (j = 0; j < n_quad; ++j) {
                         v1[jj] = d;
                         jj += n_stokes;
                    }
                    dm_v_dinv_mul(qx_v, v1, v1, n_quad_v);
                    dvec_copy(Sl_p_l[i], v1, n_quad_v);
                    dvec_copy(Sl_m_l[i], v1, n_quad_v);
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     c = 1.;
     for (i = 0; i < n_doub; ++i) {
          if (! flags_or(derivs_layers, n_derivs))
               variant_double  (R_m, T_m, S_m, Sl_m, v1, R_p, T_p, S_p, Sl_p, v2, n_quad, n_stokes, d_tran, lin_fac, F_0 != 0., thermal, i == 0, work);
          else
               variant_double_l(R_m, T_m, S_m, R_p, T_p, S_p, R_m_l, T_m_l, S_m_l, R_p_l, T_p_l, S_p_l, n_quad, n_stokes, n_derivs, d_tran, d_tran_l, derivs_layers, work);

          if (F_0 != 0.) {
               d_tran *= d_tran;

               if (flags_or(derivs_layers, n_derivs)) {
                    c *= 2.;

                    for (j = 0; j < n_derivs; ++j) {
                         if (! derivs_layers[j])
                              continue;

                         d_tran_l[j]  = d_path_l[j] * d_tran * c;
                    }
               }
          }

          if (thermal)
               lin_fac *= 2.;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (F_0 == 0. && ! thermal) {
          dvec_zero(S_p, n_quad_v);
          dvec_zero(S_m, n_quad_v);
     }
     else
     if (F_0 == 0. && thermal) {
          dvec_copy(S_p, Sl_p, n_quad_v);
          dvec_copy(S_m, Sl_m, n_quad_v);
     }
     else
     if (F_0 != 0. && thermal) {
          dvec_add (S_p, Sl_p, S_p, n_quad_v);
          dvec_add (S_m, Sl_m, S_m, n_quad_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (0) {
     double atran;
     double *atran_l;

     double *F_p;
     double *F_m;

     double **F_p_l;
     double **F_m_l;

     F_p = get_work1(&work, WORK_DX);
     F_m = get_work1(&work, WORK_DX);

     if (flags_or(derivs_beam_down, n_derivs)) {
          atran_l = get_work1(&work, WORK_DDERIVS);

          F_p_l     = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam_down);
          F_m_l     = get_work2(&work, WORK_DX, WORK_DERIVS_V, derivs_beam_down);
     }
if (1) {
     double **tpr;
     double **tmr;
     double **gamma;
     double ***tpr_l;
     double ***tmr_l;
     double ***gamma_l;

     tpr   = get_work1(&work, WORK_DXX);
     tmr   = get_work1(&work, WORK_DXX);
     gamma = get_work1(&work, WORK_DXX);

     if (flags_or(derivs_layers, n_derivs)) {
          tpr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          tmr_l   = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
          gamma_l = get_work2(&work, WORK_DXX, WORK_DERIVS_V, derivs_layers);
     }

     build_txr(n_quad, n_stokes, n_derivs, r_p, t_p, tpr, tmr, r_p_l, t_p_l, tpr_l, tmr_l, derivs_layers, work);
     build_gamma(n_quad_v, n_derivs, tpr, tmr, gamma, tpr_l, tmr_l, gamma_l, derivs_layers, work);

     build_source_vectors_solar_classic_1n (n_quad, n_stokes, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_q0_mm, P_q0_pm, tpr, tmr, gamma, F_p, F_m, P_q0_mm_l, P_q0_pm_l, tpr_l, tmr_l, gamma_l, F_p_l, F_m_l, derivs_layers, derivs_beam_down, save_tree, work);
}
else
     build_source_vectors_solar_classic_2n2(n_quad_v, n_derivs, qx_v, F_0, omega, omega_l, as_0, as_0_l, P_q0_mm, P_q0_pm, r_p, t_p, r_m, t_m, F_p, F_m, P_q0_mm_l, P_q0_pm_l, r_p_l, t_p_l, r_m_l, t_m_l, F_p_l, F_m_l, derivs_layers, derivs_beam_down, work);

     atran = exp(-ltau * as_0);

     if (flags_or(derivs_beam_down, n_derivs)) {
          for (i = 0; i < n_derivs; ++i) {
               if (! derivs_beam_down[i])
                    continue;

               atran_l[i] = -(ltau_l[i] * as_0 + ltau * as_0_l[i]) * atran;
          }
     }

     build_global_source_solar2(n_quad_v, n_derivs, atran, atran_l, R_p, T_p, R_m, T_m, F_p, F_m, S_p, S_m, R_p_l, T_p_l, R_m_l, T_m_l, F_p_l, F_m_l, S_p_l, S_m_l, derivs_layers, derivs_beam_down, work, NULL, NULL, NULL, NULL);
}

}
