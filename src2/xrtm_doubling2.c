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
void variant_double(double **R_m, double **T_m, double *Se_m, double *Sl_m, double *a_m,
                    double **R_p, double **T_p, double *Se_p, double *Sl_p, double *a_p,
                    int n_quad, int n_stokes, double atran, double lin_fac,
                    int flag1, int flag2, int flag3, work_data work) {

     layer_add2(R_m, T_m, Se_m, Sl_m, a_m, R_p, T_p, Se_p, Sl_p, a_p,
                R_m, T_m, Se_m, Sl_m, a_m, R_p, T_p, Se_p, Sl_p, a_p,
                NULL, NULL, Se_m, Sl_m, a_m, R_p, T_p, Se_p, Sl_p, a_p,
                n_quad * n_stokes, atran, lin_fac, flag1, flag2, flag3, work);

     phase_matrix_symmetry2(n_quad, n_stokes, T_p, R_p, T_m, R_m, 1.);
}



/*******************************************************************************
 *
 ******************************************************************************/
void variant_double_l(double **R_m, double **T_m, double *S_m,
                      double **R_p, double **T_p, double *S_p,
                      double ***R_m_l, double ***T_m_l, double **S_m_l,
                      double ***R_p_l, double ***T_p_l, double **S_p_l,
                      int n_quad, int n_stokes, int n_derivs,
                      double atran, double *atran_l,
                      uchar *derivs, work_data work) {

     int i;

     layer_add2_l(R_m, T_m, S_m, R_p, T_p, S_p,
                  R_m, T_m, S_m, R_p, T_p, S_p,
                  R_m, T_m, S_m, NULL, NULL, S_p,
                  R_m_l, T_m_l, S_m_l, R_p_l, T_p_l, S_p_l,
                  R_m_l, T_m_l, S_m_l, R_p_l, T_p_l, S_p_l,
                  R_m_l, T_m_l, S_m_l, NULL, NULL, S_p_l,
                  n_quad * n_stokes, n_derivs, atran, atran_l, derivs, work);

     phase_matrix_symmetry2(n_quad, n_stokes, T_m, R_m, T_p, R_p, 1.);
     for (i = 0; i < n_derivs; ++i) {
          if (! derivs[i])
               continue;

          phase_matrix_symmetry2(n_quad, n_stokes, T_m_l[i], R_m_l[i], T_p_l[i], R_p_l[i], 1.);
     }
}
