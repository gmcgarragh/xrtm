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
void radiance_slab_tl_with_ad(int n_quad, int n_derivs,
                              double **R_m, double **T_p, double *S_p,
                              double ***R_m_l, double ***T_p_l, double **S_p_l,
                              double *I1_m, double *I2_p, double **I1_m_l,
                              double **I2_p_l, double *I1_p, double **I1_p_l,
                              uchar *derivs, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;

     double **R_m_a;
     double **T_p_a;
     double *S_p_a;
     double *I1_m_a;
     double *I2_p_a;

     double *I1_p_a;

     double **a;


     if (! flags_or(derivs, n_derivs))
          return;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R_m_a  = get_work1(&work, WORK_DXX);
     T_p_a  = get_work1(&work, WORK_DXX);
     S_p_a  = get_work1(&work, WORK_DX);
     I1_m_a = get_work1(&work, WORK_DX);
     I2_p_a = get_work1(&work, WORK_DX);

     I1_p_a = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(I1_p_a, n_quad, 0.);

     n_inputs = n_quad;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          dmat_zero(R_m_a,  n_quad, n_quad);
          dmat_zero(T_p_a,  n_quad, n_quad);
          dvec_zero(S_p_a,  n_quad);
          init_array1_d(I1_m_a, n_quad, 0.);
          init_array1_d(I2_p_a, n_quad, 0.);


          if (i_input == 0)
               I1_p_a[ii_input] = 1.;


          radiance_slab_a(n_quad, R_m, T_p, S_p, R_m_a, T_p_a, S_p_a, I1_m, I2_p, I1_m_a, I2_p_a, I1_p, I1_p_a, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0) {
               if (i_input == 0)
                    a = I1_p_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    a[j][ii_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         a[j][ii_input] += I1_m_l[j][k] * I1_m_a[k];
                         a[j][ii_input] += I2_p_l[j][k] * I2_p_a[k];
                         a[j][ii_input] += S_p_l [j][k] * S_p_a [k];

                         for (l = 0; l < n_quad; ++l) {
                              a[j][ii_input] += R_m_l[j][k][l] * R_m_a[k][l];
                              a[j][ii_input] += T_p_l[j][k][l] * T_p_a[k][l];
                         }
                    }
               }
          }
     }
}
