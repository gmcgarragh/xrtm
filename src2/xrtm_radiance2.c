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
static void internal_radiance(int n_quad, int n_derivs,
          double  **R12_p, double  **T12_m, double  **R23_m, double  **T23_p,
          double  *S12_m, double  *S23_p,
          double ***R12_p_l, double ***T12_m_l, double ***R23_m_l, double ***T23_p_l,
          double **S12_m_l, double **S23_p_l,
          double *I1_m, double *I3_p, double **I1_m_l, double **I3_p_l,
          double *I2_p, double **I2_p_l, int derivs, work_data work) {

     int i;

     int *i1;

     double *v1;
     double *v2;
     double *v3;

     double **w1;

     i1 = get_work1(&work, WORK_IX);

     v1 = get_work1(&work, WORK_DX);
     v2 = get_work1(&work, WORK_DX);

     w1 = get_work1(&work, WORK_DXX);

     dmat_mul(R23_m, R12_p, n_quad, n_quad, n_quad, w1);
     dmat_i_sub(w1, w1, n_quad);
     dmat_getrf(w1, n_quad, n_quad, i1);

     dm_v_mul(T12_m, I1_m, n_quad, n_quad, v1);
     dm_v_mul(R23_m, v1, n_quad, n_quad, v2);
     dm_v_mul(T23_p, I3_p, n_quad, n_quad, v1);
     dvec_add(v2, v1, v1, n_quad);
     dm_v_mul(R23_m, S12_m, n_quad, n_quad, v2);
     dvec_add(v1, v2, v1, n_quad);
     dvec_add(v1, S23_p, I2_p, n_quad);

     dmat_getrs(w1, &I2_p, n_quad, 1, i1);

     if (n_derivs > 0) {
          v3 = get_work1(&work, WORK_DX);

          for (i = 0; i < n_derivs; ++i) {
               dm_v_mul(R12_p, I2_p, n_quad, n_quad, v1);
               dm_v_mul(R23_m_l[i], v1, n_quad, n_quad, v2);
               dm_v_mul(R12_p_l[i], I2_p, n_quad, n_quad, v1);
               dm_v_mul(R23_m, v1, n_quad, n_quad, v3);
               dvec_add(v2, v3, v1, n_quad);

               dm_v_mul(T12_m, I1_m, n_quad, n_quad, v2);
               dm_v_mul(R23_m_l[i], v2, n_quad, n_quad, v3);
               dvec_add(v1, v3, v1, n_quad);

               dm_v_mul(T12_m_l[i], I1_m, n_quad, n_quad, v2);
               dm_v_mul(T12_m, I1_m_l[i], n_quad, n_quad, v3);
               dvec_add(v2, v3, v2, n_quad);
               dm_v_mul(R23_m, v2, n_quad, n_quad, v3);
               dvec_add(v1, v3, v1, n_quad);

               dm_v_mul(T23_p_l[i], I3_p, n_quad, n_quad, v2);
               dvec_add(v1, v2, v1, n_quad);

               dm_v_mul(T23_p, I3_p_l[i], n_quad, n_quad, v2);
               dvec_add(v1, v2, v1, n_quad);

               dm_v_mul(R23_m_l[i], S12_m, n_quad, n_quad, v2);
               dvec_add(v1, v2, v1, n_quad);

               dm_v_mul(R23_m, S12_m_l[i], n_quad, n_quad, v2);
               dvec_add(v1, v2, v1, n_quad);

               dvec_add(v1, S23_p_l[i], I2_p_l[i], n_quad);

               dmat_getrs(w1, &I2_p_l[i], n_quad, 1, i1);
          }
     }
}



void radiance_internal(int n_quad, int n_derivs,
   double  **R12_p, double  **T12_m, double  **R23_m, double  **T23_p,
   double  *S12_m, double  *S23_p,
   double ***R12_p_l, double ***T12_m_l, double ***R23_m_l, double ***T23_p_l,
   double **S12_m_l, double **S23_p_l,
   double *I1_m, double *I3_p, double **I1_m_l, double **I3_p_l,
   double *I2_p, double *I2_m, double **I2_p_l, double **I2_m_l,
   int derivs, work_data work) {

     internal_radiance(n_quad, n_derivs,
                       R12_p, T12_m, R23_m, T23_p, S12_m, S23_p,
                       R12_p_l, T12_m_l, R23_m_l, T23_p_l, S12_m_l, S23_p_l,
                       I1_m, I3_p, I1_m_l, I3_p_l, I2_p, I2_p_l, derivs, work);
     internal_radiance(n_quad, n_derivs,
                       R23_m, T23_p, R12_p, T12_m, S23_p, S12_m,
                       R23_m_l, T23_p_l, R12_p_l, T12_m_l, S23_p_l, S12_m_l,
                       I3_p, I1_m, I3_p_l, I1_m_l, I2_m, I2_m_l, derivs, work);
}
