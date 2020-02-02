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
void build_phase_mats_scalar_tl_with_ad(int i_four, int n_coef, int n_mus1, int n_mus2, double **Y1, double **Y2, int lda, double *coefs_l, double **P_pp_l, double **P_mp_l, work_data work) {

     int i;
     int ii;
     int j;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double *coefs_a;

     double **P_pp_a;
     double **P_mp_a;

     double **a;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     coefs_a  = get_work_d1(&work, n_coef);

     P_pp_a = get_work1(&work, WORK_DXX);
     P_mp_a = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(P_pp_a, n_mus1, n_mus2);
     dmat_zero(P_mp_a, n_mus1, n_mus2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * n_mus1 * n_mus2;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_mus1 * n_mus2) {
               i_input = 0;
               ii = i;
               ii_input = ii / n_mus2;
               jj_input = ii % n_mus2;
          }
          else {
               i_input = 1;
               ii = i - (n_mus1 * n_mus2);
               ii_input = ii / n_mus2;
               jj_input = ii % n_mus2;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array1_d(coefs_a, n_coef, 0.);


          if (i_input == 0)
               P_pp_a[ii_input][jj_input] = 1.;
          else
               P_mp_a[ii_input][jj_input] = 1.;


          build_phase_mats_scalar_a(i_four, n_coef, n_mus1, n_mus2, Y1, Y2, lda, coefs_a, P_pp_a, P_mp_a);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = P_pp_l;
               else
                    a = P_mp_l;

               a[ii_input][jj_input] = 0.;

               for (j = i_four; j < n_coef; ++j)
                    a[ii_input][jj_input] += coefs_l[j] * coefs_a[j];
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_vecs_vector_gc_tl_with_ad(int i_four, int n_coef, int n_mus, int n_stokes, double *mu, double mu_0, double **gc, double *P_pp, double *P_mp, work_data work) {

     int i;
     int ii;
     int j;
     int k;

     int n_mus_v;

     int i_input;
     int n_inputs;

     int ii_input;

     double **gc_a;

     double *P_pp_a;
     double *P_mp_a;

     double *a;


     n_mus_v = n_mus * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     gc_a  = get_work_d2(&work, 6, n_coef);

     P_pp_a = get_work1(&work, WORK_DX);
     P_mp_a = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(P_pp_a, n_mus_v);
     dvec_zero(P_mp_a, n_mus_v);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * n_mus_v;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_mus_v) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }
          else {
               i_input = 1;
               ii = i - n_mus_v;
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array2_d(gc_a, 6, n_coef, 0.);


          if (i_input == 0)
               P_pp_a[ii_input] = 1.;
          else
               P_mp_a[ii_input] = 1.;


          build_phase_vecs_vector_gc_a(i_four, n_coef, n_mus, n_stokes, mu, mu_0, gc_a, P_pp_a, P_mp_a, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = P_pp;
               else
                    a = P_mp;

               a[ii_input] = 0.;

               for (j = 0; j < 6; ++j) {
                    for (k = i_four; k < n_coef; ++k) {
                         a[ii_input] += gc[j][k] * gc_a[j][k];
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void build_phase_mats_vector_gc_tl_with_ad(int i_four, int n_coef, int n_mus1, int n_mus2, int n_stokes, double *mu1, double *mu2, double **gc, double **P_pp, double **P_mp, double **P_mm, double **P_pm, int vector, work_data work) {

     int i;
     int ii;
     int j;
     int k;

     int n_mus_v1;
     int n_mus_v2;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double **gc_a;

     double **P_pp_a;
     double **P_mp_a;
     double **P_mm_a;
     double **P_pm_a;

     double **a;


     n_mus_v1 = n_mus1 * n_stokes;
     n_mus_v2 = n_mus2 * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     gc_a  = get_work_d2(&work, 6, n_coef);

     P_pp_a = get_work1(&work, WORK_DXX);
     P_mp_a = get_work1(&work, WORK_DXX);
     P_mm_a = get_work1(&work, WORK_DXX);
     P_pm_a = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(P_pp_a, n_mus_v1, n_mus_v2);
     dmat_zero(P_mp_a, n_mus_v1, n_mus_v2);
     dmat_zero(P_mm_a, n_mus_v1, n_mus_v2);
     dmat_zero(P_pm_a, n_mus_v1, n_mus_v2);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 4 * n_mus_v1 * n_mus_v2;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_mus_v1 * n_mus_v2) {
               i_input = 0;
               ii = i;
               ii_input = ii / n_mus_v2;
               jj_input = ii % n_mus_v2;
          }
          else
          if (i < 2 * n_mus_v1 * n_mus_v2) {
               i_input = 1;
               ii = i - (n_mus_v1 * n_mus_v2);
               ii_input = ii / n_mus_v2;
               jj_input = ii % n_mus_v2;
          }
          else
          if (i < 3 * n_mus_v1 * n_mus_v2) {
               i_input = 2;
               ii = i - (2 * n_mus_v1 * n_mus_v2);
               ii_input = ii / n_mus_v2;
               jj_input = ii % n_mus_v2;
          }
          else {
               i_input = 3;
               ii = i - (3 * n_mus_v1 * n_mus_v2);
               ii_input = ii / n_mus_v2;
               jj_input = ii % n_mus_v2;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array2_d(gc_a, 6, n_coef, 0.);


          if (i_input == 0)
               P_pp_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 1)
               P_mp_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 2)
               P_mm_a[ii_input][jj_input] = 1.;
          else
               P_pm_a[ii_input][jj_input] = 1.;


          build_phase_mats_vector_gc_a(i_four, n_coef, n_mus1, n_mus2, n_stokes, mu1, mu2, gc_a, P_pp_a, P_mp_a, P_mm_a, P_pm_a, vector, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1 || i_input == 2 || i_input == 3) {
               if (i_input == 0)
                    a = P_pp;
               else
               if (i_input == 1)
                    a = P_mp;
               else
               if (i_input == 2)
                    a = P_mm;
               else
                    a = P_pm;

               a[ii_input][jj_input] = 0.;

               for (j = 0; j < 6; ++j) {
                    for (k = i_four; k < n_coef; ++k) {
                         a[ii_input][jj_input] += gc[j][k] * gc_a[j][k];
                    }
               }
          }
     }
}
