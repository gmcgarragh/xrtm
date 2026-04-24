/*******************************************************************************
**
**    Copyright (C) 2007-2020 Greg McGarragh <greg.mcgarragh@colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


#define USE_BANDED_SOLVER


#ifdef USE_BANDED_SOLVER
#ifdef __cplusplus
extern "C" {
#endif

void xgbtrf_(int *, int *, int *, int *, TYPE *, int *, int *, int *);
void xgbtrs_(const char *, int *, int *, int *, int *, TYPE *, int *, int *, TYPE *, int *, int *);

#ifdef __cplusplus
}
#endif
#endif


/*******************************************************************************
 *
 ******************************************************************************/
static void SOLVE_BVP(int n_quad, int n_stokes, int n_derivs, int n_layers,
                      double qf, double *qx_v, double *qw_v, double F_0,
                      int n_ulevels, int *ulevels,
                      double *omega, double **omega_l, double *ltau, double **ltau_l,
                      double *btau, double **btau_l, double *btran, double **btran_l,
                      double *atran, double **atran_l,
                      double **P_q0_mm, double **P_q0_pm,
                      TYPE  **nu, TYPE ***X_p, TYPE ***X_m,
                      double  **Fs_p, double  **Fs_m,
                      double ***P_q0_mm_l, double ***P_q0_pm_l,
                      TYPE ***nu_l, TYPE ****X_p_l, TYPE ****X_m_l,
                      double ***Fs_p_l, double ***Fs_m_l,
                      double **Rs_qq, double ***Rs_qq_l,
                      TYPE *B, TYPE **B_l,
                      double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
                      int surface,
                      uchar **derivs_layers, uchar **derivs_beam, work_data work,
                      TYPE ***X_i_11, TYPE ***X_i_12, TYPE **sigma_p, TYPE **sigma_m,
                      TYPE ****X_i_11_l, TYPE ****X_i_12_l, TYPE ***sigma_p_l, TYPE ***sigma_m_l) {

     int i;
     int ii;
     int iii;
     int j;
     int jj;
     int m;

     int n_quad_v;
     int n_quad_v2;
     int n_quad_v3;
#ifdef USE_BANDED_SOLVER
     int n_diags;
#endif
     int m_comp;
     int n_comp;
#ifdef USE_BANDED_SOLVER
     int info;
     int nrhs = 1;
#endif
     int *ipiv;

     int *i1;
     int *i2;

     TYPE aa;
     TYPE bb;
     TYPE bb_l;
     TYPE cc;
     TYPE cc_l;

     TYPE *v1;
     TYPE *v2;
     TYPE *v3;
     TYPE *v4;

     TYPE **alpha;
     TYPE **beta1;
     TYPE **beta2;

     TYPE **w1;
     TYPE **w2;
     TYPE **w3;

     TYPE ***a_inv;
     TYPE ***b_inv;

     TYPE **A;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_quad_v  = n_quad * n_stokes;
     n_quad_v2 = n_quad_v * 2;
     n_quad_v3 = n_quad_v * 3;

#ifndef USE_BANDED_SOLVER
     m_comp    = 2 * n_quad_v * (n_layers + 1);
#else
     n_diags   = 3 * n_quad_v - 1;
     m_comp    = 3 * n_diags + 1;
#endif
     n_comp    = 2 * n_quad_v * (n_layers + 1);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i1    = get_work1(&work, WORK_IX);
     i2    = get_work1(&work, WORK_IX);

     v1    = get_work1(&work, WORK_XX);
     v2    = get_work1(&work, WORK_XX);

     alpha = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);
     beta1 = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);
     beta2 = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);

     w1    = get_work1(&work, WORK_XXX);

     a_inv = get_work2(&work, WORK_XXX, WORK_LAYERS_V, NULL);
     b_inv = get_work2(&work, WORK_XXX, WORK_LAYERS_V, NULL);

     ipiv  = get_work_i1(&work, n_comp);

     A     = get_work_x2(&work, n_comp, m_comp);

     if (flags_or2(derivs_layers, n_layers, n_derivs)) {
          v3    = get_work1(&work, WORK_XX);
          v4    = get_work1(&work, WORK_XX);

          w2    = get_work1(&work, WORK_XXX);
          w3    = get_work1(&work, WORK_XXX);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     xmat_zero(A, n_comp, m_comp);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     xmat_init(w1, 1., n_quad_v, n_quad_v);
#ifndef USE_BANDED_SOLVER
     insert_matrix1_x(n_quad_v, n_quad_v, 0, 1, A, 1., w1);
#else
     insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, 0, 1, A, 1., w1);
#endif

     ii = 2 * n_layers;
#ifndef USE_BANDED_SOLVER
     insert_matrix1_x(n_quad_v, n_quad_v, ii + 1, ii, A, 1., w1);
#else
     insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, ii, A, 1., w1);
#endif
     if (surface) {
          for (i = 0; i < n_quad_v; ++i) {
               for (j = 0; j < n_quad_v; ++j) {
                    w1[i][j] = Rs_qq[i][j];
               }
          }
#ifndef USE_BANDED_SOLVER
          insert_matrix1_x(n_quad_v, n_quad_v, ii + 1, ii + 1, A, -1., w1);
#else
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, ii + 1, A, -1., w1);
#endif
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_quad_v; ++i)
          B[i] = 0.;


     for (i = 0; i < n_quad_v; ++i)
          B[n_comp - n_quad_v + i] = In_p[i];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     ii = 1;
     jj = 0;

     for (i = 0; i < n_layers; ++i) {

          for (j = 0; j < n_quad_v; ++j)
               alpha[i][j] = XEXP(-nu[i][j] * ltau[i]);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          xmat_sub(X_p[i], X_m[i], a_inv[i], n_quad_v, n_quad_v);
          xmat_getrf(a_inv[i], n_quad_v, n_quad_v, i1);
          xmat_getri(a_inv[i], n_quad_v, i1);

          xmat_add(X_p[i], X_m[i], b_inv[i], n_quad_v, n_quad_v);
          xmat_getrf(b_inv[i], n_quad_v, n_quad_v, i2);
          xmat_getri(b_inv[i], n_quad_v, i2);

          xmat_add(a_inv[i], b_inv[i], X_i_11[i], n_quad_v, n_quad_v);
          xmat_scale(.5, X_i_11[i], X_i_11[i], n_quad_v, n_quad_v);
          xmat_sub(a_inv[i], b_inv[i], X_i_12[i], n_quad_v, n_quad_v);
          xmat_scale(.5, X_i_12[i], X_i_12[i], n_quad_v, n_quad_v);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
#ifndef USE_BANDED_SOLVER
          insert_matrix1_x(n_quad_v, n_quad_v, ii + 1, jj + 1, A, -1., X_i_11[i]);
          insert_matrix1_x(n_quad_v, n_quad_v, ii,     jj + 2, A, -1., X_i_11[i]);

          xmat_diag_mul(alpha[i], X_i_11[i], w1, n_quad_v, n_quad_v);
          insert_matrix1_x(n_quad_v, n_quad_v, ii,     jj,     A,  1., w1);
          insert_matrix1_x(n_quad_v, n_quad_v, ii + 1, jj + 3, A,  1., w1);

          insert_matrix1_x(n_quad_v, n_quad_v, ii + 1, jj,     A, -1., X_i_12[i]);
          insert_matrix1_x(n_quad_v, n_quad_v, ii,     jj + 3, A, -1., X_i_12[i]);

          xmat_diag_mul(alpha[i], X_i_12[i], w1, n_quad_v, n_quad_v);
          insert_matrix1_x(n_quad_v, n_quad_v, ii,     jj + 1, A,  1., w1);
          insert_matrix1_x(n_quad_v, n_quad_v, ii + 1, jj + 2, A,  1., w1);
#else
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, jj + 1, A, -1., X_i_11[i]);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii,     jj + 2, A, -1., X_i_11[i]);

          xmat_diag_mul(alpha[i], X_i_11[i], w1, n_quad_v, n_quad_v);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii,     jj,     A,  1., w1);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, jj + 3, A,  1., w1);

          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, jj,     A, -1., X_i_12[i]);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii,     jj + 3, A, -1., X_i_12[i]);

          xmat_diag_mul(alpha[i], X_i_12[i], w1, n_quad_v, n_quad_v);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii,     jj + 1, A,  1., w1);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, jj + 2, A,  1., w1);
#endif

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (j = 0; j < n_quad_v; ++j) {
               bb = nu[i][j] * ltau[i];
               beta1[i][j] = (btran[i  ] * alpha[i][j] - btran[i+1]) / (btau[i+1] - btau[i  ] - bb);
               beta2[i][j] = (btran[i+1] * alpha[i][j] - btran[i  ]) / (btau[i  ] - btau[i+1] - bb);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          aa = F_0 / (4. * PI);

          bb = aa * omega[i] * ltau[i];

          for (j = 0; j < n_quad_v; ++j) {
               cc = bb / qx_v[j];
               sigma_p[i][j] =  cc * P_q0_pm[i][j];
               sigma_m[i][j] = -cc * P_q0_mm[i][j];
               if (j % n_stokes > 1)
                    sigma_m[i][j] = -sigma_m[i][j];
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          iii = ii * n_quad_v;

          xm_v_mul(X_i_11[i], sigma_p[i], n_quad_v, n_quad_v, v1);
          xm_v_mul(X_i_12[i], sigma_m[i], n_quad_v, n_quad_v, v2);
          xvec_add(v1, v2, v1, n_quad_v);

          xm_v_diag_mul(beta1[i], v1, &B[iii], n_quad_v);


          xm_v_mul(X_i_11[i], sigma_m[i], n_quad_v, n_quad_v, v1);
          xm_v_mul(X_i_12[i], sigma_p[i], n_quad_v, n_quad_v, v2);
          xvec_add(v1, v2, v1, n_quad_v);

          xm_v_diag_mul(beta2[i], v1, &B[iii + n_quad_v], n_quad_v);
          xvec_scale(-1., &B[iii + n_quad_v], &B[iii + n_quad_v], n_quad_v);


          ii += 2;
          jj += 2;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef USE_BANDED_SOLVER
     xmat_getrf(A, n_comp, n_comp, ipiv);

     xmat_getrs(A, &B, n_comp, 1, ipiv);
#else
     xgbtrf_(&n_comp, &n_comp, &n_diags, &n_diags, (TYPE *) *A, &m_comp, ipiv, &info);
     if (info) {
          fprintf(stderr, "ERROR: xgbtrf() info = %d\n", info);
          exit(1);
     }

     xgbtrs_("N", &n_comp, &n_diags, &n_diags, &nrhs, (TYPE *) *A, &m_comp, ipiv,
             B, &n_comp, &info);
     if (info) {
          fprintf(stderr, "ERROR: xgbtrs() info = %d\n", info);
          exit(1);
     }
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (m = 0; m < n_derivs; ++m) {

          xvec_zero(B_l[m], n_comp);

          for (i = 0; i < n_quad_v; ++i)
               B_l[m][i] = 0.;

          for (i = 0; i < n_quad_v; ++i) {
               B_l[m][n_comp - n_quad_v + i] = In_p_l[m][i];
          }

          if (surface) {
               if (derivs_layers[n_layers][m]) {
                    for (i = 0; i < n_quad_v; ++i) {
                         for (j = 0; j < n_quad_v; ++j) {
                              w1[i][j] = Rs_qq_l[m][i][j];
                         }
                    }

                    xm_v_mul(w1, &B[n_comp - n_quad_v], n_quad_v, n_quad_v, v1);
                    xvec_add(&B_l[m][n_comp - n_quad_v], v1, &B_l[m][n_comp - n_quad_v], n_quad_v);
               }
          }


          /*---------------------------------------------------------------
           *
           *-------------------------------------------------------------*/
          ii = n_quad_v; jj = 0;

          for (i = 0; i < n_layers; ++i) {
               if (derivs_layers[i][m]) {

                    for (j = 0; j < n_quad_v; ++j)
                         v1[j] = (-nu_l[i][m][j] * ltau[i] - nu[i][j] * ltau_l[i][m]) * alpha[i][j];


                    /*-----------------------------------------------------
                     *
                     *---------------------------------------------------*/
                    xmat_sub(X_p_l[i][m], X_m_l[i][m], w1, n_quad_v, n_quad_v);
                    xmat_mul(a_inv[i], w1, n_quad_v, n_quad_v, n_quad_v, w3);
                    xmat_mul(w3, a_inv[i], n_quad_v, n_quad_v, n_quad_v, w1);

                    xmat_add(X_p_l[i][m], X_m_l[i][m], w2, n_quad_v, n_quad_v);
                    xmat_mul(b_inv[i], w2, n_quad_v, n_quad_v, n_quad_v, w3);
                    xmat_mul(w3, b_inv[i], n_quad_v, n_quad_v, n_quad_v, w2);

                    xmat_add(w1, w2, X_i_11_l[i][m], n_quad_v, n_quad_v);
                    xmat_scale(-.5, X_i_11_l[i][m], X_i_11_l[i][m], n_quad_v, n_quad_v);
                    xmat_sub(w1, w2, X_i_12_l[i][m], n_quad_v, n_quad_v);
                    xmat_scale(-.5, X_i_12_l[i][m], X_i_12_l[i][m], n_quad_v, n_quad_v);


                    /*-----------------------------------------------------
                     *
                     *---------------------------------------------------*/
                    xm_v_mul(X_i_11_l[i][m], &B[jj+n_quad_v2], n_quad_v, n_quad_v, v2);
                    for (j = 0; j < n_quad_v; ++j)
                         B_l[m][ii+j] -= -v2[j];

                    xm_v_mul(X_i_11[i], &B[jj], n_quad_v, n_quad_v, v2);
                    xm_v_diag_mul(v1, v2, v3, n_quad_v);
                    xm_v_mul(X_i_11_l[i][m], &B[jj], n_quad_v, n_quad_v, v2);
                    xm_v_diag_mul(alpha[i], v2, v4, n_quad_v);
                    xvec_add(v3, v4, v2, n_quad_v);
                    for (j = 0; j < n_quad_v; ++j)
                         B_l[m][ii+j] -=  v2[j];

                    xm_v_mul(X_i_12_l[i][m], &B[jj+n_quad_v3], n_quad_v, n_quad_v, v2);
                    for (j = 0; j < n_quad_v; ++j)
                         B_l[m][ii+j] -= -v2[j];

                    xm_v_mul(X_i_12[i], &B[jj+n_quad_v], n_quad_v, n_quad_v, v2);
                    xm_v_diag_mul(v1, v2, v3, n_quad_v);
                    xm_v_mul(X_i_12_l[i][m], &B[jj+n_quad_v], n_quad_v, n_quad_v, v2);
                    xm_v_diag_mul(alpha[i], v2, v4, n_quad_v);
                    xvec_add(v3, v4, v2, n_quad_v);
                    for (j = 0; j < n_quad_v; ++j)
                         B_l[m][ii+j] -=  v2[j];


                    xm_v_mul(X_i_11_l[i][m], &B[jj+n_quad_v], n_quad_v, n_quad_v, v2);
                    for (j = 0; j < n_quad_v; ++j)
                         B_l[m][ii+n_quad_v+j] -= -v2[j];

                    xm_v_mul(X_i_11[i], &B[jj+n_quad_v3], n_quad_v, n_quad_v, v2);
                    xm_v_diag_mul(v1, v2, v3, n_quad_v);
                    xm_v_mul(X_i_11_l[i][m], &B[jj+n_quad_v3], n_quad_v, n_quad_v, v2);
                    xm_v_diag_mul(alpha[i], v2, v4, n_quad_v);
                    xvec_add(v3, v4, v2, n_quad_v);
                    for (j = 0; j < n_quad_v; ++j)
                         B_l[m][ii+n_quad_v+j] -=  v2[j];

                    xm_v_mul(X_i_12_l[i][m], &B[jj], n_quad_v, n_quad_v, v2);
                    for (j = 0; j < n_quad_v; ++j)
                         B_l[m][ii+n_quad_v+j] -= -v2[j];

                    xm_v_mul(X_i_12[i], &B[jj+n_quad_v2], n_quad_v, n_quad_v, v2);
                    xm_v_diag_mul(v1, v2, v3, n_quad_v);
                    xm_v_mul(X_i_12_l[i][m], &B[jj+n_quad_v2], n_quad_v, n_quad_v, v2);
                    xm_v_diag_mul(alpha[i], v2, v4, n_quad_v);
                    xvec_add(v3, v4, v2, n_quad_v);
                    for (j = 0; j < n_quad_v; ++j)
                         B_l[m][ii+n_quad_v+j] -=  v2[j];
               }


               /*----------------------------------------------------------
                *
                *--------------------------------------------------------*/
               if (derivs_beam[i][m]) {
                    for (j = 0; j < n_quad_v; ++j) {
                         bb   = nu[i][j] * ltau[i];
                         bb_l = nu[i][j] * ltau_l[i][m];
                         if (derivs_layers[i][m])
                              bb_l += nu_l[i][m][j] * ltau[i];

                         v1[j] = ((-btau_l[i  ][m] - bb_l) * btran[i  ] * alpha[i][j] + btau_l[i+1][m] * btran[i+1] - beta1[i][j] * (btau_l[i+1][m] - btau_l[i  ][m] - bb_l)) / (btau[i+1] - btau[i  ] - bb);
                         v2[j] = ((-btau_l[i+1][m] - bb_l) * btran[i+1] * alpha[i][j] + btau_l[i  ][m] * btran[i  ] - beta2[i][j] * (btau_l[i  ][m] - btau_l[i+1][m] - bb_l)) / (btau[i  ] - btau[i+1] - bb);
                    }
               }


               /*----------------------------------------------------------
                *
                *--------------------------------------------------------*/
               if (derivs_layers[i][m]) {
                    bb   = aa * omega[i] * ltau[i];

                    bb_l = aa * (omega_l[i][m] * ltau[i] + omega[i] * ltau_l[i][m]);

                    for (j = 0; j < n_quad_v; ++j) {
                         cc   = bb   / qx_v[j];
                         cc_l = bb_l / qx_v[j];
                         sigma_p_l[i][m][j]  =  cc_l * P_q0_pm  [i]   [j];
                         sigma_m_l[i][m][j]  = -cc_l * P_q0_mm  [i]   [j];
                         sigma_p_l[i][m][j] +=  cc   * P_q0_pm_l[i][m][j];
                         sigma_m_l[i][m][j] += -cc   * P_q0_mm_l[i][m][j];

                         if (j % n_stokes > 1)
                              sigma_m_l[i][m][j] = -sigma_m_l[i][m][j];
                    }
               }


               /*----------------------------------------------------------
                *
                *--------------------------------------------------------*/
               if (derivs_beam[i][m]) {
                    xm_v_mul(X_i_11[i], sigma_p[i], n_quad_v, n_quad_v, v3);
                    xm_v_mul(X_i_12[i], sigma_m[i], n_quad_v, n_quad_v, v4);
                    xvec_add(v3, v4, v3, n_quad_v);

                    xm_v_diag_mul(v1, v3, v4, n_quad_v);
                    xvec_add(&B_l[m][ii], v4, &B_l[m][ii], n_quad_v);
               }

               if (derivs_layers[i][m]) {
                    xm_v_mul(X_i_11_l[i][m], sigma_p[i], n_quad_v, n_quad_v, v3);
                    xm_v_mul(X_i_11[i], sigma_p_l[i][m], n_quad_v, n_quad_v, v4);
                    xvec_add(v3, v4, v3, n_quad_v);
                    xm_v_mul(X_i_12_l[i][m], sigma_m[i], n_quad_v, n_quad_v, v4);
                    xvec_add(v3, v4, v3, n_quad_v);
                    xm_v_mul(X_i_12[i], sigma_m_l[i][m], n_quad_v, n_quad_v, v4);
                    xvec_add(v3, v4, v3, n_quad_v);

                    xm_v_diag_mul(beta1[i], v3, v4, n_quad_v);
                    xvec_add(&B_l[m][ii], v4, &B_l[m][ii], n_quad_v);
               }

               if (derivs_beam[i][m]) {
                    xm_v_mul(X_i_12[i], sigma_p[i], n_quad_v, n_quad_v, v3);
                    xm_v_mul(X_i_11[i], sigma_m[i], n_quad_v, n_quad_v, v4);
                    xvec_add(v3, v4, v3, n_quad_v);

                    xm_v_diag_mul(v2, v3, v4, n_quad_v);
                    xvec_sub(&B_l[m][ii + n_quad_v], v4, &B_l[m][ii + n_quad_v], n_quad_v);
               }

               if (derivs_layers[i][m]) {
                    xm_v_mul(X_i_12_l[i][m], sigma_p[i], n_quad_v, n_quad_v, v3);
                    xm_v_mul(X_i_12[i], sigma_p_l[i][m], n_quad_v, n_quad_v, v4);
                    xvec_add(v3, v4, v3, n_quad_v);
                    xm_v_mul(X_i_11_l[i][m], sigma_m[i], n_quad_v, n_quad_v, v4);
                    xvec_add(v3, v4, v3, n_quad_v);
                    xm_v_mul(X_i_11[i], sigma_m_l[i][m], n_quad_v, n_quad_v, v4);
                    xvec_add(v3, v4, v3, n_quad_v);

                    xm_v_diag_mul(beta2[i], v3, v4, n_quad_v);
                    xvec_sub(&B_l[m][ii + n_quad_v], v4, &B_l[m][ii + n_quad_v], n_quad_v);
               }

               jj += n_quad_v2;

               ii += n_quad_v2;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
#ifndef USE_BANDED_SOLVER
          xmat_getrs(A, &B_l[m], n_comp, 1, ipiv);
#else
          xgbtrs_("N", &n_comp, &n_diags, &n_diags, &nrhs, (TYPE *) *A,
                  &m_comp, ipiv, B_l[m], &n_comp, &info);
          if (info) {
               fprintf(stderr, "ERROR: xgbtrs() info = %d\n", info);
               exit(1);
          }
#endif
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void CALC_RADIANCE_LEVELS(int n_quad, int n_layers, int n_derivs,
                                 int n_ulevels, int *ulevels, TYPE *B, TYPE **B_l,
                                 double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                                 uchar **derivs_beam) {

     int i;
     int ii;
     int iii;
     int j;
     int k;

     int n_quad2;

     n_quad2 = n_quad * 2;

     if (n_derivs > 0) {
          init_array3_d(I_p_l, n_ulevels, n_derivs, n_quad, 0.);
          init_array3_d(I_m_l, n_ulevels, n_derivs, n_quad, 0.);
     }

     for (i = 0; i < n_ulevels; ++i) {
          ii = ulevels[i];

          iii = ii * n_quad2;
          for (j = 0; j < n_quad; ++j) {
               I_p[i][j] = XREAL(B[iii + j         ]);
               I_m[i][j] = XREAL(B[iii + j + n_quad]);
          }

          for (j = 0; j < n_derivs; ++j) {
               for (k = 0; k < n_quad; ++k) {
                    I_p_l[i][j][k] = XREAL(B_l[j][iii + k         ]);
                    I_m_l[i][j][k] = XREAL(B_l[j][iii + k + n_quad]);
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void CALC_RADIANCE_TAUS(int n_quad, int n_layers, int n_derivs,
                               int n_ulevels, int *ulevels, double *utaus,
                               double *ltau, double **ltau_l,
                               double *btau, double **btau_l,
                               TYPE **nu, TYPE ***X_p, TYPE ***X_m,
                               TYPE ***X_i_11, TYPE ***X_i_12, TYPE **sigma_p, TYPE **sigma_m,
                               TYPE ***nu_l, TYPE ****X_p_l, TYPE ****X_m_l,
                               TYPE ****X_i_11_l, TYPE ****X_i_12_l, TYPE ***sigma_p_l, TYPE ***sigma_m_l,
                               TYPE *B, TYPE **B_l,
                               double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                               uchar **derivs_layers, uchar **derivs_beam, work_data work) {

     int i;
     int ii;
     int iii;
     int j;
     int m;

     int n_quad2;

     TYPE aa;
     TYPE bb;
     TYPE cc;

     double x;
     double x_l;

     TYPE a1;
     TYPE a2;
     TYPE a3;
     TYPE a4;

     TYPE y1;
     TYPE y2;
     TYPE y3;

     TYPE *v1;
     TYPE *v2;
     TYPE *v3;
     TYPE *v4;
     TYPE *v5;
     TYPE *v6;

     TYPE *a;
     TYPE *b;
     TYPE *c;
     TYPE *d;

     TYPE *alpha1;
     TYPE *alpha2;
     TYPE *beta1;
     TYPE *beta2;

     TYPE *gpo;
     TYPE *hpp;


     n_quad2 = n_quad * 2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1     = get_work1(&work, WORK_XX);
     v2     = get_work1(&work, WORK_XX);
     v3     = get_work1(&work, WORK_XX);
     v4     = get_work1(&work, WORK_XX);
     v5     = get_work1(&work, WORK_XX);
     v6     = get_work1(&work, WORK_XX);

     a      = get_work1(&work, WORK_XX);
     b      = get_work1(&work, WORK_XX);
     c      = get_work1(&work, WORK_XX);
     d      = get_work1(&work, WORK_XX);

     alpha1 = get_work1(&work, WORK_XX);
     alpha2 = get_work1(&work, WORK_XX);
     beta1  = get_work1(&work, WORK_XX);
     beta2  = get_work1(&work, WORK_XX);

     gpo    = get_work1(&work, WORK_XX);
     hpp    = get_work1(&work, WORK_XX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_derivs > 0) {
          init_array3_d(I_p_l, n_ulevels, n_derivs, n_quad, 0.);
          init_array3_d(I_m_l, n_ulevels, n_derivs, n_quad, 0.);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_ulevels; ++i) {
          ii  = ulevels[i];

          iii = ii * n_quad2;

          x = utaus[i];


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (j = 0; j < n_quad; ++j) {
              alpha1[j] = XEXP(-nu[ii][j] *             x);
              alpha2[j] = XEXP(-nu[ii][j] * (ltau[ii] - x));
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          xm_v_mul(X_i_11[ii], &B[iii               ], n_quad, n_quad, v1);
          xm_v_mul(X_i_12[ii], &B[iii        +n_quad], n_quad, n_quad, v2);
          xvec_add(v1, v2, a, n_quad);

          xm_v_mul(X_i_12[ii], &B[iii+n_quad2       ], n_quad, n_quad, v1);
          xm_v_mul(X_i_11[ii], &B[iii+n_quad2+n_quad], n_quad, n_quad, v2);
          xvec_add(v1, v2, b, n_quad);


          for (j = 0; j < n_quad; ++j) {
              gpo[j]  = alpha1[j] *  a[j];
              hpp[j]  = alpha2[j] * -b[j];
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (j = 0; j < n_quad; ++j) {
              aa = nu[ii][j] *             x;
              beta1[j] = (XEXP(-(btau[ii  ] + aa)) - XEXP(-((1. - x / ltau[ii]) * btau[ii] + (x / ltau[ii]) * btau[ii+1]))) / (btau[ii+1] - btau[ii  ] - nu[ii][j] * ltau[ii]);
              aa = nu[ii][j] * (ltau[ii] - x);
              beta2[j] = (XEXP(-(btau[ii+1] + aa)) - XEXP(-((1. - x / ltau[ii]) * btau[ii] + (x / ltau[ii]) * btau[ii+1]))) / (btau[ii  ] - btau[ii+1] - nu[ii][j] * ltau[ii]);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          xm_v_mul(X_i_11[ii], sigma_p[ii], n_quad, n_quad, v1);
          xm_v_mul(X_i_12[ii], sigma_m[ii], n_quad, n_quad, v2);
          xvec_add(v1, v2, c, n_quad);

          xm_v_mul(X_i_12[ii], sigma_p[ii], n_quad, n_quad, v1);
          xm_v_mul(X_i_11[ii], sigma_m[ii], n_quad, n_quad, v2);
          xvec_add(v1, v2, d, n_quad);


          for (j = 0; j < n_quad; ++j) {
              gpo[j] += -beta1[j] *  c[j];
              hpp[j] +=  beta2[j] * -d[j];
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          xm_v_mul(X_p[ii], gpo, n_quad, n_quad, v3);
          xm_v_mul(X_m[ii], hpp, n_quad, n_quad, v4);
          for (j = 0; j < n_quad; ++j)
               I_p[i][j] =  XREAL(v3[j] + v4[j]);

          xm_v_mul(X_m[ii], gpo, n_quad, n_quad, v3);
          xm_v_mul(X_p[ii], hpp, n_quad, n_quad, v4);
          for (j = 0; j < n_quad; ++j)
               I_m[i][j] = -XREAL(v3[j] + v4[j]);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          for (m = 0; m < n_derivs; ++m) {

               x_l = utaus[i] * ltau_l[ii][m] / ltau[ii];


               xm_v_mul(X_i_11[ii],    &B_l[m][iii               ], n_quad, n_quad, v1);
               xm_v_mul(X_i_12[ii],    &B_l[m][iii        +n_quad], n_quad, n_quad, v2);
               xvec_add(v1, v2, v3, n_quad);

               if (derivs_layers[ii][m]) {
                    xm_v_mul(X_i_11_l[ii][m], &B   [iii               ], n_quad, n_quad, v1);
                    xvec_add(v3, v1, v3, n_quad);
                    xm_v_mul(X_i_12_l[ii][m], &B   [iii        +n_quad], n_quad, n_quad, v1);
                    xvec_add(v3, v1, v3, n_quad);
               }

               xm_v_mul(X_i_12[ii],    &B_l[m][iii+n_quad2       ], n_quad, n_quad, v1);
               xm_v_mul(X_i_11[ii],    &B_l[m][iii+n_quad2+n_quad], n_quad, n_quad, v2);
               xvec_add(v1, v2, v4, n_quad);

               if (derivs_layers[ii][m]) {
                    xm_v_mul(X_i_12_l[ii][m], &B   [iii+n_quad2       ], n_quad, n_quad, v1);
                    xvec_add(v4, v1, v4, n_quad);
                    xm_v_mul(X_i_11_l[ii][m], &B   [iii+n_quad2+n_quad], n_quad, n_quad, v1);
                    xvec_add(v4, v1, v4, n_quad);
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               for (j = 0; j < n_quad; ++j) {
                    v3[j]  =  alpha1[j] *  v3[j];
                    v4[j]  =  alpha2[j] * -v4[j];

                    if (derivs_layers[ii][m]) {
                         a1 = -(nu_l[ii][m][j] *             x  + nu[ii][j] *                  x_l)  * XEXP(-nu[ii][j] *             x);
                         a2 = -(nu_l[ii][m][j] * (ltau[ii] - x) + nu[ii][j] * (ltau_l[ii][m] - x_l)) * XEXP(-nu[ii][j] * (ltau[ii] - x));

                         v3[j] +=  a1        *  a [j];
                         v4[j] +=  a2        * -b [j];
                    }
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               if (derivs_layers[ii][m]) {
                    xm_v_mul(X_i_11_l[ii][m], sigma_p  [ii],    n_quad, n_quad, v1);
                    xm_v_mul(X_i_11  [ii],    sigma_p_l[ii][m], n_quad, n_quad, v2);
                    xvec_add(v1, v2, v5, n_quad);
                    xm_v_mul(X_i_12_l[ii][m], sigma_m  [ii],    n_quad, n_quad, v1);
                    xvec_add(v5, v1, v5, n_quad);
                    xm_v_mul(X_i_12  [ii],    sigma_m_l[ii][m], n_quad, n_quad, v1);
                    xvec_add(v5, v1, v5, n_quad);

                    xm_v_mul(X_i_12_l[ii][m], sigma_p  [ii],    n_quad, n_quad, v1);
                    xm_v_mul(X_i_12  [ii],    sigma_p_l[ii][m], n_quad, n_quad, v2);
                    xvec_add(v1, v2, v6, n_quad);
                    xm_v_mul(X_i_11_l[ii][m], sigma_m  [ii],    n_quad, n_quad, v1);
                    xvec_add(v6, v1, v6, n_quad);
                    xm_v_mul(X_i_11  [ii],    sigma_m_l[ii][m], n_quad, n_quad, v1);
                    xvec_add(v6, v1, v6, n_quad);

                    for (j = 0; j < n_quad; ++j) {
                         v3[j] += -beta1[j] *  v5[j];
                         v4[j] +=  beta2[j] * -v6[j];
                    }
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               if (derivs_beam[ii][m]) {
                    aa = ((-x_l / ltau[ii] + (x * ltau_l[ii][m] / (ltau[ii] * ltau[ii]))) * btau[ii] + (1. - x / ltau[ii]) * btau_l[ii][m] + (x_l / ltau[ii] - x * ltau_l[ii][m] / (ltau[ii] * ltau[ii])) * btau[ii+1] + x / ltau[ii] * btau_l[ii+1][m]) * XEXP(-((1. - x / ltau[ii]) * btau[ii] + (x / ltau[ii]) * btau[ii+1]));

                    for (j = 0; j < n_quad; ++j) {
                         y1 = nu[ii][j] * x;
                         y2 = 0.;
                         if (derivs_layers[ii][m])
                              y2 += nu_l[ii][m][j] * x + nu[ii][j] * x_l;
                         y3 = beta1[j];
                         bb = (-btau_l[ii  ][m] - y2) * XEXP(-(btau[ii  ] + y1)) + aa;
                         cc = btau_l[ii+1][m] - btau_l[ii  ][m];
                         if (derivs_layers[ii][m])
                              cc -= nu_l[ii][m][j] * ltau[ii] + nu[ii][j] * ltau_l[ii][m];
                         a3 = (bb - y3 * cc) / (btau[ii+1] - btau[ii  ] - nu[ii][j] * ltau[ii]);

                         y1 = nu[ii][j] * (ltau[ii] - x);
                         y2 = 0.;
                         if (derivs_layers[ii][m])
                              y2 += nu_l[ii][m][j] * (ltau[ii] - x) + nu[ii][j] * (ltau_l[ii][m] - x_l);
                         y3 = beta2[j];
                         bb = (-btau_l[ii+1][m] - y2) * XEXP(-(btau[ii+1] + y1)) + aa;
                         cc = btau_l[ii  ][m] - btau_l[ii+1][m];
                         if (derivs_layers[ii][m])
                              cc -= nu_l[ii][m][j] * ltau[ii] + nu[ii][j] * ltau_l[ii][m];
                         a4 = (bb - y3 * cc) / (btau[ii  ] - btau[ii+1] - nu[ii][j] * ltau[ii]);

                         v3[j] += -a3        *  c [j];
                         v4[j] +=  a4        * -d [j];
                    }
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               xm_v_mul(X_p[ii], v3, n_quad, n_quad, v1);
               xm_v_mul(X_m[ii], v4, n_quad, n_quad, v2);
               for (j = 0; j < n_quad; ++j)
                    I_p_l[i][m][j]  =  XREAL(v1[j] + v2[j]);

               if (derivs_layers[ii][m]) {
                    xm_v_mul(X_p_l[ii][m], gpo, n_quad, n_quad, v1);
                    for (j = 0; j < n_quad; ++j)
                         I_p_l[i][m][j] +=  XREAL(v1[j]);
                    xm_v_mul(X_m_l[ii][m], hpp, n_quad, n_quad, v1);
                    for (j = 0; j < n_quad; ++j)
                         I_p_l[i][m][j] +=  XREAL(v1[j]);
               }

               xm_v_mul(X_m[ii], v3, n_quad, n_quad, v1);
               xm_v_mul(X_p[ii], v4, n_quad, n_quad, v2);
               for (j = 0; j < n_quad; ++j)
                    I_m_l[i][m][j]  = -XREAL(v1[j] + v2[j]);

               if (derivs_layers[ii][m]) {
                    xm_v_mul(X_m_l[ii][m], gpo, n_quad, n_quad, v1);
                    for (j = 0; j < n_quad; ++j)
                         I_m_l[i][m][j] -=  XREAL(v1[j]);
                    xm_v_mul(X_p_l[ii][m], hpp, n_quad, n_quad, v1);
                    for (j = 0; j < n_quad; ++j)
                         I_m_l[i][m][j] -=  XREAL(v1[j]);
               }
          }
     }
}
