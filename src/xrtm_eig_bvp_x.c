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
                      double *ltau, double **ltau_l,
                      double *atran, double **atran_l,
                      TYPE **nu, TYPE ***X_p, TYPE ***X_m,
                      double **Fs_p, double **Fs_m, double **Fs1_p, double **Fs1_m,
                      double **Ft0_p, double **Ft0_m, double **Ft1_p, double **Ft1_m,
                      TYPE ***nu_l, TYPE ****X_p_l, TYPE ****X_m_l,
                      double ***Fs_p_l, double ***Fs_m_l, double ***Fs1_p_l, double ***Fs1_m_l,
                      double ***Ft0_p_l, double ***Ft0_m_l, double ***Ft1_p_l, double ***Ft1_m_l,
                      double **Rs_qq, double ***Rs_qq_l,
                      TYPE *B, TYPE **B_l,
                      double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
                      int surface, int solar, int thermal,
                      uchar **derivs_layers, uchar **derivs_beam, uchar **derivs_thermal,
                      save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int iii;
     int j;
     int jj;
     int k;
     int l;
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

     TYPE a;
     TYPE b;
     TYPE c;
     TYPE d;
     TYPE e;

     TYPE **w1;

     TYPE **lambda;
     TYPE **lambda_l;

     TYPE **A;

     FORWARD_SAVE_SOLVE_BVP_DATA *save;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_quad_v  = n_quad * n_stokes;
     n_quad_v2 = n_quad_v * 2;
     n_quad_v3 = n_quad_v * 3;

#ifndef USE_BANDED_SOLVER
     m_comp    = 2 * n_quad_v * n_layers;
#else
     n_diags   = 3 * n_quad_v - 1;
     m_comp    = 3 * n_diags + 1;
#endif
     n_comp    = 2 * n_quad_v * n_layers;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, SOLVE_BVP_SAVE_TREE_STRING);

          if (save_tree_retrieve_data(&save_tree, FORWARD_SAVE_SOLVE_BVP_DATA, &save))
               FORWARD_SAVE_SOLVE_BVP_ALLOC(save, n_layers, n_quad_v, m_comp, n_comp);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     w1     = get_work1(&work, WORK_XXX);

     lambda = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);

     ipiv   = get_work_i1(&work, n_comp);

     A      = get_work_x2(&work, n_comp, m_comp);

     if (flags_or2(derivs_layers, n_layers, n_derivs))
          lambda_l = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_quad_v; ++j) {
               lambda[i][j] = XEXP(-nu[i][j] * ltau[i]);
          }
     }

     if (save_tree.t) {
          for (i = 0; i < n_layers; ++i)
               XCAT(copy_array1_, TYPE_POSTFIX)(save->lambda[i], lambda[i], n_quad_v);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     xmat_zero(A, n_comp, m_comp);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i = 0;
#ifndef USE_BANDED_SOLVER
     insert_matrix1_x(n_quad_v, n_quad_v, 0, 0, A, -1.,            X_m[i]);
     insert_matrix2_x(n_quad_v, n_quad_v, 0, 1, A, -1., lambda[i], X_p[i]);
#else
     insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, 0, 0, A, -1.,            X_m[i]);
     insert_matrix_band_storage2_x(n_quad_v, n_quad_v, n_diags, n_diags, 0, 1, A, -1., lambda[i], X_p[i]);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     ii = 1;
     jj = 0;

     for (i = 0; i < n_layers - 1; ++i) {
#ifndef USE_BANDED_SOLVER
          insert_matrix2_x(n_quad_v, n_quad_v, ii, jj,     A,  1., lambda[i    ], X_p[i    ]);
          insert_matrix1_x(n_quad_v, n_quad_v, ii, jj + 1, A,  1.,                X_m[i    ]);
          insert_matrix1_x(n_quad_v, n_quad_v, ii, jj + 2, A, -1.,                X_p[i + 1]);
          insert_matrix2_x(n_quad_v, n_quad_v, ii, jj + 3, A, -1., lambda[i + 1], X_m[i + 1]);
#else
          insert_matrix_band_storage2_x(n_quad_v, n_quad_v, n_diags, n_diags, ii, jj,     A,  1., lambda[i    ], X_p[i    ]);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii, jj + 1, A,  1.,                X_m[i    ]);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii, jj + 2, A, -1.,                X_p[i + 1]);
          insert_matrix_band_storage2_x(n_quad_v, n_quad_v, n_diags, n_diags, ii, jj + 3, A, -1., lambda[i + 1], X_m[i + 1]);
#endif
          ii += 1;

#ifndef USE_BANDED_SOLVER
          insert_matrix2_x(n_quad_v, n_quad_v, ii, jj,     A, -1., lambda[i    ], X_m[i    ]);
          insert_matrix1_x(n_quad_v, n_quad_v, ii, jj + 1, A, -1.,                X_p[i    ]);
          insert_matrix1_x(n_quad_v, n_quad_v, ii, jj + 2, A,  1.,                X_m[i + 1]);
          insert_matrix2_x(n_quad_v, n_quad_v, ii, jj + 3, A,  1., lambda[i + 1], X_p[i + 1]);
#else
          insert_matrix_band_storage2_x(n_quad_v, n_quad_v, n_diags, n_diags, ii, jj,     A, -1., lambda[i    ], X_m[i    ]);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii, jj + 1, A, -1.,                X_p[i    ]);
          insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii, jj + 2, A,  1.,                X_m[i + 1]);
          insert_matrix_band_storage2_x(n_quad_v, n_quad_v, n_diags, n_diags, ii, jj + 3, A,  1., lambda[i + 1], X_p[i + 1]);
#endif
          ii += 1;

          jj += 2;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i  = n_layers - 1;

     ii = i * 2;
if (LEGACY_MODE) {
     for (j = 0; j < n_quad_v; ++j) {
          if (! surface) {
               for (k = 0; k < n_quad_v; ++k) {
                    w1[j][k] = X_p[i][j][k] * lambda[i][k];
               }
          }
          else {
               for (k = 0; k < n_quad_v; ++k) {
                    a = 0.;
                    for (l = 0; l < n_quad_v; ++l)
                         a += Rs_qq[j][l] * -X_m[i][l][k];

                    w1[j][k] = (X_p[i][j][k] - a) * lambda[i][k];
               }
          }
     }
#ifndef USE_BANDED_SOLVER
     insert_matrix1_x(n_quad_v, n_quad_v, ii + 1, ii, A, 1., w1);
#else
     insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, ii, A, 1., w1);
#endif
}
else {
     TYPE **pp;

     if (! surface)
          pp = X_p[i];
     else {
          pp = w1;
#ifdef REAL
          xmat_mul(Rs_qq, X_m[i], n_quad_v, n_quad_v, n_quad_v, w1);
          xmat_add(X_p[i], w1, w1, n_quad_v, n_quad_v);
#endif
     }
#ifndef USE_BANDED_SOLVER
     insert_matrix2_x(n_quad_v, n_quad_v, ii + 1, ii, A, 1., lambda[i], pp);
#else
     insert_matrix_band_storage2_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, ii, A, 1., lambda[i], pp);
#endif
}
if (LEGACY_MODE) {
     for (j = 0; j < n_quad_v; ++j) {
          if (! surface) {
               for (k = 0; k < n_quad_v; ++k) {
                    w1[j][k] = X_m[i][j][k];
               }
          }
          else {
               for (k = 0; k < n_quad_v; ++k) {
                    a = 0.;
                    for (l = 0; l < n_quad_v; ++l)
                         a += Rs_qq[j][l] * -X_p[i][l][k];

                    w1[j][k] = (X_m[i][j][k] - a);
               }
          }
     }
#ifndef USE_BANDED_SOLVER
     insert_matrix1_x(n_quad_v, n_quad_v, ii + 1, ii + 1, A, 1., w1);
#else
     insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, ii + 1, A, 1., w1);
#endif
}
else {
     TYPE **pp;

     if (! surface)
          pp = X_m[i];
     else {
          pp = w1;
#ifdef REAL
          xmat_mul(Rs_qq, X_p[i], n_quad_v, n_quad_v, n_quad_v, w1);
          xmat_add(X_m[i], w1, w1, n_quad_v, n_quad_v);
#endif
     }
#ifndef USE_BANDED_SOLVER
     insert_matrix1_x(n_quad_v, n_quad_v, ii + 1, ii + 1, A, 1., pp);
#else
     insert_matrix_band_storage1_x(n_quad_v, n_quad_v, n_diags, n_diags, ii + 1, ii + 1, A, 1., pp);
#endif
}
     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef USE_BANDED_SOLVER
     xmat_getrf(A, n_comp, n_comp, ipiv);
#else
     xgbtrf_(&n_comp, &n_comp, &n_diags, &n_diags, *A, &m_comp, ipiv, &info);
     if (info) {
          fprintf(stderr, "ERROR: xgbtrf() info = %d\n", info);
          exit(1);
     }
#endif
     if (save_tree.t) {
          copy_array1_i(save->ipiv, ipiv, n_comp);

          XCAT(copy_array2_, TYPE_POSTFIX)(save->A, A,  n_comp, m_comp);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     xvec_zero(B, n_comp);

     for (i = 0; i < n_quad_v; ++i)
          B[i] += I1_m[i];

     for (i = 0; i < n_quad_v; ++i)
          B[n_comp - n_quad_v + i] += In_p[i];


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (solar) {
          for (i = 0; i < n_quad_v; ++i)
               B[i] -= Fs_m[0][i];


          ii = n_quad_v;
          for (i = 0; i < n_layers - 1; ++i) {
               for (j = 0; j < n_quad_v; ++j)
                    B[ii++] += Fs_p[i + 1][j] - Fs1_p[i][j] * atran[i];

               for (j = 0; j < n_quad_v; ++j)
                    B[ii++] += Fs_m[i + 1][j] - Fs1_m[i][j] * atran[i];
          }


          i = n_layers - 1;

          for (j = 0; j < n_quad_v; ++j) {
               b = -Fs1_p[i][j] * atran[i];

               if (surface) {
                    c = 0.;
                    for (k = 0; k < n_quad_v; ++k)
                         c += Rs_qq[j][k] * Fs1_m[i][k];

                    b += c * atran[i];
               }

               B[n_comp - n_quad_v + j] += b;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (thermal) {
          for (i = 0; i < n_quad_v; ++i)
               B[i] -= Ft0_m[0][i];


          ii = n_quad_v;
          for (i = 0; i < n_layers - 1; ++i) {
               for (j = 0; j < n_quad_v; ++j)
                    B[ii++] += Ft0_p[i + 1][j] - Ft1_p[i][j];

               for (j = 0; j < n_quad_v; ++j)
                    B[ii++] += Ft0_m[i + 1][j] - Ft1_m[i][j];
          }


          i = n_layers - 1;

          for (j = 0; j < n_quad_v; ++j) {
               b = -Ft1_p[i][j];

               if (surface) {
                    c = 0.;
                    for (k = 0; k < n_quad_v; ++k)
                         c += Rs_qq[j][k] * Ft1_m[i][k];

                    b += c;
               }

               B[n_comp - n_quad_v + j] += b;
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
#ifndef USE_BANDED_SOLVER
     xmat_getrs(A, &B, n_comp, 1, ipiv);
#else
     xgbtrs_("N", &n_comp, &n_diags, &n_diags, &nrhs, *A, &m_comp, ipiv, B,
             &n_comp, &info);
     if (info) {
          fprintf(stderr, "ERROR: xgbtrs() info = %d\n", info);
          exit(1);
     }
#endif

     if (save_tree.t)
          XCAT(copy_array1_, TYPE_POSTFIX)(save->B, B, n_comp);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (m = 0; m < n_derivs; ++m) {

          for (i = 0; i < n_layers; ++i) {
               if (derivs_layers[i][m]) {
                    for (j = 0; j < n_quad_v; ++j) {
                         lambda_l[i][j] = -(nu_l[i][m][j] * ltau[i] + nu[i][j] * ltau_l[i][m]) * lambda[i][j];
                    }
               }
          }

          for (i = 0; i < n_quad_v; ++i) {
               a = 0.;
               if (solar && derivs_beam[0][m])
                    a -= Fs_m_l[0][m][i];
               if (thermal && derivs_thermal[0][m])
                    a -= Ft0_m_l[0][m][i];

               if (derivs_layers[0][m]) {
                    for (j = 0; j < n_quad_v; ++j) {
                         a -= B[j] * -X_m_l[0][m][i][j] + B[j + n_quad_v] * (lambda_l[0][j] * -X_p[0][i][j] + lambda[0][j] * -X_p_l[0][m][i][j]);
                    }
               }

               B_l[m][i] = a + I1_m_l[m][i];
          }


          iii = n_quad_v;
          for (i = 0; i < n_layers - 1; ++i) {
               ii = i * n_quad_v2;
               for (j = 0; j < n_quad_v; ++j) {
                    a = 0.;

                    if (solar && derivs_beam[i    ][m])
                         a -= Fs1_p_l[i][m][j] * atran[i] + Fs1_p[i][j] * atran_l[i][m];
                    if (solar && derivs_beam[i + 1][m])
                         a += Fs_p_l[i + 1][m][j];
                    if (thermal && derivs_thermal[i    ][m])
                         a -= Ft1_p_l[i][m][j];
                    if (thermal && derivs_thermal[i + 1][m])
                         a += Ft0_p_l[i + 1][m][j];

                    for (k = 0; k < n_quad_v; ++k) {
                         if (derivs_layers[i  ][m])
                              a -= B[ii + k            ] * (lambda_l[i    ][k] *  X_p[i    ][j][k] + lambda[i    ][k] * X_p_l[i    ][m][j][k]) + B[ii + k + n_quad_v ] *  X_m_l[i    ][m][j][k];
                         if (derivs_layers[i + 1][m])
                              a += B[ii + k + n_quad_v3] * (lambda_l[i + 1][k] *  X_m[i + 1][j][k] + lambda[i + 1][k] * X_m_l[i + 1][m][j][k]) + B[ii + k + n_quad_v2] *  X_p_l[i + 1][m][j][k];
                    }

                    B_l[m][iii++] = a;
               }

               for (j = 0; j < n_quad_v; ++j) {
                    a = 0.;

                    if (solar && derivs_beam[i    ][m])
                         a -= Fs1_m_l[i][m][j] * atran[i] + Fs1_m[i][j] * atran_l[i][m];
                    if (solar && derivs_beam[i + 1][m])
                         a += Fs_m_l[i + 1][m][j];
                    if (thermal && derivs_thermal[i    ][m])
                         a -= Ft1_m_l[i][m][j];
                    if (thermal && derivs_thermal[i + 1][m])
                         a += Ft0_m_l[i + 1][m][j];

                    for (k = 0; k < n_quad_v; ++k) {
                         if (derivs_layers[i  ][m])
                              a -= B[ii + k            ] * (lambda_l[i    ][k] * -X_m[i    ][j][k] + lambda[i    ][k] * -X_m_l[i    ][m][j][k]) + B[ii + k + n_quad_v ] * -X_p_l[i    ][m][j][k];
                         if (derivs_layers[i + 1][m])
                              a += B[ii + k + n_quad_v3] * (lambda_l[i + 1][k] * -X_p[i + 1][j][k] + lambda[i + 1][k] * -X_p_l[i + 1][m][j][k]) + B[ii + k + n_quad_v2] * -X_m_l[i + 1][m][j][k];
                    }

                    B_l[m][iii++] = a;
               }
          }

          i  = n_layers - 1;

          ii = n_comp - n_quad_v2;

          for (j = 0; j < n_quad_v; ++j) {
               b = 0.;
               if (solar && derivs_beam[i][m])
                    b += -Fs1_p_l[i][m][j] * atran[i] - Fs1_p[i][j] * atran_l[i][m];
               if (thermal && derivs_thermal[i][m])
                    b -= Ft1_p_l[i][m][j];

               if (surface) {
                    c = 0.;
                    for (k = 0; k < n_quad_v; ++k) {
                         if (solar && derivs_beam[i][m])
                              c += Rs_qq[j][k] * (Fs1_m_l[i][m][k] * atran[i] + Fs1_m[i][k] * atran_l[i][m]);
                         if (solar && derivs_layers[i + 1][m])
                              c += Rs_qq_l[m][j][k] * Fs1_m[i][k] * atran[i];
                         if (thermal && derivs_thermal[i   ][m])
                              c += Rs_qq[j][k] * Ft1_m_l[i][m][k];
                         if (thermal && derivs_layers[i + 1][m])
                              c += Rs_qq_l[m][j][k] * Ft1_m[i][k];
                    }

                    b += c;

                    b += In_p_l[m][j];
               }

               for (k = 0; k < n_quad_v; ++k) {
                    if (derivs_layers[i][m])
                         b -= B[ii + k] * (lambda_l[i][k] * X_p[i][j][k] + lambda[i][k] * X_p_l[i][m][j][k]) + B[ii + k + n_quad_v] * X_m_l[i][m][j][k];

                    if (surface) {
                         c = 0.;
                         d = 0.;
                         e = 0.;
                         for (l = 0; l < n_quad_v; ++l) {
                              c += Rs_qq[j][l] * -X_m[i][l][k];
/*
                              if (save_tree.t)
                                   save->c[j][k] = c;
*/
                              if (derivs_layers[i    ][m]) {
                                   d += Rs_qq[j][l] * -X_m_l[i][m][l][k];
                                   e += Rs_qq[j][l] * -X_p_l[i][m][l][k];
                              }
                              if (derivs_layers[i + 1][m]) {
                                   d += Rs_qq_l[m][j][l] * -X_m[i][l][k];
                                   e += Rs_qq_l[m][j][l] * -X_p[i][l][k];
                              }
                         }

                         if (derivs_layers[i][m])
                              b -= B[ii + k] * lambda_l[i][k] * -c;

                         b -= B[ii + k] * lambda[i][k] * -d + B[ii + k + n_quad_v] * -e;
                    }
               }

               B_l[m][n_comp - n_quad_v + j] = b;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
#ifndef USE_BANDED_SOLVER
          xmat_getrs(A, &B_l[m], n_comp, 1, ipiv);
#else
          xgbtrs_("N", &n_comp, &n_diags, &n_diags, &nrhs, *A, &m_comp,
                  ipiv, B_l[m], &n_comp, &info);
          if (info) {
               fprintf(stderr, "ERROR: xgbtrs() info = %d\n", info);
               exit(1);
          }
#endif
     }
#ifdef USE_AD_FOR_TL_CALC_SOLVE_BVP
     SOLVE_BVP_TL_WITH_AD(n_quad, n_stokes, n_derivs, n_layers, ltau, ltau_l, Rs_qq, Rs_qq_l, atran, atran_l, nu, X_p, X_m, Fs_p, Fs_m, Ft0_p, Ft0_m, Ft1_p, Ft1_m, nu_l, X_p_l, X_m_l, Fs_p_l, Fs_m_l, Ft0_p_l, Ft0_m_l, Ft1_p_l, Ft1_m_l, B, B_l, I1_m, I1_m_l, In_p, In_p_l, surface, thermal, derivs_layers, derivs_beam, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
void CALC_RADIANCE_LEVELS(int n_quad, int n_layers, int n_derivs,
                          int n_ulevels, int *ulevels,
                          double *ltau, double **ltau_l,
                          double *atran, double **atran_l,
                          TYPE **nu, TYPE ***X_p, TYPE ***X_m,
                          double **Fs_p, double **Fs_m, double **Fs1_p, double **Fs1_m,
                          double **Ft0_p, double **Ft0_m, double **Ft1_p, double **Ft1_m,
                          TYPE ***nu_l, TYPE ****X_p_l, TYPE ****X_m_l,
                          double ***Fs_p_l, double ***Fs_m_l, double ***Fs1_p_l, double ***Fs1_m_l,
                          double ***Ft0_p_l, double ***Ft0_m_l, double ***Ft1_p_l, double ***Ft1_m_l,
                          TYPE *B, TYPE **B_l,
                          double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                          int solar, int thermal,
                          uchar **derivs_layers, uchar **derivs_beam, uchar **derivs_thermal,
                          save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int iii;
     int j;
     int k;
     int m;

     int n_quad2;

     TYPE a;

     TYPE *v1;
     TYPE *v2;
     TYPE *v3;

     TYPE **lambda;

     TYPE **lambda_l;

     FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA *save;


     n_quad2 = n_quad * 2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (save_tree.t) {
          save_tree_encode_s(&save_tree, CALC_RADIANCE_LEVELS_SAVE_TREE_STRING);

          if (save_tree_retrieve_data(&save_tree, FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA, &save))
               FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC(save, n_layers, n_quad);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1 = get_work1(&work, WORK_XX);
     v2 = get_work1(&work, WORK_XX);
     v3 = get_work1(&work, WORK_XX);

     lambda = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);

     if (flags_or2(derivs_layers, n_layers, n_derivs))
          lambda_l = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_quad; ++j) {
               lambda[i][j] = XEXP(-nu[i][j] * ltau[i]);
          }
     }

     if (save_tree.t) {
          for (i = 0; i < n_layers; ++i)
               XCAT(copy_array1_, TYPE_POSTFIX)(save->lambda[i], lambda[i], n_quad);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_ulevels; ++i) {
          ii = ulevels[i];

          if (ii != n_layers) {
               iii = ii * n_quad2;
if (LEGACY_MODE) {
               for (j = 0; j < n_quad; ++j) {
                    a = 0.;
                    for (k = 0; k < n_quad; ++k) {
                         a += B[iii + k         ] *  X_p[ii][j][k];
                         a += B[iii + k + n_quad] *  X_m[ii][j][k] * lambda[ii][k];
                    }

                    I_p[i][j] = XREAL(a);
                    if (solar)
                         I_p[i][j] += Fs_p[ii][j];
                    if (thermal)
                         I_p[i][j] += Ft0_p[ii][j];
               }
}
else {
               xm_v_diag_mul(lambda[ii], &B[iii + n_quad], v1, n_quad);

               xm_v_mul(X_p[ii], &B[iii], n_quad, n_quad, v2);
               xm_v_mul(X_m[ii], v1, n_quad, n_quad, v3);

               for (j = 0; j < n_quad; ++j)
                    I_p[i][j] = XREAL(v2[j]) + XREAL(v3[j]) + Fs_p[ii][j];
}
if (LEGACY_MODE) {
               for (j = 0; j < n_quad; ++j) {
                    a = 0.;
                    for (k = 0; k < n_quad; ++k) {
                         a += B[iii + k         ] * -X_m[ii][j][k];
                         a += B[iii + k + n_quad] * -X_p[ii][j][k] * lambda[ii][k];
                    }

                    I_m[i][j] = XREAL(a);
                    if (solar)
                         I_m[i][j] += Fs_m[ii][j];
                    if (thermal)
                         I_m[i][j] += Ft0_m[ii][j];
               }
}
else {
               xm_v_mul(X_m[ii], &B[iii], n_quad, n_quad, v2);
               xm_v_mul(X_p[ii], v1, n_quad, n_quad, v3);

               for (j = 0; j < n_quad; ++j)
                    I_m[i][j] = - XREAL(v2[j]) - XREAL(v3[j]) + Fs_m[ii][j];
}
          }
          else {
               ii = n_layers - 1;

               iii = ii * n_quad2;
if (LEGACY_MODE) {
               for (j = 0; j < n_quad; ++j) {
                    a = 0.;
                    for (k = 0; k < n_quad; ++k) {
                         a += B[iii + k         ] *  X_p[ii][j][k] * lambda[ii][k];
                         a += B[iii + k + n_quad] *  X_m[ii][j][k];
                    }

                    I_p[i][j] = XREAL(a);
                    if (solar)
                         I_p[i][j] += Fs1_p[ii][j] * atran[ii];
                    if (thermal)
                         I_p[i][j] += Ft1_p[ii][j];
               }
}
else {
               xm_v_diag_mul(lambda[ii], &B[iii], v1, n_quad);

               xm_v_mul(X_p[ii], v1, n_quad, n_quad, v2);
               xm_v_mul(X_m[ii], &B[iii + n_quad], n_quad, n_quad, v3);

               for (j = 0; j < n_quad; ++j)
                    I_p[i][j] = XREAL(v2[j]) + XREAL(v3[j]) + Fs1_p[ii][j] * atran[ii];
}
if (LEGACY_MODE) {
               for (j = 0; j < n_quad; ++j) {
                    a = 0.;
                    for (k = 0; k < n_quad; ++k) {
                         a += B[iii + k         ] * -X_m[ii][j][k] * lambda[ii][k];
                         a += B[iii + k + n_quad] * -X_p[ii][j][k];
                    }

                    I_m[i][j] = XREAL(a);
                    if (solar)
                         I_m[i][j] += Fs1_m[ii][j] * atran[ii];
                    if (thermal)
                         I_m[i][j] += Ft1_m[ii][j];
               }
}
else {
               xm_v_mul(X_m[ii], v1, n_quad, n_quad, v2);
               xm_v_mul(X_p[ii], &B[iii + n_quad], n_quad, n_quad, v3);

               for (j = 0; j < n_quad; ++j)
                    I_m[i][j] = - XREAL(v2[j]) - XREAL(v3[j]) + Fs1_m[ii][j] * atran[ii];
}
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_derivs > 0) {
          for (i = 0; i < n_ulevels; ++i) {
               for (j = 0; j < n_derivs; ++j) {
                    init_array1_d(I_p_l[i][j], n_quad, 0.);
                    init_array1_d(I_m_l[i][j], n_quad, 0.);
               }
          }
     }

     for (m = 0; m < n_derivs; ++m) {

          for (i = 0; i < n_layers; ++i) {
               if (derivs_layers[i][m]) {
                    for (j = 0; j < n_quad; ++j) {
                         lambda_l[i][j] = -(nu_l[i][m][j] * ltau[i] + nu[i][j] * ltau_l[i][m]) * lambda[i][j];
                    }
               }
          }

          for (i = 0; i < n_ulevels; ++i) {
               ii = ulevels[i];

               if (ii != n_layers) {
                    iii = ii * n_quad2;

                    for (j = 0; j < n_quad; ++j) {
                         a = 0.;
                         for (k = 0; k < n_quad; ++k) {
                              a += B_l[m][iii + k         ] *  X_p[ii][j][k];
                              a += B_l[m][iii + k + n_quad] *  X_m[ii][j][k] * lambda[ii][k];

                              if (derivs_layers[ii][m]) {
                                   a += B[iii + k         ] *  X_p_l[ii][m][j][k];
                                   a += B[iii + k + n_quad] * (X_m_l[ii][m][j][k] * lambda[ii][k] + X_m[ii][j][k] * lambda_l[ii][k]);
                              }
                         }

                         I_p_l[i][m][j] = XREAL(a);

                         if (solar && derivs_beam[ii][m])
                              I_p_l[i][m][j] += Fs_p_l[ii][m][j];
                         if (thermal && derivs_thermal[ii][m])
                              I_p_l[i][m][j] += Ft0_p_l[ii][m][j];
                    }

                    for (j = 0; j < n_quad; ++j) {
                         a = 0.;
                         for (k = 0; k < n_quad; ++k) {
                              a += B_l[m][iii + k         ] *  -X_m[ii][j][k];
                              a += B_l[m][iii + k + n_quad] *  -X_p[ii][j][k] * lambda[ii][k];

                              if (derivs_layers[ii][m]) {
                                   a += B[iii + k         ] *  -X_m_l[ii][m][j][k];
                                   a += B[iii + k + n_quad] * (-X_p_l[ii][m][j][k] * lambda[ii][k] + -X_p[ii][j][k] * lambda_l[ii][k]);
                              }
                         }

                         I_m_l[i][m][j] = XREAL(a);

                         if (solar && derivs_beam[ii][m])
                              I_m_l[i][m][j] += Fs_m_l[ii][m][j];
                         if (thermal && derivs_thermal[ii][m])
                              I_m_l[i][m][j] += Ft0_m_l[ii][m][j];
                    }
               }
               else {
                    ii = n_layers - 1;

                    iii = ii * n_quad2;

                    for (j = 0; j < n_quad; ++j) {
                         a = 0.;
                         for (k = 0; k < n_quad; ++k) {
                              a += B_l[m][iii + k         ] *  X_p[ii][j][k] * lambda[ii][k];
                              a += B_l[m][iii + k + n_quad] *  X_m[ii][j][k];

                              if (derivs_layers[ii][m]) {
                                   a += B[iii + k         ] * (X_p_l[ii][m][j][k] * lambda[ii][k] + X_p[ii][j][k] * lambda_l[ii][k]);
                                   a += B[iii + k + n_quad] *  X_m_l[ii][m][j][k];
                              }
                         }

                         I_p_l[i][m][j] = XREAL(a);

                         if (solar && derivs_beam[ii][m])
                              I_p_l[i][m][j] += Fs1_p_l[ii][m][j] * atran[ii] + Fs1_p[ii][j] * atran_l[ii][m];
if (thermal && derivs_thermal[ii][m])
     I_p_l[i][m][j] += Ft1_p_l[ii][m][j];
                    }

                    for (j = 0; j < n_quad; ++j) {
                         a = 0.;
                         for (k = 0; k < n_quad; ++k) {
                              a += B_l[m][iii + k         ] *  -X_m[ii][j][k] * lambda[ii][k];
                              a += B_l[m][iii + k + n_quad] *  -X_p[ii][j][k];

                              if (derivs_layers[ii][m]) {
                                   a += B[iii + k         ] * (-X_m_l[ii][m][j][k] * lambda[ii][k] + -X_m[ii][j][k] * lambda_l[ii][k]);
                                   a += B[iii + k + n_quad] *  -X_p_l[ii][m][j][k];
                              }
                         }

                         I_m_l[i][m][j] = XREAL(a);

                         if (solar && derivs_beam[ii][m])
                              I_m_l[i][m][j] += Fs1_m_l[ii][m][j] * atran[ii] + Fs1_m[ii][j] * atran_l[ii][m];
if (thermal && derivs_thermal[ii][m])
     I_m_l[i][m][j] += Ft1_m_l[ii][m][j];
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void CALC_RADIANCE_TAUS(int n_quad, int n_layers, int n_derivs,
                        int n_ulevels, int *ulevels, double *utaus,
                        double *ltau, double **ltau_l,
                        double *as_0, double **as_0_l, double *atran, double **atran_l,
                        TYPE **nu, TYPE ***X_p, TYPE ***X_m,
                        double **Fs_p, double **Fs_m,
                        TYPE ***nu_l, TYPE ****X_p_l, TYPE ****X_m_l,
                        double ***Fs_p_l, double ***Fs_m_l,
                        TYPE *B, TYPE **B_l,
                        double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                        int solar, int thermal,
                        uchar **derivs_layers, uchar **derivs_beam, uchar **derivs_thermal,
                        save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int iii;
     int j;
     int k;
     int m;

     int n_quad2;

     TYPE a;

     double utau_l;

     TYPE *v1;
     TYPE *v2;
     TYPE *v3;
     TYPE *v4;
/*
     TYPE **lambda;

     TYPE **lambda_l;
*/
     TYPE *lambda_p;
     TYPE *lambda_m;

     n_quad2 = n_quad * 2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     v1 = get_work1(&work, WORK_XX);
     v2 = get_work1(&work, WORK_XX);
     v3 = get_work1(&work, WORK_XX);
     v4 = get_work1(&work, WORK_XX);
/*
     lambda = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);

     if (flags_or2(derivs_layers, n_layers, n_derivs))
          lambda_l = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);
*/
     lambda_p = get_work1(&work, WORK_XX);
     lambda_m = get_work1(&work, WORK_XX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     for (i = 0; i < n_layers; ++i) {
          for (j = 0; j < n_quad; ++j) {
               lambda[i][j] = XEXP(-nu[i][j] * ltau[i]);
          }
     }
*/

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_ulevels; ++i) {
          ii = ulevels[i];

          for (j = 0; j < n_quad; ++j) {
               lambda_p[j] = XEXP(-nu[ii][j] * (utaus[i] - 0.));
               lambda_m[j] = XEXP(-nu[ii][j] * (ltau[ii] - utaus[i]));
          }

          iii = ii * n_quad2;
if (LEGACY_MODE) {
          for (j = 0; j < n_quad; ++j) {
               a = 0.;
               for (k = 0; k < n_quad; ++k) {
                    a += B[iii + k         ] *  X_p[ii][j][k] * lambda_p[k];
                    a += B[iii + k + n_quad] *  X_m[ii][j][k] * lambda_m[k];
               }

               I_p[i][j] = XREAL(a);
               if (solar)
                    I_p[i][j] += Fs_p[ii][j] * XEXP(-utaus[i] * as_0[ii]);
          }
}
else {
          xm_v_diag_mul(lambda_p, &B[iii], v1, n_quad);
          xm_v_diag_mul(lambda_m, &B[iii + n_quad], v2, n_quad);

          xm_v_mul(X_p[ii], v1, n_quad, n_quad, v3);
          xm_v_mul(X_m[ii], v2, n_quad, n_quad, v4);

          for (j = 0; j < n_quad; ++j)
               I_p[i][j] = XREAL(v3[j]) + XREAL(v4[j]) + Fs_p[ii][j] * XEXP(-utaus[i] * as_0[ii]);
}
if (LEGACY_MODE) {
          for (j = 0; j < n_quad; ++j) {
               a = 0.;
               for (k = 0; k < n_quad; ++k) {
                    a += B[iii + k         ] * -X_m[ii][j][k] * lambda_p[k];
                    a += B[iii + k + n_quad] * -X_p[ii][j][k] * lambda_m[k];
               }

               I_m[i][j] = XREAL(a);
               if (solar)
                    I_m[i][j] += Fs_m[ii][j] * XEXP(-utaus[i] * as_0[ii]);
          }
}
else {
          xm_v_mul(X_m[ii], v1, n_quad, n_quad, v3);
          xm_v_mul(X_p[ii], v2, n_quad, n_quad, v4);

          for (j = 0; j < n_quad; ++j)
               I_m[i][j] = - XREAL(v3[j]) - XREAL(v4[j]) + Fs_m[ii][j] * XEXP(-utaus[i] * as_0[ii]);
}
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (n_derivs > 0) {
          for (i = 0; i < n_ulevels; ++i) {
               for (j = 0; j < n_derivs; ++j) {
                    init_array1_d(I_p_l[i][j], n_quad, 0.);
                    init_array1_d(I_m_l[i][j], n_quad, 0.);
               }
          }
     }

     for (m = 0; m < n_derivs; ++m) {
/*
          for (i = 0; i < n_layers; ++i) {
               if (derivs_layers[i][m]) {
                    for (j = 0; j < n_quad; ++j) {
                         lambda_l[i][j] = -(nu_l[i][m][j] * ltau[i] + nu[i][j] * ltau_l[i][m]) * lambda[i][j];
                    }
               }
          }
*/
          for (i = 0; i < n_ulevels; ++i) {
               ii  = ulevels[i];

               iii = ii * n_quad2;

               utau_l = utaus[i] * ltau_l[ii][m] / ltau[ii];

               for (j = 0; j < n_quad; ++j) {
                    a = 0.;
                    for (k = 0; k < n_quad; ++k) {
                         a += B_l[m][iii + k         ] *  X_p[ii][j][k] * XEXP(-nu[ii][k] * (utaus[i] - 0.));
                         a += B_l[m][iii + k + n_quad] *  X_m[ii][j][k] * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i]));

                         if (derivs_layers[ii][m]) {
                              a += B[iii + k         ] * ( X_p_l[ii][m][j][k] * XEXP(-nu[ii][k] * (utaus[i] - 0.))       +  X_p[ii][j][k] * -(nu_l[ii][m][k] * (utaus[i] - 0.)       + nu[ii][k] *                  utau_l ) * XEXP(-nu[ii][k] * (utaus[i] - 0.)));
                              a += B[iii + k + n_quad] * ( X_m_l[ii][m][j][k] * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i])) +  X_m[ii][j][k] * -(nu_l[ii][m][k] * (ltau[ii] - utaus[i]) + nu[ii][k] * (ltau_l[ii][m] - utau_l)) * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i])));
                         }
                    }

                    I_p_l[i][m][j] = XREAL(a);

                    if (solar && derivs_beam[ii][m])
                         I_p_l[i][m][j] += (double) (Fs_p_l[ii][m][j] * XEXP(-utaus[i] * as_0[ii]) + Fs_p[ii][j] * -(utau_l * as_0[ii] + utaus[i] * as_0_l[ii][m]) * XEXP(-utaus[i] * as_0[ii]));
               }

               for (j = 0; j < n_quad; ++j) {
                    a = 0.;
                    for (k = 0; k < n_quad; ++k) {
                         a += B_l[m][iii + k         ] * -X_m[ii][j][k] * XEXP(-nu[ii][k] * (utaus[i] - 0.));
                         a += B_l[m][iii + k + n_quad] * -X_p[ii][j][k] * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i]));

                         if (derivs_layers[ii][m]) {
                              a += B[iii + k         ] * (-X_m_l[ii][m][j][k] * XEXP(-nu[ii][k] * (utaus[i] - 0.))       + -X_m[ii][j][k] * -(nu_l[ii][m][k] * (utaus[i] - 0.)       + nu[ii][k] *                  utau_l ) * XEXP(-nu[ii][k] * (utaus[i] - 0.)));
                              a += B[iii + k + n_quad] * (-X_p_l[ii][m][j][k] * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i])) + -X_p[ii][j][k] * -(nu_l[ii][m][k] * (ltau[ii] - utaus[i]) + nu[ii][k] * (ltau_l[ii][m] - utau_l)) * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i])));
                         }
                    }

                    I_m_l[i][m][j] = XREAL(a);

                    if (solar && derivs_beam[ii][m])
                         I_m_l[i][m][j] += (double) (Fs_m_l[ii][m][j] * XEXP(-utaus[i] * as_0[ii]) + Fs_m[ii][j] * -(utau_l * as_0[ii] + utaus[i] * as_0_l[ii][m]) * XEXP(-utaus[i] * as_0[ii]));
               }
          }
     }
}
