/******************************************************************************%
**
**    Copyright (C) 2007-2012 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 3.  See the file COPYING for more details.
**
*******************************************************************************/


#ifdef __cplusplus
extern "C" {
#endif

void xgbtrf_(int *, int *, int *, int *, TYPE *, int *, int *, int *);
void xgbtrs_(const char *, int *, int *, int *, int *, TYPE *, int *, int *, TYPE *, int *, int *);

#ifdef __cplusplus
}
#endif


/*******************************************************************************
 *
 ******************************************************************************/
static void SOLVE_BVP(int n_quad, int n_stokes, int n_derivs, int n_layers,
            double *ltau, double **ltau_l,
            double **Rs_qq, double ***Rs_qq_l, 
            double *atran, double **atran_l,
            TYPE **nu, TYPE ***X_p, TYPE ***X_m,
            double **F_p, double **F_m,
            double **F0_p, double **F0_m, double **F1_p, double **F1_m,
            TYPE ***nu_l, TYPE ****X_p_l, TYPE ****X_m_l,
            double ***F_p_l, double ***F_m_l,
            double ***F0_p_l, double ***F0_m_l, double ***F1_p_l, double ***F1_m_l,
            TYPE *B, TYPE **B_l,
            double *I1_m, double **I1_m_l, double *In_p, double **In_p_l,
            int surface, int thermal, uchar **derivs_h, uchar **derivs_p,
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

     int n_diags;

     int m_comp;
     int n_comp;

     int info;
     int nrhs = 1;

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

     n_diags   = 3 * n_quad_v - 1;
     m_comp    = 3 * n_diags + 1;
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

     if (flags_or2(derivs_h, n_layers, n_derivs))
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

     copy_to_band_storage_x (A, X_m[i], -1., n_quad_v, n_quad_v, n_diags, n_diags, 0, 0);
     copy_to_band_storage2_x(A, X_p[i], -1., lambda[i], n_quad_v, n_quad_v, n_diags, n_diags, 0, n_quad_v);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     ii = n_quad_v; jj = 0;

     for (i = 0; i < n_layers - 1; ++i) {
          copy_to_band_storage2_x(A, X_p[i  ],  1., lambda[i  ], n_quad_v, n_quad_v, n_diags, n_diags, ii, jj);
          copy_to_band_storage_x (A, X_m[i  ],  1., n_quad_v, n_quad_v, n_diags, n_diags, ii, jj+n_quad_v);
          copy_to_band_storage_x (A, X_p[i+1], -1., n_quad_v, n_quad_v, n_diags, n_diags, ii, jj+n_quad_v2);
          copy_to_band_storage2_x(A, X_m[i+1], -1., lambda[i+1], n_quad_v, n_quad_v, n_diags, n_diags, ii, jj+n_quad_v3);
          ii += n_quad_v;


          copy_to_band_storage2_x(A, X_m[i  ], -1., lambda[i  ], n_quad_v, n_quad_v, n_diags, n_diags, ii, jj);
          copy_to_band_storage_x (A, X_p[i  ], -1., n_quad_v, n_quad_v, n_diags, n_diags, ii, jj+n_quad_v);
          copy_to_band_storage_x (A, X_m[i+1],  1., n_quad_v, n_quad_v, n_diags, n_diags, ii, jj+n_quad_v2);
          copy_to_band_storage2_x(A, X_p[i+1],  1., lambda[i+1], n_quad_v, n_quad_v, n_diags, n_diags, ii, jj+n_quad_v3);
          ii += n_quad_v;

          jj += n_quad_v2;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     i  = n_layers - 1;

     ii = i * n_quad_v2;

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

     copy_to_band_storage_x(A, w1, 1., n_quad_v, n_quad_v, n_diags, n_diags, ii + n_quad_v, ii);


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

     copy_to_band_storage_x(A, w1, 1., n_quad_v, n_quad_v, n_diags, n_diags, ii + n_quad_v, ii + n_quad_v);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     xgbtrf_(&n_comp, &n_comp, &n_diags, &n_diags, *A, &m_comp, ipiv, &info);
     if (info) {
          eprintf("ERROR: xgbtrf() info = %d\n", info);
          exit(1);
     }

     if (save_tree.t) {
          copy_array1_i(save->ipiv, ipiv, n_comp);

          XCAT(copy_array2_, TYPE_POSTFIX)(save->A,    A,    n_comp, m_comp);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_quad_v; ++i)
          B[i] = I1_m[i] - F_m[0][i];


     ii = n_quad_v;
     for (i = 0; i < n_layers - 1; ++i) {
          for (j = 0; j < n_quad_v; ++j)
               B[ii++] = F_p[i+1][j] - F_p[i][j] * atran[i];

          for (j = 0; j < n_quad_v; ++j)
               B[ii++] = F_m[i+1][j] - F_m[i][j] * atran[i];
     }


     i = n_layers - 1;

     for (j = 0; j < n_quad_v; ++j) {
          b = -F_p[i][j] * atran[i];

          if (surface) {
               c = 0.;
               for (k = 0; k < n_quad_v; ++k)
                    c += Rs_qq[j][k] * F_m[i][k];

               b += In_p[j] + c * atran[i];
          }

          B[n_comp - n_quad_v + j] = b;
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
if (thermal) {
     for (i = 0; i < n_quad_v; ++i)
          B[i] += -F0_m[0][i];


     ii = n_quad_v;
     for (i = 0; i < n_layers - 1; ++i) {
          for (j = 0; j < n_quad_v; ++j)
               B[ii++] += F0_p[i+1][j] - F1_p[i][j];

          for (j = 0; j < n_quad_v; ++j)
               B[ii++] += F0_m[i+1][j] - F1_m[i][j];
     }


     i = n_layers - 1;

     for (j = 0; j < n_quad_v; ++j) {
          b = -F1_p[i][j];

          if (surface) {
               c = 0.;
               for (k = 0; k < n_quad_v; ++k)
                    c += Rs_qq[j][k] * F1_m[i][k];

               b += c;
          }

          B[n_comp - n_quad_v + j] += b;
     }
}


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     xgbtrs_("N", &n_comp, &n_diags, &n_diags, &nrhs, *A, &m_comp, ipiv, B,
             &n_comp, &info);
     if (info) {
          eprintf("ERROR: xgbtrs() info = %d\n", info);
          exit(1);
     }

     if (save_tree.t)
          XCAT(copy_array1_, TYPE_POSTFIX)(save->B, B, n_comp);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flags_or2(derivs_p, n_layers, n_derivs)) {

          for (m = 0; m < n_derivs; ++m) {

               for (i = 0; i < n_layers; ++i) {
                    if (derivs_h[i][m]) {
                         for (j = 0; j < n_quad_v; ++j) {
                              lambda_l[i][j] = -(nu_l[i][m][j] * ltau[i] + nu[i][j] * ltau_l[i][m]) * lambda[i][j];
                         }
                    }
               }

               for (i = 0; i < n_quad_v; ++i) {
                    a = 0.;
                    if (derivs_p[0][m])
                         a = -F_m_l[0][m][i];
                    if (derivs_h[0][m]) {
                         for (j = 0; j < n_quad_v; ++j) {
                              a -= B[j] * -X_m_l[0][m][i][j] + B[j + n_quad_v] * (lambda_l[0][j] * -X_p[0][i][j] + lambda[0][j] * -X_p_l[0][m][i][j]);
                         }
                    }

                    B_l[m][i] = a;
               }


               iii = n_quad_v;
               for (i = 0; i < n_layers - 1; ++i) {
                    ii = i * n_quad_v2;
                    for (j = 0; j < n_quad_v; ++j) {
                         a = 0.;

                         if (derivs_p[i  ][m])
                              a += -(F_p_l[i][m][j] * atran[i] + F_p[i][j] * atran_l[i][m]);
                         if (derivs_p[i+1][m])
                              a +=   F_p_l[i+1][m][j];

                         for (k = 0; k < n_quad_v; ++k) {
                              if (derivs_h[i  ][m])
                                   a -= B[ii + k            ] * (lambda_l[i  ][k] *  X_p[i  ][j][k] + lambda[i  ][k] * X_p_l[i  ][m][j][k]) + B[ii + k + n_quad_v ] *  X_m_l[i  ][m][j][k];
                              if (derivs_h[i+1][m])
                                   a += B[ii + k + n_quad_v3] * (lambda_l[i+1][k] *  X_m[i+1][j][k] + lambda[i+1][k] * X_m_l[i+1][m][j][k]) + B[ii + k + n_quad_v2] *  X_p_l[i+1][m][j][k];
                         }

                         B_l[m][iii++] = a;
                    }

                    for (j = 0; j < n_quad_v; ++j) {
                         a = 0.;

                         if (derivs_p[i  ][m])
                              a += -(F_m_l[i][m][j] * atran[i] + F_m[i][j] * atran_l[i][m]);
                         if (derivs_p[i+1][m])
                              a +=   F_m_l[i+1][m][j];

                         for (k = 0; k < n_quad_v; ++k) {
                              if (derivs_h[i  ][m])
                                   a -= B[ii + k            ] * (lambda_l[i  ][k] * -X_m[i  ][j][k] + lambda[i  ][k] * -X_m_l[i  ][m][j][k]) + B[ii + k + n_quad_v ] * -X_p_l[i  ][m][j][k];
                              if (derivs_h[i+1][m])
                                   a += B[ii + k + n_quad_v3] * (lambda_l[i+1][k] * -X_p[i+1][j][k] + lambda[i+1][k] * -X_p_l[i+1][m][j][k]) + B[ii + k + n_quad_v2] * -X_m_l[i+1][m][j][k];
                         }

                         B_l[m][iii++] = a;
                    }
               }

               i  = n_layers - 1;

               ii = n_comp - n_quad_v2;

               for (j = 0; j < n_quad_v; ++j) {
                    b = 0.;
                    if (derivs_p[i][m])
                         b = -F_p_l[i][m][j] * atran[i] - F_p[i][j] * atran_l[i][m];

                    if (surface) {
                         c = 0.;
                         for (k = 0; k < n_quad_v; ++k) {
                              if (derivs_p[i  ][m])
                                   c += Rs_qq[j][k] * (F_m_l[i][m][k] * atran[i] + F_m[i][k] * atran_l[i][m]);
                              if (derivs_h[i+1][m])
                                   c += Rs_qq_l[m][j][k] * F_m[i][k] * atran[i];
                         }

                         b += c;

                         if (derivs_p[i+1][m])
                             b += In_p_l[m][j];
                    }

                    for (k = 0; k < n_quad_v; ++k) {
                         if (derivs_h[i][m])
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
                                   if (derivs_h[i  ][m]) {
                                        d += Rs_qq[j][l] * -X_m_l[i][m][l][k];
                                        e += Rs_qq[j][l] * -X_p_l[i][m][l][k];
                                   }
                                   if (derivs_h[i+1][m]) {
                                        d += Rs_qq_l[m][j][l] * -X_m[i][l][k];
                                        e += Rs_qq_l[m][j][l] * -X_p[i][l][k];
                                   }
                              }

                              if (derivs_h[i][m])
                                   b -= B[ii + k] * lambda_l[i][k] * -c;

                              b -= B[ii + k] * lambda[i][k] * -d + B[ii + k + n_quad_v] * -e;
                         }
                    }

                    B_l[m][n_comp - n_quad_v + j] = b;
               }


               /*---------------------------------------------------------------
                *
                *-------------------------------------------------------------*/
               xgbtrs_("N", &n_comp, &n_diags, &n_diags, &nrhs, *A, &m_comp,
                       ipiv, B_l[m], &n_comp, &info);
               if (info) {
                    eprintf("ERROR: xgbtrs() info = %d\n", info);
                    exit(1);
               }
          }
     }
#ifdef USE_AD_FOR_TL_CALC_SOLVE_BVP
     SOLVE_BVP_TL_WITH_AD(n_quad, n_stokes, n_derivs, n_layers, ltau, ltau_l, Rs_qq, Rs_qq_l, atran, atran_l, nu, X_p, X_m, F_p, F_m, F0_p, F0_m, F1_p, F1_m, nu_l, X_p_l, X_m_l, F_p_l, F_m_l, F0_p_l, F0_m_l, F1_p_l, F1_m_l, B, B_l, I1_m, I1_m_l, In_p, In_p_l, surface, thermal, derivs_h, derivs_p, save_tree, work);
#endif
}



/*******************************************************************************
 *
 ******************************************************************************/
static void CALC_RADIANCE_LEVELS(int n_quad, int n_layers, int n_derivs, int n_ulevels,
                       int *ulevels, double *ltau, double **ltau_l,
                       double *atran, double **atran_l,
                       TYPE **nu, TYPE ***X_p, TYPE ***X_m,
                       double **F_p, double **F_m,
                       double **F0_p, double **F0_m, double **F1_p, double **F1_m,
                       TYPE ***nu_l, TYPE ****X_p_l, TYPE ****X_m_l,
                       double ***F_p_l, double ***F_m_l,
                       double ***F0_p_l, double ***F0_m_l, double ***F1_p_l, double ***F1_m_l,
                       TYPE *B, TYPE **B_l,
                       double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                       int thermal, uchar **derivs_h, uchar **derivs_p,
                       save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int iii;
     int j;
     int k;
     int m;

     int n_quad2;

     TYPE a;

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
     lambda = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);

     if (flags_or2(derivs_h, n_layers, n_derivs))
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

          for (j = 0; j < n_quad; ++j) {
               a = 0.;
               for (k = 0; k < n_quad; ++k) {
                    a += B[iii + k         ] *  X_p[ii][j][k];
                    a += B[iii + k + n_quad] *  X_m[ii][j][k] * lambda[ii][k];
               }
if (! thermal)
               I_p[i][j] = XREAL(a) + F_p[ii][j];
else
               I_p[i][j] = XREAL(a) + F_p[ii][j] + F0_p[ii][j];

          }

          for (j = 0; j < n_quad; ++j) {
               a = 0.;
               for (k = 0; k < n_quad; ++k) {
                    a += B[iii + k         ] * -X_m[ii][j][k];
                    a += B[iii + k + n_quad] * -X_p[ii][j][k] * lambda[ii][k];
               }
if (! thermal)
               I_m[i][j] = XREAL(a) + F_m[ii][j];
else
               I_m[i][j] = XREAL(a) + F_m[ii][j] + F0_m[ii][j];

          }
}
else {
          ii = n_layers - 1;

          iii = ii * n_quad2;

          for (j = 0; j < n_quad; ++j) {
               a = 0.;
               for (k = 0; k < n_quad; ++k) {
                    a += B[iii + k         ] *  X_p[ii][j][k] * lambda[ii][k];
                    a += B[iii + k + n_quad] *  X_m[ii][j][k];
               }
if (! thermal)
               I_p[i][j] = XREAL(a) + F_p[ii][j] * atran[ii];
else
               I_p[i][j] = XREAL(a) + F_p[ii][j] * atran[ii] + F1_p[ii][j];

          }

          for (j = 0; j < n_quad; ++j) {
               a = 0.;
               for (k = 0; k < n_quad; ++k) {
                    a += B[iii + k         ] * -X_m[ii][j][k] * lambda[ii][k];
                    a += B[iii + k + n_quad] * -X_p[ii][j][k];
               }
if (! thermal)
               I_m[i][j] = XREAL(a) + F_m[ii][j] * atran[ii];
else
               I_m[i][j] = XREAL(a) + F_m[ii][j] * atran[ii] + F1_m[ii][j];

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
     if (flags_or2(derivs_p, n_layers, n_derivs)) {

          for (m = 0; m < n_derivs; ++m) {

               for (i = 0; i < n_layers; ++i) {
                    if (derivs_h[i][m]) {
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

                              if (derivs_h[ii][m]) {
                                   a += B[iii + k         ] *  X_p_l[ii][m][j][k];
                                   a += B[iii + k + n_quad] * (X_m_l[ii][m][j][k] * lambda[ii][k] + X_m[ii][j][k] * lambda_l[ii][k]);
                              }
                         }

                         I_p_l[i][m][j] = XREAL(a);

                         if (derivs_p[ii][m])
                              I_p_l[i][m][j] += F_p_l[ii][m][j];
                    }

                    for (j = 0; j < n_quad; ++j) {
                         a = 0.;
                         for (k = 0; k < n_quad; ++k) {
                              a += B_l[m][iii + k         ] *  -X_m[ii][j][k];
                              a += B_l[m][iii + k + n_quad] *  -X_p[ii][j][k] * lambda[ii][k];

                              if (derivs_h[ii][m]) {
                                   a += B[iii + k         ] *  -X_m_l[ii][m][j][k];
                                   a += B[iii + k + n_quad] * (-X_p_l[ii][m][j][k] * lambda[ii][k] + -X_p[ii][j][k] * lambda_l[ii][k]);
                              }
                         }

                         I_m_l[i][m][j] = XREAL(a);

                         if (derivs_p[ii][m])
                              I_m_l[i][m][j] += F_m_l[ii][m][j];
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

                              if (derivs_h[ii][m]) {
                                   a += B[iii + k         ] * (X_p_l[ii][m][j][k] * lambda[ii][k] + X_p[ii][j][k] * lambda_l[ii][k]);
                                   a += B[iii + k + n_quad] *  X_m_l[ii][m][j][k];
                              }
                         }

                         I_p_l[i][m][j] = XREAL(a);

                         if (derivs_p[ii][m])
                              I_p_l[i][m][j] += F_p_l[ii][m][j] * atran[ii] + F_p[ii][j] * atran_l[ii][m];
                    }

                    for (j = 0; j < n_quad; ++j) {
                         a = 0.;
                         for (k = 0; k < n_quad; ++k) {
                              a += B_l[m][iii + k         ] *  -X_m[ii][j][k] * lambda[ii][k];
                              a += B_l[m][iii + k + n_quad] *  -X_p[ii][j][k];

                              if (derivs_h[ii][m]) {
                                   a += B[iii + k         ] * (-X_m_l[ii][m][j][k] * lambda[ii][k] + -X_m[ii][j][k] * lambda_l[ii][k]);
                                   a += B[iii + k + n_quad] *  -X_p_l[ii][m][j][k];
                              }
                         }

                         I_m_l[i][m][j] = XREAL(a);

                         if (derivs_p[ii][m])
                              I_m_l[i][m][j] += F_m_l[ii][m][j] * atran[ii] + F_m[ii][j] * atran_l[ii][m];
                    }
}
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void CALC_RADIANCE_TAUS(int n_quad, int n_layers, int n_derivs, int n_ulevels,
                     int *ulevels, double *utaus, double *ltau, double **ltau_l,
                     double *as_0, double **as_0_l, double *atran, double **atran_l,
                     TYPE **nu, TYPE ***X_p, TYPE ***X_m,
                     double **F_p, double **F_m,
                     TYPE ***nu_l, TYPE ****X_p_l, TYPE ****X_m_l,
                     double ***F_p_l, double ***F_m_l,
                     TYPE *B, TYPE **B_l,
                     double **I_p, double **I_m, double ***I_p_l, double ***I_m_l,
                     uchar **derivs_h, uchar **derivs_p,
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
/*
     TYPE **lambda;

     TYPE **lambda_l;
*/

     n_quad2 = n_quad * 2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
/*
     lambda = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);

     if (flags_or2(derivs_h, n_layers, n_derivs))
          lambda_l = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);
*/

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

          iii = ii * n_quad2;

          for (j = 0; j < n_quad; ++j) {
               a = 0.;
               for (k = 0; k < n_quad; ++k) {
                    a += B[iii + k         ] *  X_p[ii][j][k] * XEXP(-nu[ii][k] * (utaus[i] - 0.));
                    a += B[iii + k + n_quad] *  X_m[ii][j][k] * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i]));
               }

               I_p[i][j] = XREAL(a + F_p[ii][j] * XEXP(-utaus[i] * as_0[ii]));
          }

          for (j = 0; j < n_quad; ++j) {
               a = 0.;
               for (k = 0; k < n_quad; ++k) {
                    a += B[iii + k         ] * -X_m[ii][j][k] * XEXP(-nu[ii][k] * (utaus[i] - 0.));
                    a += B[iii + k + n_quad] * -X_p[ii][j][k] * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i]));
               }

               I_m[i][j] = XREAL(a + F_m[ii][j] * XEXP(-utaus[i] * as_0[ii]));
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
     if (flags_or2(derivs_p, n_layers, n_derivs)) {

          for (m = 0; m < n_derivs; ++m) {
/*
               for (i = 0; i < n_layers; ++i) {
                    if (derivs_h[i][m]) {
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

                              if (derivs_h[ii][m]) {
                                   a += B[iii + k         ] * ( X_p_l[ii][m][j][k] * XEXP(-nu[ii][k] * (utaus[i] - 0.))       +  X_p[ii][j][k] * -(nu_l[ii][m][k] * (utaus[i] - 0.)       + nu[ii][k] *                  utau_l ) * XEXP(-nu[ii][k] * (utaus[i] - 0.)));
                                   a += B[iii + k + n_quad] * ( X_m_l[ii][m][j][k] * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i])) +  X_m[ii][j][k] * -(nu_l[ii][m][k] * (ltau[ii] - utaus[i]) + nu[ii][k] * (ltau_l[ii][m] - utau_l)) * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i])));
                              }
                         }

                         I_p_l[i][m][j] = XREAL(a);

                         if (derivs_p[ii][m])
                              I_p_l[i][m][j] += (double) (F_p_l[ii][m][j] * XEXP(-utaus[i] * as_0[ii]) + F_p[ii][j] * -(utau_l * as_0[ii] + utaus[i] * as_0_l[ii][m]) * XEXP(-utaus[i] * as_0[ii]));
                    }

                    for (j = 0; j < n_quad; ++j) {
                         a = 0.;
                         for (k = 0; k < n_quad; ++k) {
                              a += B_l[m][iii + k         ] * -X_m[ii][j][k] * XEXP(-nu[ii][k] * (utaus[i] - 0.));
                              a += B_l[m][iii + k + n_quad] * -X_p[ii][j][k] * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i]));

                              if (derivs_h[ii][m]) {
                                   a += B[iii + k         ] * (-X_m_l[ii][m][j][k] * XEXP(-nu[ii][k] * (utaus[i] - 0.))       + -X_m[ii][j][k] * -(nu_l[ii][m][k] * (utaus[i] - 0.)       + nu[ii][k] *                  utau_l ) * XEXP(-nu[ii][k] * (utaus[i] - 0.)));
                                   a += B[iii + k + n_quad] * (-X_p_l[ii][m][j][k] * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i])) + -X_p[ii][j][k] * -(nu_l[ii][m][k] * (ltau[ii] - utaus[i]) + nu[ii][k] * (ltau_l[ii][m] - utau_l)) * XEXP(-nu[ii][k] * (ltau[ii] - utaus[i])));
                              }
                         }

                         I_m_l[i][m][j] = XREAL(a);

                         if (derivs_p[ii][m])
                              I_m_l[i][m][j] += (double) (F_m_l[ii][m][j] * XEXP(-utaus[i] * as_0[ii]) + F_m[ii][j] * -(utau_l * as_0[ii] + utaus[i] * as_0_l[ii][m]) * XEXP(-utaus[i] * as_0[ii]));
                    }
               }
          }
     }
}
