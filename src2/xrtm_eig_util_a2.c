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
void eig_1n_gen_real_tl_with_ad(int n_quad, int n_derivs,
                                double  **gamma, double  *evals, double  **evecs,
                                double ***gamma_l, double **evals_l, double ***evecs_l,
                                int eigen_solver, uchar *derivs, save_tree_data save_tree,
                                work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double **gamma_a;

     double *evals_a;
     double **evecs_a;

     double **a;
     double ***b;


     if (! flags_or(derivs, n_derivs))
          return;


     save_tree_decode_s(&save_tree, "eig_1n_gen_real");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     gamma_a = get_work1(&work, WORK_DXX);

     evals_a = get_work1(&work, WORK_DX);
     evecs_a = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(evals_a, n_quad);
     dmat_zero(evecs_a, n_quad, n_quad);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = n_quad + n_quad * n_quad;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }
          else {
               i_input = 1;
               ii = i - n_quad;
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          dmat_zero(gamma_a, n_quad, n_quad);


          if (i_input == 0)
               evals_a[ii_input] = 1.;
          else
               evecs_a[ii_input][jj_input] = 1.;


          eig_1n_gen_real_a(n_quad, gamma, evals, evecs, gamma_a, evals_a, evecs_a, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0) {
               if (i_input == 0)
                    a = evals_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    a[j][ii_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         for (l = 0; l < n_quad; ++l) {
                              a[j][ii_input] += gamma_l[j][k][l] * gamma_a[k][l];
                         }
                    }
               }
          }

          if (i_input == 1) {
               if (i_input == 1)
                    b = evecs_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    b[j][ii_input][jj_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         for (l = 0; l < n_quad; ++l) {
                              b[j][ii_input][jj_input] += gamma_l[j][k][l] * gamma_a[k][l];
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_gen_complex_tl_with_ad(int n_quad, int n_derivs, double **gamma, double  *evals_r, double *evals_i, double **evecs, double ***gamma_l, double **evals_r_l, double **evals_i_l, double ***evecs_l, int eigen_solver, uchar *derivs, save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double **gamma_a;

     double *evals_r_a;
     double *evals_i_a;
     double **evecs_a;

     double **a;
     double ***b;


     if (! flags_or(derivs, n_derivs))
          return;


     save_tree_decode_s(&save_tree, "eig_1n_gen_complex");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     gamma_a = get_work1(&work, WORK_DXX);

     evals_r_a = get_work1(&work, WORK_DX);
     evals_i_a = get_work1(&work, WORK_DX);
     evecs_a   = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(evals_r_a, n_quad);
     dvec_zero(evals_i_a, n_quad);
     dmat_zero(evecs_a, n_quad, n_quad);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * n_quad + n_quad * n_quad;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }
          else
          if (i < 2 * n_quad) {
               i_input = 1;
               ii = i - n_quad;
               ii_input = ii;
          }
          else {
               i_input = 2;
               ii = i - 2 * n_quad;
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          dmat_zero(gamma_a, n_quad, n_quad);


          if (i_input == 0)
               evals_r_a[ii_input] = 1.;
          else
          if (i_input == 1)
               evals_i_a[ii_input] = 1.;
          else
               evecs_a[ii_input][jj_input] = 1.;


          for (j = 0; j < n_quad; ++j) {
               if (evals_i[j] > 0.) {
                     evals_i_a[j    ] = -evals_i_a[j    ];
                     evals_i_a[j + 1] = -evals_i_a[j + 1];

                     for (k = 0; k < n_quad; ++k)
                          evecs_a[k][j + 1] = -evecs_a[k][j + 1];
               }
          }

          eig_1n_gen_complex_a(n_quad, gamma, evals_r, evals_i, evecs, gamma_a, evals_r_a, evals_i_a, evecs_a, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = evals_r_l;
               else
                    a = evals_i_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    a[j][ii_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         for (l = 0; l < n_quad; ++l) {
                              a[j][ii_input] += gamma_l[j][k][l] * gamma_a[k][l];
                         }
                    }
               }
          }

          if (i_input == 2) {
               if (i_input == 2)
                    b = evecs_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    b[j][ii_input][jj_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         for (l = 0; l < n_quad; ++l) {
                              b[j][ii_input][jj_input] += gamma_l[j][k][l] * gamma_a[k][l];
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_to_2n_real_tl_with_ad(int n_quad, int n_derivs,
                                  double  **tpr, double  **tmr,
                                  double  *evals, double  **evecs,
                                  double  *nu, double  **X_p, double  **X_m,
                                  double ***tpr_l, double ***tmr_l,
                                  double **evals_l, double ***evecs_l,
                                  double **nu_l, double ***X_p_l, double ***X_m_l,
                                  double **aux, double ***aul,
                                  uchar *derivs, save_tree_data save_tree, work_data work) {

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

     double  *evals_a;
     double **evecs_a;

     double *nu_a;
     double **X_p_a;
     double **X_m_a;

     double **a;
     double ***b;


     if (! flags_or(derivs, n_derivs))
          return;


     save_tree_decode_s(&save_tree, "eig_1n_to_2n_real");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     tpr_a   = get_work1(&work, WORK_DXX);
     tmr_a   = get_work1(&work, WORK_DXX);
     evals_a = get_work1(&work, WORK_DX);
     evecs_a = get_work1(&work, WORK_DXX);

     nu_a    = get_work1(&work, WORK_DX);
     X_p_a   = get_work1(&work, WORK_DXX);
     X_m_a   = get_work1(&work, WORK_DXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dvec_zero(nu_a, n_quad);
     dmat_zero(X_p_a, n_quad, n_quad);
     dmat_zero(X_m_a, n_quad, n_quad);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = n_quad + 2 * n_quad * n_quad;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }
          else
          if (i < n_quad + n_quad * n_quad) {
               i_input = 1;
               ii = i - n_quad;
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }
          else {
               i_input = 2;
               ii = i - (n_quad + n_quad * n_quad);
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          dvec_zero(evals_a, n_quad);
          dmat_zero(evecs_a, n_quad, n_quad);
          dmat_zero(tpr_a, n_quad, n_quad);
          dmat_zero(tmr_a, n_quad, n_quad);


          if (i_input == 0)
               nu_a[ii_input] = 1.;
          else
          if (i_input == 1)
               X_p_a[ii_input][jj_input] = 1.;
          else
               X_m_a[ii_input][jj_input] = 1.;


          eig_1n_to_2n_real_a(n_quad, tpr, tmr, evals, evecs, nu, X_p, X_m, tpr_a, tmr_a, evals_a, evecs_a, nu_a, X_p_a, X_m_a, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0) {
               if (i_input == 0)
                    a = nu_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    a[j][ii_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         a[j][ii_input] += evals_l[j][k] * evals_a[k];

                         for (l = 0; l < n_quad; ++l) {
                              a[j][ii_input] += tpr_l[j][k][l] * tpr_a[k][l];
                              a[j][ii_input] += tmr_l[j][k][l] * tmr_a[k][l];
                              a[j][ii_input] += evecs_l[j][k][l] * evecs_a[k][l];
                         }
                    }
               }
          }

          if (i_input == 1 || i_input == 2) {
               if (i_input == 1)
                    b = X_p_l;
               else
                    b = X_m_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    b[j][ii_input][jj_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         b[j][ii_input][jj_input] += evals_l[j][k] * evals_a[k];

                         for (l = 0; l < n_quad; ++l) {
                              b[j][ii_input][jj_input] += tpr_l[j][k][l] * tpr_a[k][l];
                              b[j][ii_input][jj_input] += tmr_l[j][k][l] * tmr_a[k][l];
                              b[j][ii_input][jj_input] += evecs_l[j][k][l] * evecs_a[k][l];
                         }
                    }
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void eig_1n_to_2n_complex_tl_with_ad(int n_quad, int n_derivs,
                                     double **tpr, double **tmr,
                                     double  *evals_r, double  *evals_i, double  **evecs,
                                     dcomplex  *nu, dcomplex  **X_p, dcomplex  **X_m,
                                     double ***tpr_l, double ***tmr_l,
                                     double **evals_r_l, double **evals_i_l, double ***evecs_l,
                                     dcomplex **nu_l, dcomplex ***X_p_l, dcomplex ***X_m_l,
                                     uchar *derivs, save_tree_data save_tree, work_data work) {

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

     double  *evals_r_a;
     double  *evals_i_a;
     double **evecs_a;

     dcomplex *nu_a;
     dcomplex **X_p_a;
     dcomplex **X_m_a;

     dcomplex **a;
     dcomplex ***b;


     if (! flags_or(derivs, n_derivs))
          return;


     save_tree_decode_s(&save_tree, "eig_1n_to_2n_complex");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     tpr_a     = get_work1(&work, WORK_DXX);
     tmr_a     = get_work1(&work, WORK_DXX);
     evals_r_a = get_work1(&work, WORK_DX);
     evals_i_a = get_work1(&work, WORK_DX);
     evecs_a   = get_work1(&work, WORK_DXX);

     nu_a      = get_work1(&work, WORK_ZX);
     X_p_a     = get_work1(&work, WORK_ZXX);
     X_m_a     = get_work1(&work, WORK_ZXX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     zvec_zero(nu_a, n_quad);
     zmat_zero(X_p_a, n_quad, n_quad);
     zmat_zero(X_m_a, n_quad, n_quad);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 2 * n_quad + 4 * n_quad * n_quad;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_quad) {
               i_input = 0;
               ii = i;
               ii_input = ii;
          }
          else
          if (i < 2 * n_quad) {
               i_input = 1;
               ii = i - n_quad;
               ii_input = ii;
          }
          else
          if (i < 2 * n_quad + n_quad * n_quad) {
               i_input = 2;
               ii = i - 2 * n_quad;
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }
          else
          if (i < 2 * n_quad + 2 * n_quad * n_quad) {
               i_input = 3;
               ii = i - (2 * n_quad + n_quad * n_quad);
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }
          else
          if (i < 2 * n_quad + 3 * n_quad * n_quad) {
               i_input = 4;
               ii = i - (2 * n_quad + 2 * n_quad * n_quad);
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }
          else {
               i_input = 5;
               ii = i - (2 * n_quad + 3 * n_quad * n_quad);
               ii_input = ii / n_quad;
               jj_input = ii % n_quad;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          dvec_zero(evals_r_a, n_quad);
          dvec_zero(evals_i_a, n_quad);
          dmat_zero(evecs_a, n_quad, n_quad);
          dmat_zero(tpr_a, n_quad, n_quad);
          dmat_zero(tmr_a, n_quad, n_quad);


          if (i_input == 0)
               nu_a[ii_input] = 1.;
          else
          if (i_input == 1)
               nu_a[ii_input] = _Complex_I * 1.;
          else
          if (i_input == 2)
               X_p_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 3)
               X_p_a[ii_input][jj_input] = _Complex_I * 1.;
          else
          if (i_input == 4)
               X_m_a[ii_input][jj_input] = 1.;
          else
               X_m_a[ii_input][jj_input] = _Complex_I * 1.;


          for (j = 0; j < n_quad; ++j) {
                nu_a[j] = conj(nu_a[j]);

                for (k = 0; k < n_quad; ++k) {
                     X_p_a[j][k] = conj(X_p_a[j][k]);
                     X_m_a[j][k] = conj(X_m_a[j][k]);
                }
          }

          eig_1n_to_2n_complex_a(n_quad, tpr, tmr, evals_r, evals_i, evecs, nu, X_p, X_m, tpr_a, tmr_a, evals_r_a, evals_i_a, evecs_a, nu_a, X_p_a, X_m_a, save_tree, work);

          for (j = 0; j < n_quad; ++j) {
               if (evals_i[j] > 0.) {
                     evals_i_a[j    ] = -evals_i_a[j    ];
                     evals_i_a[j + 1] = -evals_i_a[j + 1];

                     for (k = 0; k < n_quad; ++k)
                          evecs_a[k][j + 1] = -evecs_a[k][j + 1];
               }
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0 || i_input == 1) {
               if (i_input == 0)
                    a = nu_l;
               else
                    a = nu_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    if (i_input == 0)
                         a[j][ii_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         if (i_input == 0) {
                              a[j][ii_input] += evals_r_l[j][k] * evals_r_a[k];
                              a[j][ii_input] += evals_i_l[j][k] * evals_i_a[k];
                         }
                         else {
                              a[j][ii_input] += _Complex_I * evals_r_l[j][k] * evals_r_a[k];
                              a[j][ii_input] += _Complex_I * evals_i_l[j][k] * evals_i_a[k];
                         }

                         for (l = 0; l < n_quad; ++l) {
                              if (i_input == 0) {
                                   a[j][ii_input] += tpr_l[j][k][l] * tpr_a[k][l];
                                   a[j][ii_input] += tmr_l[j][k][l] * tmr_a[k][l];
                              }
                              else {
                                   a[j][ii_input] += _Complex_I * tpr_l[j][k][l] * tpr_a[k][l];
                                   a[j][ii_input] += _Complex_I * tmr_l[j][k][l] * tmr_a[k][l];
                              }

                              if (i_input == 0)
                                   a[j][ii_input] += evecs_l[j][k][l] * evecs_a[k][l];
                              else
                                   a[j][ii_input] += _Complex_I * evecs_l[j][k][l] * evecs_a[k][l];
                         }
                    }
               }
          }

          if (i_input == 2 || i_input == 3 || i_input == 4 || i_input == 5) {
               if (i_input == 2)
                    b = X_p_l;
               else
               if (i_input == 3)
                    b = X_p_l;
               else
               if (i_input == 4)
                    b = X_m_l;
               else
                    b = X_m_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    if (i_input == 2 || i_input == 4)
                         b[j][ii_input][jj_input] = 0.;

                    for (k = 0; k < n_quad; ++k) {
                         if (i_input == 2 || i_input == 4) {
                              b[j][ii_input][jj_input] += evals_r_l[j][k] * evals_r_a[k];
                              b[j][ii_input][jj_input] += evals_i_l[j][k] * evals_i_a[k];
                         }
                         else {
                              b[j][ii_input][jj_input] += _Complex_I * evals_r_l[j][k] * evals_r_a[k];
                              b[j][ii_input][jj_input] += _Complex_I * evals_i_l[j][k] * evals_i_a[k];
                         }

                         for (l = 0; l < n_quad; ++l) {
                              if (i_input == 2 || i_input == 4) {
                                   b[j][ii_input][jj_input] += tpr_l[j][k][l] * tpr_a[k][l];
                                   b[j][ii_input][jj_input] += tmr_l[j][k][l] * tmr_a[k][l];
                              }
                              else {
                                   b[j][ii_input][jj_input] += _Complex_I * tpr_l[j][k][l] * tpr_a[k][l];
                                   b[j][ii_input][jj_input] += _Complex_I * tmr_l[j][k][l] * tmr_a[k][l];
                              }

                              if (i_input == 2 || i_input == 4)
                                   b[j][ii_input][jj_input] += evecs_l[j][k][l] * evecs_a[k][l];
                              else
                                   b[j][ii_input][jj_input] += _Complex_I * evecs_l[j][k][l] * evecs_a[k][l];
                         }
                    }
               }
          }
     }
}
