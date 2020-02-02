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
void layer_add_all_tl_with_ad(double **R12, double **T12, double *S12,
                              double **R21, double **T21, double *S21,
                              double **R23, double **T23, double *S23,
                              double **R32, double **T32, double *S32,
                              double **R13, double **T13, double *S13,
                              double **R31, double **T31, double *S31,
                              double ***R12_l, double ***T12_l, double **S12_l,
                              double ***R21_l, double ***T21_l, double **S21_l,
                              double ***R23_l, double ***T23_l, double **S23_l,
                              double ***R32_l, double ***T32_l, double **S32_l,
                              double ***R13_l, double ***T13_l, double **S13_l,
                              double ***R31_l, double ***T31_l, double **S31_l,
                              int n, int n_derivs, double atran, double *atran_l,
                              uchar *derivs, int flag, save_tree_data save_tree,
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

     double atran_a;

     double **R12_a;
     double **T12_a;
     double  *S12_a;
     double **R21_a;
     double **T21_a;
     double  *S21_a;

     double **R23_a;
     double **T23_a;
     double  *S23_a;
     double **R32_a;
     double **T32_a;
     double  *S32_a;

     double **R13_a;
     double **T13_a;
     double  *S13_a;
     double **R31_a;
     double **T31_a;
     double  *S31_a;

     double **a;
     double ***b;


     if (! flags_or(&derivs, n_derivs))
          return;


     save_tree_decode_s(&save_tree, "layer_add_all");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     R12_a = get_work1(&work, WORK_DXX);
     T12_a = get_work1(&work, WORK_DXX);
     S12_a = get_work1(&work, WORK_DX);
     R21_a = get_work1(&work, WORK_DXX);
     T21_a = get_work1(&work, WORK_DXX);
     S21_a = get_work1(&work, WORK_DX);

     R23_a = get_work1(&work, WORK_DXX);
     T23_a = get_work1(&work, WORK_DXX);
     S23_a = get_work1(&work, WORK_DX);
     R32_a = get_work1(&work, WORK_DXX);
     T32_a = get_work1(&work, WORK_DXX);
     S32_a = get_work1(&work, WORK_DX);

     R13_a = get_work1(&work, WORK_DXX);
     T13_a = get_work1(&work, WORK_DXX);
     S13_a = get_work1(&work, WORK_DX);

     R31_a = get_work1(&work, WORK_DXX);
     T31_a = get_work1(&work, WORK_DXX);
     S31_a = get_work1(&work, WORK_DX);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     dmat_zero(R13_a, n, n);
     dmat_zero(T13_a, n, n);
     dvec_zero(S13_a, n);

     dmat_zero(R31_a, n, n);
     dmat_zero(T31_a, n, n);
     dvec_zero(S31_a, n);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = 4 * n * n + 2 * n;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n * n) {
               i_input = 0;
               ii = i;
               ii_input = ii / n;
               jj_input = ii % n;
          }
          else
          if (i < 2 * n * n) {
               i_input = 1;
               ii = i - n * n;
               ii_input = ii / n;
               jj_input = ii % n;
          }
          else
          if (i < 3 * n * n) {
               i_input = 2;
               ii = i - 2 * n * n;
               ii_input = ii / n;
               jj_input = ii % n;
          }
          else
          if (i < 4 * n * n) {
               i_input = 3;
               ii = i - 3 * n * n;
               ii_input = ii / n;
               jj_input = ii % n;
          }
          else
          if (i < 4 * n * n + n) {
               i_input = 4;
               ii = i - 4 * n * n;
               ii_input = ii;
          }
          else {
               i_input = 5;
               ii = i - (4 * n * n + n);
               ii_input = ii;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          atran_a = 0.;

          dmat_zero(R12_a, n, n);
          dmat_zero(T12_a, n, n);
          dvec_zero(S12_a, n);

          dmat_zero(R21_a, n, n);
          dmat_zero(T21_a, n, n);
          dvec_zero(S21_a, n);

          dmat_zero(R23_a, n, n);
          dmat_zero(T23_a, n, n);
          dvec_zero(S23_a, n);

          dmat_zero(R32_a, n, n);
          dmat_zero(T32_a, n, n);
          dvec_zero(S32_a, n);


          if (i_input == 0)
               R13_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 1)
               T13_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 2)
               R31_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 3)
               T31_a[ii_input][jj_input] = 1.;
          else
          if (i_input == 4)
               S13_a[ii_input] = 1.;
          else
               S31_a[ii_input] = 1.;


          layer_add_all_a(R12, T12, S12, R21, T21, S21, R23, T23, S23, R32, T32, S32, R13, T13, S13, R31, T31, S31, R12_a, T12_a, S12_a, R21_a, T21_a, S21_a, R23_a, T23_a, S23_a, R32_a, T32_a, S32_a, R13_a, T13_a, S13_a, R31_a, T31_a, S31_a, n, atran, &atran_a, derivs[0], save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 4 || i_input == 5) {
               if (i_input == 4)
                    a = S13_l;
               else
                    a = S31_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    a[j][ii_input]  = 0.;

                    a[j][ii_input] += atran_l[j] * atran_a;

                    for (k = 0; k < n; ++k) {
                         if (derivs[j] == ADDING_S_S || derivs[j] == ADDING_B_S || derivs[j] == ADDING_B_B) {
                              a[j][ii_input] += S12_l[j][k] * S12_a[k];
                              a[j][ii_input] += S21_l[j][k] * S21_a[k];
                         }
                         if (derivs[j] == ADDING_U_S || derivs[j] == ADDING_S_S || derivs[j] == ADDING_U_B || derivs[j] == ADDING_B_S || derivs[j] == ADDING_B_B) {
                              a[j][ii_input] += S23_l[j][k] * S23_a[k];
                              a[j][ii_input] += S32_l[j][k] * S32_a[k];
                         }

                         for (l = 0; l < n; ++l) {
                              if (derivs[j] == ADDING_B_S || derivs[j] == ADDING_B_B) {
                                   a[j][ii_input] += R12_l[j][k][l] * R12_a[k][l];
                                   a[j][ii_input] += T12_l[j][k][l] * T12_a[k][l];
                                   a[j][ii_input] += R21_l[j][k][l] * R21_a[k][l];
                                   a[j][ii_input] += T21_l[j][k][l] * T21_a[k][l];
                              }
                              if (derivs[j] == ADDING_U_B || derivs[j] == ADDING_B_B) {
                                   a[j][ii_input] += R23_l[j][k][l] * R23_a[k][l];
                                   a[j][ii_input] += T23_l[j][k][l] * T23_a[k][l];
                                   a[j][ii_input] += R32_l[j][k][l] * R32_a[k][l];
                                   a[j][ii_input] += T32_l[j][k][l] * T32_a[k][l];
                              }
                         }
                    }
               }
          }

          if (i_input == 0 || i_input == 1 || i_input == 2 || i_input == 3) {
               if (i_input == 0)
                    b = R13_l;
               else
               if (i_input == 1)
                    b = T13_l;
               else
               if (i_input == 2)
                    b = R31_l;
               else
                    b = T31_l;

               for (j = 0; j < n_derivs; ++j) {
                    if (! derivs[j])
                         continue;

                    b[j][ii_input][jj_input]  = 0.;

                    b[j][ii_input][jj_input] += atran_l[j] * atran_a;

                    for (k = 0; k < n; ++k) {
                         if (derivs[j] == ADDING_S_S || derivs[j] == ADDING_B_S || derivs[j] == ADDING_B_B) {
                              b[j][ii_input][jj_input] += S12_l[j][k] * S12_a[k];
                              b[j][ii_input][jj_input] += S21_l[j][k] * S21_a[k];
                         }
                         if (derivs[j] == ADDING_U_S || derivs[j] == ADDING_S_S || derivs[j] == ADDING_U_B || derivs[j] == ADDING_B_S || derivs[j] == ADDING_B_B) {
                              b[j][ii_input][jj_input] += S23_l[j][k] * S23_a[k];
                              b[j][ii_input][jj_input] += S32_l[j][k] * S32_a[k];
                         }

                         for (l = 0; l < n; ++l) {
                              if (derivs[j] == ADDING_B_S || derivs[j] == ADDING_B_B) {
                                   b[j][ii_input][jj_input] += R12_l[j][k][l] * R12_a[k][l];
                                   b[j][ii_input][jj_input] += T12_l[j][k][l] * T12_a[k][l];
                                   b[j][ii_input][jj_input] += R21_l[j][k][l] * R21_a[k][l];
                                   b[j][ii_input][jj_input] += T21_l[j][k][l] * T21_a[k][l];
                              }

                              if (derivs[j] == ADDING_U_B || derivs[j] == ADDING_B_B) {
                                   b[j][ii_input][jj_input] += R23_l[j][k][l] * R23_a[k][l];
                                   b[j][ii_input][jj_input] += T23_l[j][k][l] * T23_a[k][l];
                                   b[j][ii_input][jj_input] += R32_l[j][k][l] * R32_a[k][l];
                                   b[j][ii_input][jj_input] += T32_l[j][k][l] * T32_a[k][l];
                              }
                         }
                    }
               }
          }
     }
}
