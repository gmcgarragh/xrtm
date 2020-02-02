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
void ssr_up_layer_tl_with_ad(int i_layer, int n_stokes, int n_derivs, double utau, double *umus, int n_umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double **P, double ***P_l, double *I, double **I_l, double *I_ss, double **I_ss_l, uchar **derivs_layers, work_data work) {

     int i;
     int j;
     int k;

     int i_input;
     int n_inputs;

     int ii_input;

     double a;

     double *omega_a;
     double *ltau_a;
     double *btran_a;
     double *as_0_a;
     double **P_a;
     double *I_a;

     double *I_ss_a;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     omega_a = get_work_d1(&work, 128);
     ltau_a  = get_work_d1(&work, 128);
     btran_a = get_work_d1(&work, 128 + 1);
     as_0_a  = get_work_d1(&work, 128);
     P_a     = get_work_d2(&work, 128, n_umus * n_stokes);
     I_a     = get_work_d1(&work,      n_umus * n_stokes);

     I_ss_a  = get_work_d1(&work,      n_umus * n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(I_ss_a, n_umus * n_stokes, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = n_umus * n_stokes;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_umus * n_stokes) {
               i_input = 0;
               ii_input = i;
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array1_d(omega_a, 128,     0.);
          init_array1_d(ltau_a,  128,     0.);
          init_array1_d(btran_a, 128 + 1, 0.);
          init_array1_d(as_0_a,  128,     0.);
          init_array2_d(P_a,     128, n_umus * n_stokes, 0.);
          init_array1_d(I_a,          n_umus * n_stokes, 0.);


          if (i_input == 0)
               I_ss_a[ii_input] = 1.;

          ssr_up_layer_a(i_layer, n_stokes, utau, umus, n_umus, omega, omega_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, P, P_a, I, I_a, I_ss, I_ss_a, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0) {
               for (j = 0; j < n_derivs; ++j) {
                    a  = 0.;

                    a += omega_l[i_layer][j] * omega_a[i_layer];
                    a += ltau_l [i_layer][j] * ltau_a [i_layer];
                    a += btran_l[i_layer][j] * btran_a[i_layer];
                    a += as_0_l [i_layer][j] * as_0_a [i_layer];

                    for (k = 0; k < n_umus * n_stokes; ++k) {
                        if (derivs_layers[i_layer][j])
                             a += P_l[i_layer][j][k] * P_a[i_layer][k];
                        a += I_l[j][k] * I_a[k];
                    }

                    if (i_input == 0)
                         I_ss_l[j][ii_input] = a;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void single_scattered_radiance_up_tl_with_ad(int n_stokes, int n_derivs, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P, double ***P_l, double *I_in, double **I_in_l, double **I_ss, double ***I_ss_l, int utau_output, uchar **derivs_layers, uchar **derivs_beam_down, save_tree_data save_tree, work_data work) {

     int i;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double ***a;

     double *omega_a;
     double *ltau_a;
     double *btran_a;
     double *as_0_a;
     double **P_a;
     double *I_in_a;

     double **I_ss_a;


     save_tree_decode_s(&save_tree, "single_scattered_radiance_up");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     omega_a = get_work_d1(&work, n_layers);
     ltau_a  = get_work_d1(&work, n_layers);
     btran_a = get_work_d1(&work, n_layers + 1);
     as_0_a  = get_work_d1(&work, n_layers);
     P_a     = get_work_d2(&work, n_layers, n_umus * n_stokes);
     I_in_a  = get_work_d1(&work,           n_umus * n_stokes);

     I_ss_a  = get_work_d2(&work, n_ulevels, n_umus * n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_ss_a, n_ulevels, n_umus * n_stokes, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = n_ulevels * n_umus * n_stokes;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_ulevels * n_umus * n_stokes) {
               i_input = 0;
               ii_input = i / (n_umus * n_stokes);
               jj_input = i % (n_umus * n_stokes);
          }

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array1_d(omega_a, n_layers,     0.);
          init_array1_d(ltau_a,  n_layers,     0.);
          init_array1_d(btran_a, n_layers + 1, 0.);
          init_array1_d(as_0_a,  n_layers,     0.);
          init_array2_d(P_a,     n_layers, n_umus * n_stokes, 0.);
          init_array1_d(I_in_a,            n_umus * n_stokes, 0.);


          if (i_input == 0)
               I_ss_a[ii_input][jj_input] = 1.;

          single_scattered_radiance_up_a(n_stokes, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega, omega_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, atran, NULL, P, P_a, I_in, I_in_a, I_ss, I_ss_a, utau_output, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0) {
               if (i_input == 0)
                    a = I_ss_l;

               for (j = 0; j < n_derivs; ++j) {
                    a[j][ii_input][jj_input] = 0.;

                    for (k = 0; k < n_layers; ++k) {
                         a[j][ii_input][jj_input] += omega_l [k][j] * omega_a[k];
                         a[j][ii_input][jj_input] += ltau_l  [k][j] * ltau_a [k];
                         a[j][ii_input][jj_input] += btran_l [k][j] * btran_a[k];
                         a[j][ii_input][jj_input] += as_0_l  [k][j] * as_0_a [k];

                         for (l = 0; l < n_umus * n_stokes; ++l) {
                             if (derivs_layers[k][j])
                                  a[j][ii_input][jj_input] += P_l[k][j][l] * P_a[k][l];
                         }
                    }

                    a[j][ii_input][jj_input] += btran_l[k][j] * btran_a[k];

                    for (k = 0; k < n_umus * n_stokes; ++k)
                        a[j][ii_input][jj_input] += I_in_l[j][k] * I_in_a[k];
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void ssr_dn_layer_tl_with_ad(int i_layer, int n_stokes, int n_derivs, double utau, double *umus, int n_umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double **P, double ***P_l, double *I, double **I_l, double *I_ss, double **I_ss_l, uchar **derivs_layers, work_data work) {

     int i;
     int j;
     int k;

     int i_input;
     int n_inputs;

     int ii_input;

     double a;

     double *omega_a;
     double *ltau_a;
     double *btran_a;
     double *as_0_a;
     double **P_a;
     double *I_a;

     double *I_ss_a;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     omega_a = get_work_d1(&work, 128);
     ltau_a  = get_work_d1(&work, 128);
     btran_a = get_work_d1(&work, 128 + 1);
     as_0_a  = get_work_d1(&work, 128);
     P_a     = get_work_d2(&work, 128, n_umus * n_stokes);
     I_a     = get_work_d1(&work,      n_umus * n_stokes);

     I_ss_a  = get_work_d1(&work,      n_umus * n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array1_d(I_ss_a, n_umus * n_stokes, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = n_umus * n_stokes;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_umus * n_stokes) {
               i_input = 0;
               ii_input = i;
          }

          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array1_d(omega_a, 128,     0.);
          init_array1_d(ltau_a,  128,     0.);
          init_array1_d(btran_a, 128 + 1, 0.);
          init_array1_d(as_0_a,  128,     0.);
          init_array2_d(P_a,     128, n_umus * n_stokes, 0.);
          init_array1_d(I_a,          n_umus * n_stokes, 0.);


          if (i_input == 0)
               I_ss_a[ii_input] = 1.;


          ssr_dn_layer_a(i_layer, n_stokes, utau, umus, n_umus, omega, omega_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, P, P_a, I, I_a, I_ss, I_ss_a, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0) {
               for (j = 0; j < n_derivs; ++j) {
                    a  = 0.;

                    a += omega_l[i_layer][j] * omega_a[i_layer];
                    a += ltau_l [i_layer][j] * ltau_a [i_layer];
                    a += btran_l[i_layer][j] * btran_a[i_layer];
                    a += as_0_l [i_layer][j] * as_0_a [i_layer];

                    for (k = 0; k < n_umus * n_stokes; ++k) {
                        if (derivs_layers[i_layer][j])
                             a += P_l[i_layer][j][k] * P_a[i_layer][k];
                        a += I_l[j][k] * I_a[k];
                    }

                    if (i_input == 0)
                         I_ss_l[j][ii_input] = a;
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
void single_scattered_radiance_dn_tl_with_ad(int n_stokes, int n_derivs, int n_layers, double F_0, int n_ulevels, int *ulevels, double *utaus, double *umus, int n_umus, double *omega, double **omega_l, double *ltau, double **ltau_l, double *btran, double **btran_l, double *as_0, double **as_0_l, double *atran, double **atran_l, double **P, double ***P_l, double *I_in, double **I_in_l, double **I_ss, double ***I_ss_l, int utau_output, uchar **derivs_layers, uchar **derivs_beam_down, save_tree_data save_tree, work_data work) {

     int i;
     int j;
     int k;
     int l;

     int i_input;
     int n_inputs;

     int ii_input;
     int jj_input;

     double ***a;

     double *omega_a;
     double *ltau_a;
     double *btran_a;
     double *as_0_a;
     double **P_a;
     double *I_in_a;

     double **I_ss_a;


     save_tree_decode_s(&save_tree, "single_scattered_radiance_dn");


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     omega_a = get_work_d1(&work, n_layers);
     ltau_a  = get_work_d1(&work, n_layers);
     btran_a = get_work_d1(&work, n_layers + 1);
     as_0_a  = get_work_d1(&work, n_layers);
     P_a     = get_work_d2(&work, n_layers, n_umus * n_stokes);
     I_in_a  = get_work_d1(&work,           n_umus * n_stokes);

     I_ss_a  = get_work_d2(&work, n_ulevels, n_umus * n_stokes);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     init_array2_d(I_ss_a, n_ulevels, n_umus * n_stokes, 0.);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     n_inputs = n_ulevels * n_umus * n_stokes;

     for (i = 0; i < n_inputs; ++i) {
          if (i < n_ulevels * n_umus * n_stokes) {
               i_input = 0;
               ii_input = i / (n_umus * n_stokes);
               jj_input = i % (n_umus * n_stokes);
          }


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          init_array1_d(omega_a, n_layers,     0.);
          init_array1_d(ltau_a,  n_layers,     0.);
          init_array1_d(btran_a, n_layers + 1, 0.);
          init_array1_d(as_0_a,  n_layers,     0.);
          init_array2_d(P_a,     n_layers, n_umus * n_stokes, 0.);
          init_array1_d(I_in_a,            n_umus * n_stokes, 0.);


          if (i_input == 0)
               I_ss_a[ii_input][jj_input] = 1.;


          single_scattered_radiance_dn_a(n_stokes, n_layers, F_0, n_ulevels, ulevels, utaus, umus, n_umus, omega, omega_a, ltau, ltau_a, btran, btran_a, as_0, as_0_a, atran, NULL, P, P_a, I_in, I_in_a, I_ss, I_ss_a, utau_output, save_tree, work);


          /*--------------------------------------------------------------------
           *
           *------------------------------------------------------------------*/
          if (i_input == 0) {
               if (i_input == 0)
                    a = I_ss_l;

               for (j = 0; j < n_derivs; ++j) {
                    a[j][ii_input][jj_input] = 0.;

                    for (k = 0; k < n_layers; ++k) {
                         a[j][ii_input][jj_input] += omega_l [k][j] * omega_a[k];
                         a[j][ii_input][jj_input] += ltau_l  [k][j] * ltau_a [k];
                         a[j][ii_input][jj_input] += btran_l [k][j] * btran_a[k];
                         a[j][ii_input][jj_input] += as_0_l  [k][j] * as_0_a [k];

                         for (l = 0; l < n_umus * n_stokes; ++l) {
                              if (derivs_layers[k][j])
                                  a[j][ii_input][jj_input] += P_l[k][j][l] * P_a[k][l];
                         }
                    }

                    a[j][ii_input][jj_input] += btran_l[k][j] * btran_a[k];

                    for (k = 0; k < n_umus * n_stokes; ++k)
                        a[j][ii_input][jj_input] += I_in_l[j][k] * I_in_a[k];
               }
          }
     }
}
