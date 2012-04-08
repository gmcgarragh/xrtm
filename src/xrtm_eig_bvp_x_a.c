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
void FORWARD_SAVE_SOLVE_BVP_FREE(FORWARD_SAVE_SOLVE_BVP_DATA *d);

int FORWARD_SAVE_SOLVE_BVP_ALLOC(FORWARD_SAVE_SOLVE_BVP_DATA *d, int n_layers, int n_quad_v, int m_comp, int n_comp) {

     d->free = (void (*)(void *)) FORWARD_SAVE_SOLVE_BVP_FREE;

     d->ipiv   = alloc_array1_i(n_comp);

     d->lambda = XCAT(alloc_array2_, TYPE_POSTFIX)(n_layers, n_quad_v);

     d->A      = XCAT(alloc_array2_, TYPE_POSTFIX)(n_comp, m_comp);
     d->B      = XCAT(alloc_array1_, TYPE_POSTFIX)(n_comp);
/*
     d->c      = XCAT(alloc_array2_, TYPE_POSTFIX)(n_quad_v, n_quad_v);
*/
     return 0;
}



void FORWARD_SAVE_SOLVE_BVP_FREE(FORWARD_SAVE_SOLVE_BVP_DATA *d) {

     free_array1_i(d->ipiv);

     XCAT(free_array2_, TYPE_POSTFIX)(d->lambda);

     XCAT(free_array2_, TYPE_POSTFIX)(d->A);
     XCAT(free_array1_, TYPE_POSTFIX)(d->B);
/*
     XCAT(free_array2_, TYPE_POSTFIX)(d->c);
*/
}



/*******************************************************************************
 *
 ******************************************************************************/
static void SOLVE_BVP_A(int n_quad, int n_stokes, int n_layers,
            double *ltau, double *ltau_a,
            double **Rs_qq, double **Rs_qq_a, 
            double *atran, double *atran_a,
            TYPE **nu, TYPE ***X_p, TYPE ***X_m,
            double **F_p, double **F_m,
            double **F0_p, double **F0_m, double **F1_p, double **F1_m,
            TYPE **nu_a, TYPE ***X_p_a, TYPE ***X_m_a,
            double **F_p_a, double **F_m_a,
            double **F0_p_a, double **F0_m_a, double **F1_p_a, double **F1_m_a,
            TYPE *B, TYPE *B_a,
            double *I1_m, double *I1_m_a, double *In_p, double *In_p_a,
            int surface, int thermal, uchar *derivs_h, uchar *derivs_p,
            save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int iii;
     int j;
     int k;
     int l;

     int n_quad_v;
     int n_quad_v2;
     int n_quad_v3;

     int n_diags;

     int m_comp;
     int n_comp;

     int info;
     int nrhs = 1;

     TYPE a;
     TYPE b;
     TYPE c;
     TYPE d;
     TYPE e;

     TYPE **lambda_a;

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
     save_tree_encode_s(&save_tree, SOLVE_BVP_SAVE_TREE_STRING);

     save_tree_retrieve_data(&save_tree, FORWARD_SAVE_SOLVE_BVP_DATA, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     lambda_a = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     xgbtrs_("T", &n_comp, &n_diags, &n_diags, &nrhs, *(save->A), &m_comp,
             save->ipiv, B_a, &n_comp, &info);
     if (info) {
          eprintf("ERROR: xgbtrs() info = %d\n", info);
          exit(1);
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          if (derivs_h[i]) {
               for (j = 0; j < n_quad_v; ++j) {
                    lambda_a[i][j] = 0.;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_quad_v; ++i) {
          a = B_a[i];

          if (derivs_p[0])
               F_m_a[0][i] -= XREAL(a);
          if (derivs_h[0]) {
               for (j = 0; j < n_quad_v; ++j) {
                    X_m_a[0][i][j] += a * save->B[j];

                    lambda_a[0][j] += a * save->B[j + n_quad_v] * X_p[0][i][j];

                    X_p_a[0][i][j] += a * save->B[j + n_quad_v] * save->lambda[0][j];
               }
          }
     }

     iii = n_quad_v;
     for (i = 0; i < n_layers - 1; ++i) {
          ii = i * n_quad_v2;
          for (j = 0; j < n_quad_v; ++j) {
               a = B_a[iii++];

               if (derivs_p[i  ]) {
                    F_p_a[i][j] -= XREAL(a) * atran[i];

                    atran_a [i] -= XREAL(a) * F_p[i][j];
               }
               if (derivs_p[i+1])
                    F_p_a[i+1][j] += XREAL(a);

               for (k = 0; k < n_quad_v; ++k) {
                    if (derivs_h[i  ]) {
                         lambda_a[i][k] -= a * save->B[ii + k] * X_p[i][j][k];
                         X_p_a[i][j][k] -= a * save->B[ii + k] * save->lambda[i][k];
                         X_m_a[i][j][k] -= a * save->B[ii + k + n_quad_v];
                    }
                    if (derivs_h[i+1]) {
                         lambda_a[i+1][k] += a * save->B[ii + k + n_quad_v3] * X_m[i+1][j][k];
                         X_m_a[i+1][j][k] += a * save->B[ii + k + n_quad_v3] * save->lambda[i+1][k];
                         X_p_a[i+1][j][k] += a * save->B[ii + k + n_quad_v2];
                    }
               }
          }

          for (j = 0; j < n_quad_v; ++j) {
               a = B_a[iii++];

               if (derivs_p[i  ]) {
                    F_m_a[i][j] -= XREAL(a) * atran[i];

                    atran_a [i] -= XREAL(a) * F_m[i][j];
               }
               if (derivs_p[i+1])
                    F_m_a[i+1][j] += XREAL(a);

               for (k = 0; k < n_quad_v; ++k) {
                    if (derivs_h[i  ]) {
                         lambda_a[i][k] += a * save->B[ii + k] * X_m[i][j][k];
                         X_m_a[i][j][k] += a * save->B[ii + k] * save->lambda[i][k];
                         X_p_a[i][j][k] += a * save->B[ii + k + n_quad_v];
                    }
                    if (derivs_h[i+1]) {
                         lambda_a[i+1][k] -= a * save->B[ii + k + n_quad_v3] * X_p[i+1][j][k];
                         X_p_a[i+1][j][k] -= a * save->B[ii + k + n_quad_v3] * save->lambda[i+1][k];
                         X_m_a[i+1][j][k] -= a * save->B[ii + k + n_quad_v2];
                    }
               }
          }
     }


     i  = n_layers - 1;

     ii = n_comp - n_quad_v2;

     for (j = 0; j < n_quad_v; ++j) {
          b = B_a[n_comp - n_quad_v + j];

          if (derivs_p[i]) {
               F_p_a[i][j] -= XREAL(b) * atran[i];

               atran_a [i] -= XREAL(b) * F_p[i][j];
          }

          if (surface) {
               c = b;

               for (k = 0; k < n_quad_v; ++k) {
                    if (derivs_p[i  ])  {
                         F_m_a[i][k] += XREAL(c) * Rs_qq[j][k] * atran[i];
                         atran_a[i]  += XREAL(c) * Rs_qq[j][k] * F_m[i][k];
                    }
                    if (derivs_h[i+1])
                         Rs_qq_a[j][k] += XREAL(c) * F_m[i][k] * atran[i];
               }

               if (derivs_p[i+1])
                    In_p_a[j] += XREAL(b);
          }

          for (k = 0; k < n_quad_v; ++k) {
               if (derivs_h[i]) {
                    lambda_a[i][k] -= b * save->B[ii + k] * X_p[i][j][k];
                    X_p_a[i][j][k] -= b * save->B[ii + k] * save->lambda[i][k];
                    X_m_a[i][j][k] -= b * save->B[ii + k + n_quad_v];
               }

               if (surface) {
                   c = 0.;
                   for (l = 0; l < n_quad_v; ++l)
                         c += Rs_qq[j][l] * -X_m[i][l][k];

                    if (derivs_h[i])
                         lambda_a[i][k] -= b * B[ii + k] * -      c;
/*
                         lambda_a[i][k] -= b * B[ii + k] * -save->c[j][k];
*/
                    d = b * B[ii + k] * save->lambda[i][k];
                    e = b * B[ii + k + n_quad_v];

                    for (l = 0; l < n_quad_v; ++l) {
                         if (derivs_h[i  ]) {
                              X_m_a[i][l][k] -= Rs_qq[j][l] * d;
                              X_p_a[i][l][k] -= Rs_qq[j][l] * e;
                         }
                         if (derivs_h[i+1]) {
                              Rs_qq_a[j][l] -= XREAL(X_m[i][l][k] * d);
                              Rs_qq_a[j][l] -= XREAL(X_p[i][l][k] * e);
                         }
                    }
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          if (derivs_h[i]) {
               for (j = 0; j < n_quad_v; ++j) {
                    nu_a[i][j] -=       lambda_a[i][j] * ltau[i] * save->lambda[i][j];

                    ltau_a[i]  -= XREAL(lambda_a[i][j] * nu[i][j] * save->lambda[i][j]);
               }
          }
     }


     XCAT(init_array1_, TYPE_POSTFIX)(B_a, n_comp, 0.);
}



/*******************************************************************************
 *
 ******************************************************************************/
void FORWARD_SAVE_CALC_RADIANCE_LEVELS_FREE(FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA *d);

int FORWARD_SAVE_CALC_RADIANCE_LEVELS_ALLOC(FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA *d, int n_layers, int n_quad_v) {

     d->free = (void (*)(void *)) FORWARD_SAVE_CALC_RADIANCE_LEVELS_FREE;

     d->lambda = XCAT(alloc_array2_, TYPE_POSTFIX)(n_layers, n_quad_v);

     return 0;
}



void FORWARD_SAVE_CALC_RADIANCE_LEVELS_FREE(FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA *d) {

     XCAT(free_array2_, TYPE_POSTFIX)(d->lambda);
}



/*******************************************************************************
 *
 ******************************************************************************/
static void CALC_RADIANCE_LEVELS_A(int n_quad, int n_layers, int n_ulevels,
                       int *ulevels, double *ltau, double *ltau_a,
                       double *atran, double *atran_a,
                       TYPE **nu, TYPE ***X_p, TYPE ***X_m,
                       double **F_p, double **F_m,
                       double **F0_p, double **F0_m, double **F1_p, double **F1_m,
                       TYPE **nu_a, TYPE ***X_p_a, TYPE ***X_m_a,
                       double **F_p_a, double **F_m_a,
                       double **F0_p_a, double **F0_m_a, double **F1_p_a, double **F1_m_a,
                       TYPE *B, TYPE *B_a,
                       double **I_p, double **I_m, double **I_p_a, double **I_m_a,
                       uchar *derivs_h, uchar *derivs_p, int thermal,
                       save_tree_data save_tree, work_data work) {

     int i;
     int ii;
     int iii;
     int j;
     int k;

     int n_quad2;

     TYPE a;

     TYPE **lambda_a;

     FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA *save;


     n_quad2 = n_quad * 2;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     save_tree_encode_s(&save_tree, CALC_RADIANCE_LEVELS_SAVE_TREE_STRING);

     save_tree_retrieve_data(&save_tree, FORWARD_SAVE_CALC_RADIANCE_LEVELS_DATA, &save);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     lambda_a = get_work2(&work, WORK_XX, WORK_LAYERS_V, NULL);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          if (derivs_h[i]) {
               for (j = 0; j < n_quad; ++j) {
                    lambda_a[i][j] = 0.;
               }
          }
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_ulevels; ++i) {
          ii = ulevels[i];
if (ii != n_layers) {
          iii = ii * n_quad2;

          for (j = 0; j < n_quad; ++j) {
               a = I_p_a[i][j];

               for (k = 0; k < n_quad; ++k) {
                    B_a[iii + k         ] += a * X_p[ii][j][k];
                    B_a[iii + k + n_quad] += a * X_m[ii][j][k] * save->lambda[ii][k];

                    if (derivs_h[ii]) {
                         X_p_a[ii][j][k] += a * B[iii + k         ];
                         X_m_a[ii][j][k] += a * B[iii + k + n_quad] * save->lambda[ii][k];
                         lambda_a[ii][k] += a * B[iii + k + n_quad] * X_m[ii][j][k];
                    }
               }

               if (derivs_p[ii])
                    F_p_a[ii][j] += I_p_a[i][j];
          }

          for (j = 0; j < n_quad; ++j) {
               a = I_m_a[i][j];

               for (k = 0; k < n_quad; ++k) {
                    B_a[iii + k         ] -= a * X_m[ii][j][k];
                    B_a[iii + k + n_quad] -= a * X_p[ii][j][k] * save->lambda[ii][k];

                    if (derivs_h[ii]) {
                         X_m_a[ii][j][k] -= a * B[iii + k         ];
                         X_p_a[ii][j][k] -= a * B[iii + k + n_quad] * save->lambda[ii][k];
                         lambda_a[ii][k] -= a * B[iii + k + n_quad] * X_p[ii][j][k];
                    }
               }

               if (derivs_p[ii])
                    F_m_a[ii][j] += I_m_a[i][j];
          }
}
else {
          ii = n_layers - 1;

          iii = ii * n_quad2;

          for (j = 0; j < n_quad; ++j) {
               a = I_p_a[i][j];

               for (k = 0; k < n_quad; ++k) {
                    B_a[iii + k         ] += a * X_p[ii][j][k] * save->lambda[ii][k];
                    B_a[iii + k + n_quad] += a * X_m[ii][j][k];

                    if (derivs_h[ii]) {
                         X_p_a[ii][j][k] += a * B[iii + k         ] * save->lambda[ii][k];
                         lambda_a[ii][k] += a * B[iii + k         ] * X_p[ii][j][k];
                         X_m_a[ii][j][k] += a * B[iii + k + n_quad];
                    }
               }

               if (derivs_p[ii]) {
                    F_p_a[ii][j] += I_p_a[i][j] * atran[ii];
                    atran_a[ii]  += I_p_a[i][j] * F_p[ii][j];
               }
          }

          for (j = 0; j < n_quad; ++j) {
               a = I_m_a[i][j];

               for (k = 0; k < n_quad; ++k) {
                    B_a[iii + k         ] -= a * X_m[ii][j][k] * save->lambda[ii][k];
                    B_a[iii + k + n_quad] -= a * X_p[ii][j][k];

                    if (derivs_h[ii]) {
                         X_m_a[ii][j][k] -= a * B[iii + k         ] * save->lambda[ii][k];
                         lambda_a[ii][k] -= a * B[iii + k         ] * X_m[ii][j][k];
                         X_p_a[ii][j][k] -= a * B[iii + k + n_quad];
                    }
               }

               if (derivs_p[ii]) {
                    F_m_a[ii][j] += I_m_a[i][j] * atran[ii];
                    atran_a[ii]  += I_m_a[i][j] * F_m[ii][j];
               }
          }
}
     }


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     for (i = 0; i < n_layers; ++i) {
          if (derivs_h[i]) {
               for (j = 0; j < n_quad; ++j) {
                    nu_a[i][j] -=       lambda_a[i][j] * ltau[i]  * save->lambda[i][j];
                    ltau_a[i]  -= XREAL(lambda_a[i][j] * nu[i][j] * save->lambda[i][j]);
               }
          }
     }
}



/*******************************************************************************
 *
 ******************************************************************************/
static void CALC_RADIANCE_TAUS_A(int n_quad, int n_layers, int n_ulevels,
                     int *ulevels, double *utaus, double *ltau, double *ltau_a,
                     double *as_0, double *as_0_a, double *atran, double *atran_a,
                     TYPE **nu, TYPE ***X_p, TYPE ***X_m,
                     double **F_p, double **F_m,
                     TYPE **nu_a, TYPE ***X_p_a, TYPE ***X_m_a,
                     double **F_p_a, double **F_m_a,
                     TYPE *B, TYPE *B_aa,
                     double **I_p, double **I_m, double **I_p_a, double **I_m_a,
                     uchar *derivs_h, uchar *derivs_p,
                     save_tree_data save_tree, work_data work) {

}



/*******************************************************************************
 *
 ******************************************************************************/
#ifdef INCLUDE_DEV_SOURCE
#include "../src2/xrtm_eig_bvp_x_a2.c"
#endif
