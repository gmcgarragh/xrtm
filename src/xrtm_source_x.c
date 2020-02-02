/*******************************************************************************
 *
 ******************************************************************************/
void BUILD_SOURCE_VECTORS_SOLAR_GREEN_S_1N
     (int n_quad, int n_stokes, int n_derivs,
      double *qx_v, double *qw_v, double F_0,
      double omega, double *omega_l, double ltau, double *ltau_l,
      double as_0, double *as_0_l, double atran, double *atran_l,
      double  *P_q0_mm, double  *P_q0_pm,
      TYPE *nu, TYPE **X_p, TYPE **X_m,
      double  *F0_p, double  *F0_m, double  *F1_p, double  *F1_m,
      double **P_q0_mm_l, double **P_q0_pm_l,
      TYPE **nu_l, TYPE ***X_p_l, TYPE ***X_m_l,
      double **F0_p_l, double **F0_m_l, double **F1_p_l, double **F1_m_l,
      uchar *derivs_layers, uchar *derivs_beam,
      save_tree_data save_tree, work_data work) {

     int i;
     int j;
     int k;

     int n_quad_v;

     double a;
#ifdef REAL
     double *v1;
#endif
     double *c;
     double *e;
     double *d;

     TYPE *Lambda;

     TYPE *b_p;
     TYPE *b_m;
     TYPE *f_p;
     TYPE *f_m;
     TYPE *g;
     TYPE *h;
     TYPE *q;
     TYPE *r;

     TYPE *b_p_l;
     TYPE *b_m_l;
     TYPE *d_l;
     TYPE *e_l;
     TYPE *f_p_l;
     TYPE *f_m_l;
     TYPE *g_l;
     TYPE *h_l;
     TYPE *q_l;
     TYPE *r_l;
#ifndef REAL
     dcomplex *z1;
#endif
     n_quad_v = n_quad * n_stokes;


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     Lambda = get_work1(&work, WORK_XX);
     b_p    = get_work1(&work, WORK_XX);
     b_m    = get_work1(&work, WORK_XX);
     c      = get_work1(&work, WORK_DX);
     d      = get_work1(&work, WORK_DX);
     e      = get_work1(&work, WORK_DX);
     f_p    = get_work1(&work, WORK_XX);
     f_m    = get_work1(&work, WORK_XX);
     g      = get_work1(&work, WORK_XX);
     h      = get_work1(&work, WORK_XX);
     q      = get_work1(&work, WORK_XX);
     r      = get_work1(&work, WORK_XX);
#ifndef REAL
     z1     = get_work1(&work, WORK_XX);
#endif

     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     a = F_0 / (4. * PI);

     for (i = 0; i < n_quad_v; ++i)
          Lambda[i] = XEXP(-ltau * nu[i]);

     for (i = 0; i < n_quad_v; ++i) {
          b_m[i] = (Lambda[i] - atran            ) / (as_0 - nu[i]);
          b_p[i] = (1.        - atran * Lambda[i]) / (as_0 + nu[i]);
/*
          double ptau;

          ptau = ltau;
          b_m[i] = (XEXP(-ptau * nu[i]) - XEXP(-ptau * as_0)                               ) / (as_0 - nu[i]);

          ptau = 0.;
          b_p[i] = (XEXP(-ptau * as_0)  - XEXP(-ltau * as_0) * XEXP(-(ltau - ptau) * nu[i])) / (as_0 + nu[i]);
*/
     }

     for (i = 0; i < n_quad_v; ++i)
          c[i] = qx_v[i] * qw_v[i];

     for (i = 0; i < n_quad_v; ++i) {
          d[i] = 0.;
          for (j = 0; j < n_quad_v; ++j)
               d[i] += c[j] * (XCONJ(X_p[j][i]) * X_p[j][i] - XCONJ(X_m[j][i]) * X_m[j][i]);
     }

     for (i = 0; i < n_quad_v; ++i)
          e[i] = a * omega / d[i];

     for (i = 0; i < n_quad_v; ++i) {
          f_m[i] = e[i] * b_m[i];
          f_p[i] = e[i] * b_p[i];
     }

     dm_v_mul_D_A(n_quad, n_stokes, P_q0_mm, P_q0_mm);

     for (i = 0; i < n_quad_v; ++i) {
          g[i] = 0.;
          h[i] = 0.;
          for (j = 0; j < n_quad_v; ++j) {
               g[i] += qw_v[j] * (P_q0_pm[j] * -XCONJ(X_p[j][i]) + P_q0_mm[j] *  XCONJ(X_m[j][i]));
               h[i] += qw_v[j] * (P_q0_mm[j] *  XCONJ(X_p[j][i]) + P_q0_pm[j] * -XCONJ(X_m[j][i]));
          }
     }

     dm_v_mul_D_A(n_quad, n_stokes, P_q0_mm, P_q0_mm);

     for (i = 0; i < n_quad_v; ++i) {
          q[i] = f_m[i] * g[i];
          r[i] = f_p[i] * h[i];
     }
#ifdef REAL
     dm_v_mul(X_m, r, n_quad_v, n_quad_v, F0_p);
     dm_v_mul(X_p, q, n_quad_v, n_quad_v, F1_p);

     dm_v_mul(X_p, r, n_quad_v, n_quad_v, F0_m);
     dvec_scale(-1., F0_m, F0_m, n_quad_v);
     dm_v_mul(X_m, q, n_quad_v, n_quad_v, F1_m);
     dvec_scale(-1., F1_m, F1_m, n_quad_v);
#else
     zm_v_mul(X_m, r, n_quad_v, n_quad_v, z1);
     zm_v_mul(X_p, q, n_quad_v, n_quad_v, z1);

     zm_v_mul(X_p, r, n_quad_v, n_quad_v, z1);
     dvec_scale(-1., F0_m, F0_m, n_quad_v);
     zm_v_mul(X_m, q, n_quad_v, n_quad_v, z1);
     dvec_scale(-1., F1_m, F1_m, n_quad_v);
#endif
     dm_v_mul_D_A(n_quad, n_stokes, F0_m, F0_m);
     dm_v_mul_D_A(n_quad, n_stokes, F1_m, F1_m);


     /*-------------------------------------------------------------------------
      *
      *-----------------------------------------------------------------------*/
     if (flags_or(derivs_layers, n_derivs)) {
#ifdef REAL
          v1    = get_work1(&work, WORK_DX);
#endif
          b_p_l = get_work1(&work, WORK_XX);
          b_m_l = get_work1(&work, WORK_XX);
          e_l   = get_work1(&work, WORK_XX);
          d_l   = get_work1(&work, WORK_XX);
          f_p_l = get_work1(&work, WORK_XX);
          f_m_l = get_work1(&work, WORK_XX);
          g_l   = get_work1(&work, WORK_XX);
          h_l   = get_work1(&work, WORK_XX);
          q_l   = get_work1(&work, WORK_XX);
          r_l   = get_work1(&work, WORK_XX);
     }

     if (flags_or(derivs_beam, n_derivs)) {
          b_p_l = get_work1(&work, WORK_XX);
          b_m_l = get_work1(&work, WORK_XX);
          f_p_l = get_work1(&work, WORK_XX);
          f_m_l = get_work1(&work, WORK_XX);
          q_l   = get_work1(&work, WORK_XX);
          r_l   = get_work1(&work, WORK_XX);
     }

     for (i = 0; i < n_derivs; ++i) {
          if (derivs_layers[i]) {
               for (j = 0; j < n_quad_v; ++j) {
                    TYPE Lambda_l = (-ltau_l[i] * nu[j] - ltau * nu_l[i][j]) * Lambda[j];

                    b_m_l[j] = (Lambda_l - atran_l[i]                                - b_m[j] * (as_0_l[i] - nu_l[i][j])) / (as_0 - nu[j]);
                    b_p_l[j] = (         - atran_l[i] * Lambda[j] - atran * Lambda_l - b_p[j] * (as_0_l[i] + nu_l[i][j])) / (as_0 + nu[j]);
/*
                    double ptau;
                    double ptau_l;

                    ptau   = ltau;
                    ptau_l = ltau_l[i];
                    b_m_l[j] = ((-ptau_l * nu[j] - ptau * nu_l[i][j]) * exp(-ptau * nu[j]) - (-ptau_l * as_0 - ptau * as_0_l[i]) * exp(-ptau * as_0)                                                                                                                                                   - b_m[j] * (as_0_l[i] - nu_l[i][j])) / (as_0 - nu[j]);

                    ptau   = 0.;
                    ptau_l = 0.;
                    b_p_l[j] = ((-ptau_l * as_0 - ptau * as_0_l[i]) * exp(-ptau * as_0)    - (-ltau_l[i] * as_0 - ltau * as_0_l[i]) * exp(-ltau * as_0) * exp(-(ltau - ptau) * nu[j]) - exp(-ltau * as_0) * (-(ltau_l[i] - ptau_l) * nu[j] - (ltau - ptau) * nu_l[i][j]) * exp(-(ltau - ptau) * nu[j]) - b_p[j] * (as_0_l[i] + nu_l[i][j])) / (as_0 + nu[j]);
*/
               }

               for (j = 0; j < n_quad_v; ++j) {
                    d_l[j] = 0.;
                    for (k = 0; k < n_quad_v; ++k)
                         d_l[j] += c[k] * 2. * (X_p_l[i][k][j] * X_p[k][j] - X_m_l[i][k][j] * X_m[k][j]);
               }

               for (j = 0; j < n_quad_v; ++j)
                    e_l[j] = (a * omega_l[i] - e[j] * d_l[j]) / d[j];

               for (j = 0; j < n_quad_v; ++j) {
                    f_p_l[j] = e_l[j] * b_p[j] + e[j] * b_p_l[j];
                    f_m_l[j] = e_l[j] * b_m[j] + e[j] * b_m_l[j];
               }

               for (j = 0; j < n_quad_v; ++j) {
                    g_l[j] = 0.;
                    h_l[j] = 0.;
                    for (k = 0; k < n_quad_v; ++k) {
                         g_l[j] += qw_v[k] * (P_q0_pm_l[i][k] * -X_p[k][j] + P_q0_pm[k] * -X_p_l[i][k][j] + P_q0_mm_l[i][k] *  X_m[k][j] + P_q0_mm[k] *  X_m_l[i][k][j]);
                         h_l[j] += qw_v[k] * (P_q0_mm_l[i][k] *  X_p[k][j] + P_q0_mm[k] *  X_p_l[i][k][j] + P_q0_pm_l[i][k] * -X_m[k][j] + P_q0_pm[k] * -X_m_l[i][k][j]);
                    }
               }

               for (j = 0; j < n_quad_v; ++j) {
                    q_l[j] = f_m_l[j] * g[j] + f_m[j] * g_l[j];
                    r_l[j] = f_p_l[j] * h[j] + f_p[j] * h_l[j];
               }
#ifdef REAL
               dm_v_mul(X_m_l[i], r, n_quad_v, n_quad_v, F0_p_l[i]);
               dm_v_mul(X_m, r_l, n_quad_v, n_quad_v, v1);
               dvec_add(F0_p_l[i], v1, F0_p_l[i], n_quad_v);

               dm_v_mul(X_p_l[i], q, n_quad_v, n_quad_v, F1_p_l[i]);
               dm_v_mul(X_p, q_l, n_quad_v, n_quad_v, v1);
               dvec_add(F1_p_l[i], v1, F1_p_l[i], n_quad_v);

               dm_v_mul(X_p_l[i], r, n_quad_v, n_quad_v, F0_m_l[i]);
               dm_v_mul(X_p, r_l, n_quad_v, n_quad_v, v1);
               dvec_add(F0_m_l[i], v1, F0_m_l[i], n_quad_v);
               dvec_scale(-1., F0_m_l[i], F0_m_l[i], n_quad_v);

               dm_v_mul(X_m_l[i], q, n_quad_v, n_quad_v, F1_m_l[i]);
               dm_v_mul(X_m, q_l, n_quad_v, n_quad_v, v1);
               dvec_add(F1_m_l[i], v1, F1_m_l[i], n_quad_v);
               dvec_scale(-1., F1_m_l[i], F1_m_l[i], n_quad_v);
#else
#endif
          }
          else
          if (derivs_beam[i] && as_0_l[i] != 0.) {
               for (j = 0; j < n_quad_v; ++j) {

                    b_m_l[j] = (- atran_l[i]             - b_m[j] * as_0_l[i]) / (as_0 - nu[j]);
                    b_p_l[j] = (- atran_l[i] * Lambda[j] - b_p[j] * as_0_l[i]) / (as_0 + nu[j]);
/*
                    double ptau;

                    ptau   = ltau;
                    b_m_l[j] = (                                          - (-ptau * as_0_l[i]) * exp(-ptau * as_0)                               - b_m[j] * as_0_l[i]) / (as_0 - nu[j]);

                    ptau   = 0.;
                    b_p_l[j] = (  (-ptau * as_0_l[i]) * exp(-ptau * as_0) - (-ltau * as_0_l[i]) * exp(-ltau * as_0) * exp(-(ltau - ptau) * nu[j]) - b_p[j] * as_0_l[i]) / (as_0 + nu[j]);
*/
               }

               for (j = 0; j < n_quad_v; ++j) {
                    f_p_l[j] = e[j] * b_p_l[j];
                    f_m_l[j] = e[j] * b_m_l[j];
               }

               for (j = 0; j < n_quad_v; ++j) {
                    q_l[j] = f_m_l[j] * g[j];
                    r_l[j] = f_p_l[j] * h[j];
               }
#ifdef REAL
               dm_v_mul(X_m, r_l, n_quad_v, n_quad_v, F0_p_l[i]);

               dm_v_mul(X_p, q_l, n_quad_v, n_quad_v, F1_p_l[i]);

               dm_v_mul(X_p, r_l, n_quad_v, n_quad_v, F0_m_l[i]);
               dvec_scale(-1., F0_m_l[i], F0_m_l[i], n_quad_v);

               dm_v_mul(X_m, q_l, n_quad_v, n_quad_v, F1_m_l[i]);
               dvec_scale(-1., F1_m_l[i], F1_m_l[i], n_quad_v);
#else
#endif
          }
          else
          if (derivs_beam[i] && as_0_l[i] == 0.) {
               dvec_zero(F0_p_l[i], n_quad_v);
               dvec_zero(F0_m_l[i], n_quad_v);

               dvec_zero(F1_p_l[i], n_quad_v);
               dvec_zero(F1_m_l[i], n_quad_v);
          }
     }
}
