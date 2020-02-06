/* xrtm_source.c */
void build_source_vectors_solar_classic_1n(int n_quad, int n_stokes, int n_derivs, double *qx_v, double F_0, double omega, double *omega_l, double as_0, double *as_0_l, double *P_q0_mm, double *P_q0_pm, double **tpr, double **tmr, double **gamma, double *F_p, double *F_m, double **P_q0_mm_l, double **P_q0_pm_l, double ***tpr_l, double ***tmr_l, double ***gamma_l, double **F_p_l, double **F_m_l, uchar *derivs_layers, uchar *derivs_beam, save_tree_data save_tree, work_data work);
void build_source_vectors_solar_classic_2n(int n_quad_v, int n_derivs, double *qx_v, double F_0, double omega, double *omega_l, double as_0, double *as_0_l, double *P_q0_mm, double *P_q0_pm, double **r, double **t, double *F_p, double *F_m, double **P_q0_mm_l, double **P_q0_pm_l, double ***r_l, double ***t_l, double **F_p_l, double **F_m_l, uchar *derivs_layers, uchar *derivs_beam, work_data work);
void build_source_vectors_solar_classic_2n2(int n_quad_v, int n_derivs, double *qx_v, double F_0, double omega, double *omega_l, double as_0, double *as_0_l, double *P_q0_mm, double *P_q0_pm, double **r_p, double **t_p, double **r_m, double **t_m, double *F_p, double *F_m, double **P_q0_mm_l, double **P_q0_pm_l, double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l, double **F_p_l, double **F_m_l, uchar *derivs_layers, uchar *derivs_beam, work_data work);
void build_thermal_source_coefficients(int n_quad, int n_stokes, int n_derivs, double *qx_v, double planck0, double planck1, double *planck0_l, double *planck1_l, double omega, double *omega_l, double ltau, double *ltau_l, double **r_p, double **t_p, double **r_m, double **t_m, double **At_p, double **At_m, double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l, double ***At_p_l, double ***At_m_l, uchar *derivs_layers, uchar *derivs_thermal, work_data work);
void build_source_vectors_thermal(int n_quad, int n_stokes, int n_derivs, double ltau, double *ltau_l, double **At_p, double **At_m, double *F0_p, double *F0_m, double *F1_p, double *F1_m, double ***At_p_l, double ***At_m_l, double **F0_p_l, double **F0_m_l, double **F1_p_l, double **F1_m_l, uchar *derivs_layers, uchar *derivs_thermal, work_data work);
void build_source_vectors_thermal2(int n_quad, int n_stokes, int n_derivs, double *qx_v, double planck0, double planck1, double *planck0_l, double *planck1_l, double omega, double *omega_l, double ltau, double *ltau_l, double **r_p, double **t_p, double **r_m, double **t_m, double *F0_p, double *F0_m, double *F1_p, double *F1_m, double ***r_p_l, double ***t_p_l, double ***r_m_l, double ***t_m_l, double **F0_p_l, double **F0_m_l, double **F1_p_l, double **F1_m_l, uchar *derivs_layers, uchar *derivs_thermal, work_data work);
void build_global_source_solar(int n_quad_v, int n_derivs, double atran, double *atran_l, double **R, double **T, double *F_p, double *F_m, double *S_p, double *S_m, double ***R_l, double ***T_l, double **F_p_l, double **F_m_l, double **S_p_l, double **S_m_l, uchar *derivs_layers, uchar *derivs_beam, work_data work, double *F1_p, double *F1_m, double **F1_p_l, double **F1_m_l);
void build_global_source_solar2(int n_quad_v, int n_derivs, double atran, double *atran_l, double **R_p, double **T_p, double **R_m, double **T_m, double *F_p, double *F_m, double *S_p, double *S_m, double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l, double **F_p_l, double **F_m_l, double **S_p_l, double **S_m_l, uchar *derivs_layers, uchar *derivs_beam, work_data work, double *F1_p, double *F1_m, double **F1_p_l, double **F1_m_l);
void build_global_source_thermal(int n_quad_v, int n_derivs, double **R_p, double **T_p, double **R_m, double **T_m, double *F0_p, double *F0_m, double *F1_p, double *F1_m, double *S_p, double *S_m, double ***R_p_l, double ***T_p_l, double ***R_m_l, double ***T_m_l, double **F0_p_l, double **F0_m_l, double **F1_p_l, double **F1_m_l, double **S_p_l, double **S_m_l, uchar *derivs_layers, uchar *derivs_thermal, work_data work);
void scale_source_vectors_solar(int n_quad_v, int n_derivs, double btran, double *btran_l, double *F_p1, double *F_m1, double *F_p2, double *F_m2, double **F_p_l1, double **F_m_l1, double **F_p_l2, double **F_m_l2, uchar *derivs_layers, uchar *derivs_beam, work_data work);
void scale_add_source_vectors(int n_quad_v, int n_derivs, double btran, double *btran_l, double *F_p1, double *F_m1, double *Fl_p1, double *Fl_m1, double *F_p2, double *F_m2, double **F_p_l1, double **F_m_l1, double **Fl_p_l1, double **Fl_m_l1, double **F_p_l2, double **F_m_l2, int solar, int thermal, uchar *derivs_layers, uchar *derivs_beam, uchar *derivs_thermal, uchar *derivs_sources, work_data work);
void build_source_vectors_solar_green_s_1n(int n_quad, int n_stokes, int n_derivs, double *qx_v, double *qw_v, double F_0, double omega, double *omega_l, double ltau, double *ltau_l, double as_0, double *as_0_l, double atran, double *atran_l, double *P_q0_mm, double *P_q0_pm, double *nu, double **X_p, double **X_m, double *F0_p, double *F0_m, double *F1_p, double *F1_m, double **P_q0_mm_l, double **P_q0_pm_l, double **nu_l, double ***X_p_l, double ***X_m_l, double **F0_p_l, double **F0_m_l, double **F1_p_l, double **F1_m_l, uchar *derivs_layers, uchar *derivs_beam, save_tree_data save_tree, work_data work);
void build_source_vectors_solar_green_s_1n2(int n_quad, int n_stokes, int n_derivs, double *qx_v, double *qw_v, double F_0, double omega, double *omega_l, double ltau, double *ltau_l, double as_0, double *as_0_l, double atran, double *atran_l, double *P_q0_mm, double *P_q0_pm, dcomplex *nu, dcomplex **X_p, dcomplex **X_m, double *F0_p, double *F0_m, double *F1_p, double *F1_m, double **P_q0_mm_l, double **P_q0_pm_l, dcomplex **nu_l, dcomplex ***X_p_l, dcomplex ***X_m_l, double **F0_p_l, double **F0_m_l, double **F1_p_l, double **F1_m_l, uchar *derivs_layers, uchar *derivs_beam, save_tree_data save_tree, work_data work);