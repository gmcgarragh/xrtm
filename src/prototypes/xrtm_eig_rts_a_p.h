/* xrtm_eig_rts_a.c */
int forward_save_calc_global_r_and_t_alloc(forward_save_calc_global_r_and_t_data *d, int n_quad);
void forward_save_calc_global_r_and_t_free(forward_save_calc_global_r_and_t_data *d);
void calc_global_r_and_t_a(int n_quad, double ltau, double *ltau_a, double *nu, double **X_p, double **X_m, double **R, double **T, double *nu_a, double **X_p_a, double **X_m_a, double **R_a, double **T_a, int symmetric, save_tree_data save_tree, work_data work);
int forward_save_calc_global_r_and_t_alloc2(forward_save_calc_global_r_and_t_data2 *d, int n_quad);
void forward_save_calc_global_r_and_t_free2(forward_save_calc_global_r_and_t_data2 *d);
void calc_global_r_and_t_a2(int n_quad, double ltau, double *ltau_a, dcomplex *nu, dcomplex **X_p, dcomplex **X_m, double **R, double **T, dcomplex *nu_a, dcomplex **X_p_a, dcomplex **X_m_a, double **R_a, double **T_a, int symmetric, save_tree_data save_tree, work_data work);
int forward_save_rtm_eig_rts_alloc(forward_save_rtm_eig_rts_data *d, int n_quad_v, int flag);
void rtm_eig_rts_a(int n_quad, int n_stokes, double F_0, double *qx_v, double *qw_v, double planck0, double planck1, double omega, double *omega_a, double ltau, double *ltau_a, double as_0, double *as_0_a, double atran, double *atran_a, double *P_0p, double *P_0m, double **r_p, double **t_p, double **r_m, double **t_m, double **R_p, double **T_p, double **R_m, double **T_m, double *S_p, double *S_m, double *P_0p_a, double *P_0m_a, double **r_p_a, double **t_p_a, double **r_m_a, double **t_m_a, double **R_p_a, double **T_p_a, double **R_m_a, double **T_m_a, double *S_p_a, double *S_m_a, int symmetric, int thermal, int vector, int eigen_solver_real, int eigen_solver_complex, uchar derivs_h, uchar derivs_p, save_tree_data save_tree, work_data work);
