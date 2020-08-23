/* xrtm_eig_bvp_a.c */
int forward_save_solve_bvp_alloc(forward_save_solve_bvp_data *d, int n_layers, int n_quad_v, int m_comp, int n_comp);
void forward_save_solve_bvp_free(forward_save_solve_bvp_data *d);
int forward_save_calc_radiance_levels_alloc(forward_save_calc_radiance_levels_data *d, int n_layers, int n_quad_v);
void forward_save_calc_radiance_levels_free(forward_save_calc_radiance_levels_data *d);
int forward_save_solve_bvp_alloc2(forward_save_solve_bvp_data2 *d, int n_layers, int n_quad_v, int m_comp, int n_comp);
void forward_save_solve_bvp_free2(forward_save_solve_bvp_data2 *d);
int forward_save_calc_radiance_levels_alloc2(forward_save_calc_radiance_levels_data2 *d, int n_layers, int n_quad_v);
void forward_save_calc_radiance_levels_free2(forward_save_calc_radiance_levels_data2 *d);
int forward_save_rtm_eig_bvp_alloc(forward_save_rtm_eig_bvp_data *d, int n_layers, int n_quad_v, int n_comp, int thermal, int vector);
void rtm_eig_bvp_a(int i_four, int n_quad, int n_stokes, int n_layers, double qf, double *qx_v, double *qw_v, double F_0, double mu_0, int *ulevels, double *utaus, int n_ulevels, double *umus, int n_umus, double *planck, double *omega, double *omega_a, double *ltau, double *ltau_a, double **Rs_qq, double **Rs_qq_a, double *Rs_u0, double *Rs_u0_a, double **Rs_uq, double **Rs_uq_a, double *btran, double *btran_a, double *as_0, double *as_0_a, double *atran, double *atran_a, double **P_q0_mm, double **P_q0_pm, double **P_u0_mm, double **P_u0_pm, double ***P_uq_pp, double ***P_uq_mp, double ***P_uq_mm, double ***P_uq_pm, double ***r_p, double ***t_p, double ***r_m, double ***t_m, double **P_q0_mm_a, double **P_q0_pm_a, double **P_u0_mm_a, double **P_u0_pm_a, double ***P_uq_pp_a, double ***P_uq_mp_a, double ***P_uq_mm_a, double ***P_uq_pm_a, double ***r_p_a, double ***t_p_a, double ***r_m_a, double ***t_m_a, double *I1_m, double *I1_m_a, double *In_p, double *In_p_a, double **I_p, double **I_m, double **I_p_a, double **I_m_a, int sfi, int surface, int thermal, int upwelling, int downwelling, int utau_output, int vector, uchar *derivs_layers, uchar *derivs_beam, save_tree_data save_tree, work_data work);
