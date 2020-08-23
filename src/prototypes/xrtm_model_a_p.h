/* xrtm_model_a.c */
void init_opt_props0_all_a(xrtm_data *d, int n_coef, work_data work);
void init_opt_props_all_a(xrtm_data *d, int n_coef, work_data work);
void init_beam_params_all_a(xrtm_data *d, work_data work);
int forward_save_get_local_r_t_u_w_alloc(forward_save_get_local_r_t_u_w_data *d, int options, int n_quad_v_x);
int forward_save_get_layer_R_T_S_U_W_V_alloc(forward_save_get_layer_R_T_S_U_W_V_data *d, int options, int n_quad_v_x);
int forward_save_get_total_TOA_R_S_U_V_alloc(forward_save_get_total_TOA_R_S_U_V_data *d, int n_quad_v_x, int n_layers);
int forward_save_get_total_R_T_S_U_W_V_alloc(forward_save_get_total_R_T_S_U_W_V_data *d, int n_quad_v_x, int n_layers);
void init_diff_bound_input_all_a(xrtm_data *d, int n_four, work_data work);
int forward_save_update_diff_bound_input_alloc(forward_save_update_diff_bound_input_data *d, int n_quad_v, int n_umus_v);
int forward_save_get_single_alloc(forward_save_get_single_data *d, int options, int n_coef, int n_quad_v, int n_umus_v, int n_mus2, int n_phis, int n_layers, int n_stokes);
int forward_save_fourier_get_add_bottom_up_alloc(forward_save_fourier_get_add_bottom_up_data *d, int n_quad_v_x);
int forward_save_fourier_get_add_both_ways_alloc(forward_save_fourier_get_add_both_ways_data *d, int n_quad_v_x);
int forward_save_fourier_get_bvp_alloc(forward_save_fourier_get_bvp_data *d, int options, int solver, int surface, int n_layers, int n_quad_v_x, int n_umus_v);
int forward_save_apply_corrections_radiance_alloc(forward_save_apply_corrections_radiance_data *d, int options, int n_coef, int n_quad_v_x, int n_mus2, int n_phis, int n_layers, int n_stokes);
int forward_save_apply_corrections_fourier_alloc(forward_save_apply_corrections_fourier_data *d, int n_quad_v, int n_quad_v_x, int n_umus_v, int n_layers);
int forward_save_get_solution_internal_alloc(forward_save_get_solution_internal_data *d, int n_four, int n_ulevels, int n_quad_v_x);
int get_solution_internal_a(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double ****I_p_a, double ****I_m_a, double *mean_p, double *mean_m, double *mean_p_a, double *mean_m_a, double *flux_p, double *flux_m, double *flux_p_a, double *flux_m_a, double *flux_div, double *flux_div_a, save_tree_data save_tree, work_data work);
int xrtm_solution_a(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double ****I_p_a, double ****I_m_a, double *mean_p, double *mean_m, double *mean_p_a, double *mean_m_a, double *flux_p, double *flux_m, double *flux_p_a, double *flux_m_a, double *flux_div, double *flux_div_a);
