/* xrtm_model.c */
int xrtm_phase_type_to_coef_type(xrtm_data *d);
int solar_surface_for_i_four(xrtm_data *d, int i_four);
int thermal_surface_for_i_four(xrtm_data *d, int i_four);
void get_delta_m_f(xrtm_data *d, int i_layer, int i_deriv, int n_derivs, double a, double *f, double *f_l);
void get_chapman(xrtm_data *d, double ***chapman);
void get_Y_p(xrtm_data *d, int i_four, double ***Y_p);
void get_Y_0(xrtm_data *d, int i_four, double **Y_0);
void get_Y_u(xrtm_data *d, int i_four, double ***Y_u);
int get_solution(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_phis, double **phis, double ****I_p, double ****I_m, double *****I_p_l, double *****I_m_l, double *mean_p, double *mean_m, double **mean_p_l, double **mean_m_l, double *flux_p, double *flux_m, double **flux_p_l, double **flux_m_l, double *flux_div, double **flux_div_l, save_tree_data save_tree, work_data work);
