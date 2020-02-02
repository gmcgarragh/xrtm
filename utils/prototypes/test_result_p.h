/* test_result.c */
int test_result_init(test_result_data *d, int n_solvers);
void test_result_free(test_result_data *d);
int test_cmd_init(test_cmd_data *d, int n_solvers);
void test_cmd_free(test_cmd_data *d);
int save_result_header(test_result_data *d, int index, int index2, int i_out_level, int i_deriv, double out_theta, double out_phi, int i_stokes);
int save_result_solver_mask(test_result_data *d, int i_solver, enum xrtm_solver_mask mask);
int save_result_solver_values(test_result_data *d, int i_solver, double x_p, double x_m);
int save_result_solver_deltas(test_result_data *d, int i_solver, double x_p, double x_m);
