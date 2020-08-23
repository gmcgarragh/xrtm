/* xrtm_eig_util_a.c */
int forward_save_eig_1n_gen_real_alloc(forward_save_eig_1n_gen_real_data *d, int n_quad);
void eig_1n_gen_real_a(int n_quad, double **gamma, double *evals, double **evecs, double **gamma_a, double *evals_a, double **evecs_a, save_tree_data save_tree, work_data work);
int forward_save_eig_1n_gen_complex_alloc(forward_save_eig_1n_gen_complex_data *d, int n_quad);
void eig_1n_gen_complex_a(int n_quad, double **gamma, double *evals_r, double *evals_i, double **evecs, double **gamma_a, double *evals_r_a, double *evals_i_a, double **evecs_a, save_tree_data save_tree, work_data work);
int forward_save_eig_1n_to_2n_real_alloc(forward_save_eig_1n_to_2n_real_data *d, int n_quad);
void eig_1n_to_2n_real_a(int n_quad, double **tpr, double **tmr, double *evals, double **evecs, double *nu, double **X_p, double **X_m, double **tpr_a, double **tmr_a, double *evals_a, double **evecs_a, double *nu_a, double **X_p_a, double **X_m_a, save_tree_data save_tree, work_data work);
int forward_save_eig_1n_to_2n_complex_alloc(forward_save_eig_1n_to_2n_complex_data *d, int n_quad);
void eig_1n_to_2n_complex_a(int n_quad, double **tpr, double **tmr, double *evals_r, double *evals_i, double **evecs, dcomplex *nu, dcomplex **X_p, dcomplex **X_m, double **tpr_a, double **tmr_a, double *evals_r_a, double *evals_i_a, double **evecs_a, dcomplex *nu_a, dcomplex **X_p_a, dcomplex **X_m_a, save_tree_data save_tree, work_data work);
int forward_save_eig_2n_gen_real_alloc(forward_save_eig_2n_gen_real_data *d, int n_quad);
void eig_2n_gen_real_a(int n_quad, double **tpr, double **tmr, double **gamma, double *nu, double **X_p, double **X_m, double **tpr_a, double **tmr_a, double **gamma_a, double *nu_a, double **X_p_a, double **X_m_a, save_tree_data save_tree, work_data work);
int forward_save_eig_2n_gen_complex_alloc(forward_save_eig_2n_gen_complex_data *d, int n_quad);
void eig_2n_gen_complex_a(int n_quad, double **tpr, double **tmr, double **gamma, dcomplex *nu, dcomplex **X_p, dcomplex **X_m, double **tpr_a, double **tmr_a, double **gamma_a, dcomplex *nu_a, dcomplex **X_p_a, dcomplex **X_m_a, save_tree_data save_tree, work_data work);
