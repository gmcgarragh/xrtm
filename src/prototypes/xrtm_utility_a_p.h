/* xrtm_utility_a.c */
void delta_m_coef_a(int n_coef, double f, double *f_a, double **coef, double **coef_a, double **coef_prime_a, int allsix, int coef_type);
void delta_m_omega_a(double f, double *f_a, double omega, double *omega_a, double *omega_prime_a);
void delta_m_ltau_a(double f, double *f_a, double omega, double *omega_a, double ltau, double *ltau_a, double *ltau_prime_a);
void delta_m_a(int n_coef, double f, double *f_a, double **coef, double omega, double ltau, double **coef_a, double **coef_prime_a, double *omega_a, double *omega_prime_a, double *ltau_a, double *ltau_prime_a, int allsix, int coef_type);
void n_t_tms_scaling_a(double f, double *f_a, double omega0, double omega_tms, double *omega0_a, double omega_tms_a);
void build_local_r_and_t_a(int i_four, int n_quad, double *qx_v, double *qw_v, double omega, double *omega_a, double **P_p, double **P_m, double **r, double **t, double **P_p_a, double **P_m_a, double **r_a, double **t_a, work_data work);
void build_txr_a(int n_quad, int n_stokes, double **r_p, double **t_p, double **tpr, double **tmr, double **r_p_a, double **t_p_a, double **tpr_a, double **tmr_a, work_data work);
void build_gamma_a(int n_quad, double **tpr, double **tmr, double **gamma, double **tpr_a, double **tmr_a, double **gamma_a, work_data work);
void build_local_r_and_t_tl_with_ad(int i_four, int n_quad, int n_derivs, double *qx_v, double *qw_v, double omega, double *omega_l, double **P_p, double **P_m, double **r, double **t, double ***P_p_l, double ***P_m_l, double ***r_l, double ***t_l, uchar *derivs, work_data work);
void build_txr_tl_with_ad(int n_quad, int n_stokes, int n_derivs, double **r_p, double **t_p, double **tpr, double **tmr, double ***r_p_l, double ***t_p_l, double ***tpr_l, double ***tmr_l, uchar *derivs, work_data work);
void build_gamma_tl_with_ad(int n_quad, int n_derivs, double **tpr, double **tmr, double **gamma, double ***tpr_l, double ***tmr_l, double ***gamma_l, uchar *derivs, work_data work);
