/* xrtm_scatter.c */
double scat_angle(double mu_0, double phi_0, double mu, double phi);
void phase_func(int n_coef, double *p, double *chi, double *P);
void build_phase_vecs_scalar(int i_four, int n_coef, int n_mus, double **Y1, double *Y2, int lda, double *chi, double *P_pp, double *P_mp);
void build_phase_mats_scalar(int i_four, int n_coef, int n_mus1, int n_mus2, double **Y1, double **Y2, int lda, double *chi, double **P_pp, double **P_mp);
void build_scat_vector_gc(int n_coef, int n_stokes, double **gsf, double **chi, double *F);
void build_scat_matrix_gc(int n_coef, double **gsf, double **chi, double **F, int flag);
void build_scat_vector_lc(int n_coef, int n_stokes, double *p, double **chi, double *F);
void build_scat_matrix_lc(int n_coef, double *p, double **chi, double **F, int flag);
void basic_matrix(int i_four, int n_coef, int n_stokes, double mu, double ***A);
void build_phase_vecs_vector_gc(int i_four, int n_coef, int n_mus, int n_stokes, double *mu, double mu_0, double **gc, double *P_pp, double *P_mp, work_data work);
void build_phase_mats_vector_gc(int i_four, int n_coef, int n_mus1, int n_mus2, int n_stokes, double *mu1, double *mu2, double **gc, double **P_pp, double **P_mp, double **P_mm, double **P_pm, int vector, work_data work);
void scat_vector_rotate(int n_stokes, double mu1, double mu2, double MU, double d_phi, double *P1, double *P2);
void scat_vector_rotate_a(int n_stokes, double mu1, double mu2, double MU, double d_phi, double *P1_a, double *P2_a);
void phase_matrix_symmetry2(int n_quad, int n_stokes, double **P_pp, double **P_mp, double **P_mm, double **P_pm, double f);
void phase_matrix_symmetry_ldx2(int n_quad, int n_stokes, double *P_pp, double *P_mp, double *P_mm, double *P_pm, int ldx, double f);
void phase_matrix_symmetry3(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2, double **P_pp, double **P_mp, double **P_mm, double **P_pm, double f);
void phase_matrix_symmetry_ldx3(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2, double *P_pp, double *P_mp, double *P_mm, double *P_pm, int ldx, double f);
