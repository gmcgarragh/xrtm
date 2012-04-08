/* xrtm_scatter_a.c */
void phase_func_a(int n_coef, double *p, double *chi_a, double P_a);
void build_phase_vecs_scalar_a(int i_four, int n_coef, int n_mus, double **Y1, double *Y2, int lda, double *chi_a, double *P_pp_a, double *P_mp_a);
void build_phase_mats_scalar_a(int i_four, int n_coef, int n_mus1, int n_mus2, double **Y1, double **Y2, int lda, double *chi_a, double **P_pp_a, double **P_mp_a);
void build_scat_vector_gc_a(int n_coef, int n_stokes, double **gsf, double **chi_a, double *F_a);
void build_phase_vecs_vector_gc_a(int i_four, int n_coef, int n_mus, int n_stokes, double *mu, double mu_0, double **gc_a, double *P_pp_a, double *P_mp_a, work_data work);
void build_phase_mats_vector_gc_a(int i_four, int n_coef, int n_mus1, int n_mus2, int n_stokes, double *mu1, double *mu2, double **gc_a, double **P_pp_a, double **P_mp_a, double **P_mm_a, double **P_pm_a, int vector, work_data work);
void phase_matrix_symmetry_a3(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2, double **P_pp_a, double **P_mp_a, double **P_mm_a, double **P_pm_a, double f);
void phase_matrix_symmetry_a_ldx3(int n_quad1, int n_stokes1, int n_quad2, int n_stokes2, double *P_pp_a, double *P_mp_a, double *P_mm_a, double *P_pm_a, int ldx, double f);
