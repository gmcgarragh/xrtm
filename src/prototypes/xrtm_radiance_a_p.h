/* xrtm_radiance_a.c */
void radiance_slab_a(int n_quad, double **R_m, double **T_p, double *S_p, double **R_m_a, double **T_p_a, double *S_p_a, double *I1_m, double *I2_p, double *I1_m_a, double *I2_p_a, double *I1_p, double *I1_p_a, work_data work);
void radiance_toa_ref_a(int n_quad, double **R_m, double *S_p, double **R_m_a, double *S_p_a, double *I1_m, double *I1_m_a, double *I1_p, double *I1_p_a, work_data work);
void radiance_toa_all_a(int n_quad, double **R_m, double **T_p, double *S_p, double **R_m_a, double **T_p_a, double *S_p_a, double *I1_m, double *I2_p, double *I1_m_a, double *I2_p_a, double *I1_p, double *I1m_, double *I1_p_a, double *I1m_a_, work_data work);
int forward_save_radiance_boa_all_alloc(forward_save_radiance_boa_all_data *d, int n_quad);
void radiance_boa_all_a(int n_quad, double **R12_p, double **T12_m, double **R23_m, double *S12_m, double **R12_p_a, double **T12_m_a, double **R23_m_a, double *S12_m_a, double *I3_p, double *I1_m, double *I3_p_a, double *I1_m_a, double *I2_p, double *I2_m, double *I2_p_a, double *I2_m_a, save_tree_data save_tree, work_data work);
void radiance_slab_tl_with_ad(int n_quad, int n_derivs, double **R_m, double **T_p, double *S_p, double ***R_m_l, double ***T_p_l, double **S_p_l, double *I1_m, double *I2_p, double **I1_m_l, double **I2_p_l, double *I1_p, double **I1_p_l, uchar *derivs, work_data work);
