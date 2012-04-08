/* xrtm_doubling.c */
void layer_double(double **R, double **T, double *Se_m, double *Se_p, double *Sl_m, double *Sl_p, int n, double atran, double lin_fac, int flag1, int flag2, int flag3, work_data work);
void layer_double_s(double **R, double **T, double *S_m, double *S_p, int n, double atran, work_data work);
void layer_double_l(double **R, double **T, double *S_m, double *S_p, double ***R_l, double ***T_l, double **S_m_l, double **S_p_l, int n, int n_derivs, double atran, double *atran_l, uchar *derivs_h, uchar *derivs_p, work_data work);
void layer_double_s_l(double **R, double **T, double *S_m, double *S_p, double ***R_l, double ***T_l, double **S_m_l, double **S_p_l, int n, int n_derivs, double atran, double *atran_l, uchar *derivs_h, uchar *derivs_p, work_data work);
