/* xrtm_support.c */
int kernel_code(char *name);
const char *kernel_name(enum xrtm_kernel_type code);
int eigen_solver_gen_real_code(char *name);
const char *eigen_solver_gen_real_name(enum eigen_solver_gen_real_type code);
int eigen_solver_gen_complex_code(char *name);
const char *eigen_solver_gen_complex_name(enum eigen_solver_gen_complex_type code);
void *jalloc_array1(jmp_buf env, long m, int size);
void *jalloc_array2(jmp_buf env, long m, long n, int size);
void *jalloc_array3(jmp_buf env, long m, long n, long o, int size);
void *jalloc_array4(jmp_buf env, long m, long n, long o, long p, int size);
void *jalloc_array5(jmp_buf env, long m, long n, long o, long p, long q, int size);
int xrtm_set_layer_x(xrtm_data *d, double *x, double *x_, int i_layer, int n_layers, uchar *set_flags);
int xrtm_set_layer_x_l(xrtm_data *d, double **x, void *x_, int i_layer, int n_layers, int i_deriv, int n_derivs, int type);
