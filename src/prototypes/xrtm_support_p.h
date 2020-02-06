/* xrtm_support.c */
int eigen_solver_gen_real_n(void);
int eigen_solver_gen_real_max_name_length(void);
const char *eigen_solver_gen_real_index_to_name(int index);
int eigen_solver_gen_real_name_to_index(const char *name);
enum eigen_solver_gen_real_type eigen_solver_gen_real_index_to_value(int index);
int eigen_solver_gen_real_value_to_index(enum eigen_solver_gen_real_type type);
enum eigen_solver_gen_real_type eigen_solver_gen_real_name_to_value(const char *name);
const char *eigen_solver_gen_real_value_to_name(enum eigen_solver_gen_real_type type);
int eigen_solver_gen_complex_n(void);
int eigen_solver_gen_complex_max_name_length(void);
const char *eigen_solver_gen_complex_index_to_name(int index);
int eigen_solver_gen_complex_name_to_index(const char *name);
enum eigen_solver_gen_complex_type eigen_solver_gen_complex_index_to_value(int index);
int eigen_solver_gen_complex_value_to_index(enum eigen_solver_gen_complex_type type);
enum eigen_solver_gen_complex_type eigen_solver_gen_complex_name_to_value(const char *name);
const char *eigen_solver_gen_complex_value_to_name(enum eigen_solver_gen_complex_type type);
void *jalloc_array1(jmp_buf env, long m, int size);
void *jalloc_array2(jmp_buf env, long m, long n, int size);
void *jalloc_array3(jmp_buf env, long m, long n, long o, int size);
void *jalloc_array4(jmp_buf env, long m, long n, long o, long p, int size);
void *jalloc_array5(jmp_buf env, long m, long n, long o, long p, long q, int size);
int xrtm_set_layer_x(xrtm_data *d, double *x, double *x_, int i_layer, int n_layers, uchar *set_flags);
int xrtm_set_layer_x_l(xrtm_data *d, double **x, void *x_, int i_layer, int n_layers, int i_deriv, int n_derivs, int type);