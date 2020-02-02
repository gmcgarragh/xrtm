/* test_xrtm.c */
int xrtm_input_init(xrtm_input_data *d);
void xrtm_input_free(xrtm_input_data *d);
int test_xrtm_init(test_xrtm_data *d);
void test_xrtm_free(test_xrtm_data *d);
int test_xrtm_fill(test_xrtm_data *d, int options, int solvers, int max_coef, int n_quad, int n_stokes, int n_derivs, int n_layers, int n_kernels, int n_kernel_quad, enum xrtm_kernel_type *kernels, int n_out_levels, int n_out_thetas);
