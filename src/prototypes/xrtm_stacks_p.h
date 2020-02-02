/* xrtm_stacks.c */
void print_stacks(stack_data *stacks, int n_stacks);
int build_stack_chain(int n_layers, int n_derivs, uchar **derivs, stack_data *stacks);
int stack_chain_alloc(int n_four, int n_quad, int n_derivs, int n_layers, int n_stacks0, stack_data *stacks0, int n_stacks, stack_data *stacks, uchar **derivs_layers, uchar **derivs_s, uchar **derivs_d);
int stack_chain_free(int n_four, int n_quad, int n_derivs, int n_layers, int n_stacks, stack_data *stacks, uchar **derivs_layers, uchar **derivs_s, uchar **derivs_d);
int stack_chain_find(int n_stacks, stack_data *stacks, int i1, int i2);
int stack_chain_find2(int n_stacks, stack_data *stacks, int i11, int i12, int i21, int i22);
