/* xrtm_save_tree.c */
int save_tree_init(save_tree_data *d);
void save_tree_free(save_tree_data *d);
void save_tree_encode_s(save_tree_data *d, const char *s);
void save_tree_encode_i(save_tree_data *d, int i);
void save_tree_encode_i_j(save_tree_data *d, int i, int j);
void save_tree_decode_s(save_tree_data *d, const char *s);
void save_tree_decode_i(save_tree_data *d);
void save_tree_decode_i_j(save_tree_data *d);
void save_tree_recode_s(save_tree_data *d, const char *s);
void save_tree_recode_i(save_tree_data *d, int i, int flag);
void save_tree_recode_i_j(save_tree_data *d, int i, int j, int flag);
int __save_tree_retrieve_data(save_tree_data *d, size_t size, void **v);
int __save_tree_retrieve_proxy(save_tree_data *d, const char *s, size_t size, void **v);
