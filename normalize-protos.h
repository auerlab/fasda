/* normalize.c */
int main(int argc, char *argv[]);
int mrn(int argc, char *argv[], int arg, FILE *norm_all_stream);
void skip_headers(char *abundance_files[], FILE *abundance_streams[], dsv_line_t dsv_line[], size_t sample_count);
void check_all_eof(char *abundance_files[], FILE *abundance_streams[], size_t sample, size_t sample_count);
void usage(char *argv[]);
