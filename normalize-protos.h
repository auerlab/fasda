/* normalize.c */
int main(int argc, const char *argv[]);
int mrn(const char *abundance_files[], FILE *norm_all_stream);
void close_all_streams(FILE *abundance_streams[], FILE *norm_sample_streams[], FILE *norm_all_stream, size_t sample_count);
void skip_headers(const char *abundance_files[], FILE *abundance_streams[], dsv_line_t dsv_line[], size_t sample_count);
void check_all_eof(const char *abundance_files[], FILE *abundance_streams[], size_t sample, size_t sample_count);
void usage(const char *argv[]);
