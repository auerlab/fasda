/* abundance.c */
int main(int argc, char *argv[]);
int exact_abundance(const char *feature_file, char *sam_files[], FILE *abundance_streams[], int file_count, const char *feature_type, int flags);
int stringtie_abundance(const char *feature_file, char *sam_files[], char *gtf_files[], FILE *abundance_streams[], int file_count, const char *feature_type, unsigned read_length, int flags);
int bl_gff_find_overlapping_alignment(bl_gff_t *feature, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_sam_t *alignment, bl_alignment_stats_t *alignment_stats, int flags);
double count_coverage(bl_gff_t *feature, bl_sam_t *alignment, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_alignment_stats_t *alignment_stats, int flags);
int print_abundance(FILE *abundance_stream, bl_gff_t *feature, int64_t length, double est_counts, int flags);
int fasda_chrom_name_cmp(const char *s1, const char *s2, int flags);
char *bl_sam_filename_extension(const char *filename);
void sort_abundance(char *abundance_files[], int file_count);
void usage(char *argv[]);
