/* abundance.c */
int main(int argc, char *argv[]);
int abundance(FILE *feature_stream, FILE *sam_streams[], FILE *abundance_streams[], int file_count, const char *feature_type, int flags);
int bl_gff_find_overlapping_alignment(bl_gff_t *feature, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_sam_t *alignment, bl_alignment_stats_t *alignment_stats, int flags);
double count_coverage(bl_gff_t *feature, bl_sam_t *alignment, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_alignment_stats_t *alignment_stats, int flags);
int print_abundance(FILE *abundance_stream, bl_gff_t *feature, int64_t length, double est_counts, int flags);
int fasda_chrom_name_cmp(const char *s1, const char *s2, int flags);
void usage(char *argv[]);
