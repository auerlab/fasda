/* abundance.c */
int main(int argc, char *argv[]);
int abundance(FILE *feature_stream, FILE *sam_streams[], FILE *abundance_streams[], int file_count, int flags);
int bl_gff_find_overlapping_alignment(bl_gff_t *feature, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_sam_t *alignment, bl_alignment_stats_t *alignment_stats);
double count_coverage(bl_gff_t *feature, bl_sam_t *alignment, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_alignment_stats_t *alignment_stats);
int print_abundance(FILE *abundance_stream, bl_gff_t *feature, double counts, int flags);
void usage(char *argv[]);
/* alignment-stats-mutators.c */
int bl_alignment_stats_set_total(bl_alignment_stats_t *bl_alignment_stats_ptr, unsigned long new_total);
int bl_alignment_stats_set_overlapping(bl_alignment_stats_t *bl_alignment_stats_ptr, unsigned long new_overlapping);
/* diffanal.c */
int main(int argc, char *argv[]);
int diffanal(FILE *feature_stream, FILE *sam_streams[], int conditions, int flags);
int bl_gff_find_overlapping_alignment(bl_gff_t *feature, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_sam_t *alignment, bl_alignment_stats_t *alignment_stats);
double count_coverage(bl_gff_t *feature, bl_sam_t *alignment, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_alignment_stats_t *alignment_stats);
void print_header(int conditions);
void print_fold_change(bl_gff_t *feature, double coverage[], int conditions, int flags);
void usage(char *argv[]);
/* fold-change.c */
int main(int argc, char *argv[]);
int diffanal(FILE *feature_stream, FILE *sam_streams[], int conditions, int flags);
int bl_gff_find_overlapping_alignment(bl_gff_t *feature, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_sam_t *alignment, bl_alignment_stats_t *alignment_stats);
double count_coverage(bl_gff_t *feature, bl_sam_t *alignment, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_alignment_stats_t *alignment_stats);
void print_header(int conditions);
void print_fold_change(bl_gff_t *feature, double coverage[], int conditions, int flags);
void usage(char *argv[]);
/* normalize.c */
int main(int argc, char *argv[]);
void usage(char *argv[]);
