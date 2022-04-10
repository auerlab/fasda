/* diffanal.c */
int main(int argc, char *argv[]);
int diffanal(FILE *feature_stream, FILE *condition_streams[], int conditions);
int bl_gff_find_overlapping_alignment(bl_gff_t *feature, FILE *alignment_stream, FILE *buffer_stream, char *previous_alignment_chrom, bl_sam_t *alignment);
double count_coverage(bl_gff_t *feature, bl_sam_t *alignment, FILE *sam_stream, FILE *buffer_stream, char *previous_alignment_chrom);
void print_header(int conditions);
void print_fold_change(bl_gff_t *feature, double coverage[], int conditions);
void usage(char *argv[]);
