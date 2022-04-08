/* diffanal.c */
int main(int argc, char *argv[]);
int bl_gff_find_overlapping_alignment(bl_gff_t *feature, FILE *stream, char *last_chrom, bl_sam_t *alignment);
double count_coverage(bl_gff_t *feature, bl_sam_t *alignment, FILE *sam_stream, char *last_chrom);
void usage(char *argv[]);
