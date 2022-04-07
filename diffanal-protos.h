/* diffanal.c */
int main(int argc, char *argv[]);
int bl_sam_gff_cmp(bl_sam_t *alignment, bl_gff_t *feature);
int bl_gff_sam_cmp(bl_gff_t *feature, bl_sam_t *alignment);
int gff_find_overlapping_alignment(bl_gff_t *feature, FILE *stream, bl_sam_t *alignment);
double count_coverage(bl_gff_t *feature, bl_sam_t *alignment, FILE *sam_stream);
int64_t bl_gff_sam_overlap(bl_gff_t *feature, bl_sam_t *alignment);
void usage(char *argv[]);
