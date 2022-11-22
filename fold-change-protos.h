/* fold-change.c */
int main(int argc, char *argv[]);
int fold_change(FILE *condition_streams[], int conditions, FILE *diff_stream, unsigned flags);
void print_header(FILE *diff_stream, int conditions);
void print_fold_change(FILE *diff_stream, const char *id, double cond_tot_counts[], double condition_stddevs[], int conditions, double *rep_counts[], size_t num_repls[], unsigned flags);
unsigned agreement(int c1, int c2, size_t num_repls[], double *rep_counts[]);
double dsv_total_counts(dsv_line_t *dsv_line, double rep_counts[], double *condition_stddevs);
void usage(char *argv[]);
