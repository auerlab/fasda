/* fold-change.c */
int main(int argc, char *argv[]);
int fold_change(FILE *condition_streams[], int conditions, FILE *diff_stream, unsigned flags);
void print_header(FILE *diff_stream, int conditions);
void print_fold_change(FILE *diff_stream, const char *id, double condition_counts[], double condition_stddevs[], int conditions, double *rep_counts[], size_t num_reps[], unsigned flags);
double dsv_total_counts(dsv_line_t *dsv_line, double rep_counts[], double *condition_stddevs);
void usage(char *argv[]);
