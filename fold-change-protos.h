/* fold-change.c */
int main(int argc, char *argv[]);
int fold_change(FILE *condition_streams[], int conditions);
void print_header(int conditions);
void print_fold_change(const char *id, double condition_counts[], int conditions, double *rep_counts[], size_t num_reps[]);
double mann_whitney_p_val(double rep_counts1[], double rep_counts2[], size_t num_reps1, size_t num_reps2);
double dsv_total_counts(dsv_line_t *dsv_line, double rep_counts[]);
void usage(char *argv[]);
