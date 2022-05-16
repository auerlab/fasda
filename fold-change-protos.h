/* fold-change.c */
int main(int argc, char *argv[]);
int fold_change(FILE *condition_streams[], int conditions);
void print_header(int conditions);
void print_fold_change(const char *id, double condition_counts[], int conditions);
double dsv_total_counts(dsv_line_t *dsv_line);
void usage(char *argv[]);
