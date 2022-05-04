/* fold-change.c */
int main(int argc, char *argv[]);
int fold_change(FILE *condition_streams[], int conditions);
void print_header(int conditions);
void print_fold_change(const char *id, double coverage[], int conditions);
void usage(char *argv[]);
