/* exact-p-val.c */
double fc_exact_pval(count_pair_t count_pairs[], size_t pair_count, size_t replicates, double observed_fc);
double near_exact_pval(double counts1[], double counts2[], size_t replicates);
void adjust_counts(double counts1[], double counts2[], unsigned long replicates);
