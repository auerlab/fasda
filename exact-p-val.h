typedef struct
{
    double  c1_count;
    double  c2_count;
}   count_pair_t;

unsigned long   extreme_fcs_count(count_pair_t count_pairs[], unsigned long pair_count,
		      unsigned long replicates, double observed_fc,
		      unsigned long *fc_count);
double  fc_exact_pval(count_pair_t count_pairs[], size_t pair_count,
			    size_t replicates, double observed_fc);
double  near_exact_pval(double rep_counts1[], double rep_counts2[],
			   size_t replicates);
double  mann_whitney_pval(double rep_counts1[], double rep_counts2[],
			   size_t num_reps1, size_t num_reps2);

