typedef struct
{
    double  c1_count;
    double  c2_count;
}   count_pair_t;

unsigned long   extreme_fcs_count(count_pair_t count_pairs[], unsigned long pair_count,
		      unsigned long replicates, double observed_fc,
		      unsigned long *fc_count);
void    fc_exact_p_val(count_pair_t count_pairs[], size_t pair_count,
			    size_t replicates, double observed_fc);
double  near_exact_p_val(double rep_counts1[], double rep_counts2[],
			   size_t replicates);

