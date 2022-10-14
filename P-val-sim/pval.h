typedef struct
{
    unsigned long   c1_count;
    unsigned long   c2_count;
}   count_pair_t;

unsigned long   fc_ge_count(count_pair_t count_pairs[], unsigned long pair_count,
		      unsigned long replicates, double observed_fc,
		      unsigned long *fc_count);
void    fc_mean_exact_p_val(count_pair_t count_pairs[], size_t pair_count,
			    size_t replicates, double observed_fc);
