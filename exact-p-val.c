#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <xtend/math.h>
#include <xtend/array.h>
#include <xtend/mem.h>
#include "exact-p-val.h"

const extern int    Debug;

/*
 *  Find all possible means of fold-change triplets
 *  Count the number of means >= observed mean
 */

double  fc_exact_p_val(count_pair_t count_pairs[], size_t pair_count,
			    size_t replicates, double observed_fc)

{
    unsigned long   fc_count, extreme_fcs, actual_fc_count, c;
    double          p_val, p_val_sum;

    if ( replicates <= 10 )
	fc_count = xt_n_choose_k(pair_count, replicates);
    else
	fc_count = 0;

    if ( Debug )
    {
	printf("\n%zu choose %zu = %lu possible FCs from %lu pairs.\n",
	    pair_count, replicates, fc_count, replicates);
	printf("P-value = likelihood of FC from %lu pairs >= %0.3f or <= %0.3f\n\n",
		replicates, observed_fc, 1.0 / observed_fc);
    }
    
    // Run several reps with the same fold-changes for down-sampled FCs
    // to check stability
    p_val_sum = 0.0;
    for (c = 0; c < (replicates >= 5 ? 1 : 1); ++c)
    {
	// printf("%lu %lu\n", c, replicates);
	extreme_fcs = extreme_fcs_count(count_pairs, pair_count, replicates,
			    observed_fc, &actual_fc_count);
	
	p_val = (double)extreme_fcs / actual_fc_count;
	// printf("\nLower FC, higher stddev, and outlier counts cause higher P-values.\n");
	if ( Debug )
	    printf("FCs sampled = %lu  P-value = %lu / %lu = %0.3f\n\n",
		    actual_fc_count, extreme_fcs, actual_fc_count, p_val);
	p_val_sum += p_val;
    }
    
    return p_val_sum / c;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *      Compute exact or near-exact P-value for paired lists
 *      counts1 and counts2.  Must accept lists in the same format as
 *      mann_whitney_p-val().
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  Examples:
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-10-19  Jason Bacon Begin
 ***************************************************************************/

double  near_exact_p_val(double counts1[], double counts2[],
			   size_t replicates)

{
    unsigned long   c, c1, c2, pair_count, samples, extreme_fcs,
		    half_pair_count;

    double  c1_sum, c2_sum, observed_fc, *counts;
    count_pair_t    *count_pairs;
    if ( replicates > 10 )
    {
	fprintf(stderr, "near_exact_p_val(): Use mann_whitney_p_val() for reps > 10.\n");
	return 1.0;
    }

    c1_sum = c2_sum = 0;
    for (c = 0; c < replicates; ++c)
    {
	c1_sum += counts1[c];
	c2_sum += counts2[c];
    }
    observed_fc = c2_sum / c1_sum;
    if ( observed_fc < 1.0 )
	observed_fc = 1.0 / observed_fc;
    if ( Debug )
	printf("Observed: FC = %0.3f  1 / Observed FC = %0.3f\n",
		observed_fc, 1.0 / observed_fc);
	
    /*
     *  Compute fold-change for every possible pairing.
     *  #samples choose 2 * 2 (FC and 1/FC), since this affects
     *  the distribution of means of N samples
     *  Satisfies the null hypothesis P(n1 > n2) = P(n1 < N2)
     */
    
    samples = replicates * 2;
    if ( Debug )
    {
	printf("\n%lu choose %d = %lu combinations of counts\n",
		samples, 2, xt_n_choose_k(samples, 2));
	puts("2 ordered pairs for each combination:");
    }
    half_pair_count = xt_n_choose_k(samples, 2);
    pair_count = half_pair_count * 2;   // FC and 1/FC
    count_pairs = xt_malloc(pair_count, sizeof(*count_pairs));
    
    /*
     *  Include both FC and 1/FC for each count combination in the set
     *  so that the distribution of FCs covers both cases.
     */
    
    counts = xt_malloc(replicates * 2, sizeof(*counts));
    for (c = 0; c < replicates; ++c)
    {
	counts[c] = counts1[c];
	counts[c + replicates] = counts2[c];
    }

    c = extreme_fcs = 0;
    for (c1 = 0; c1 < replicates * 2; ++c1)
    {
	for (c2 = c1 + 1; c2 < replicates * 2; ++c2)
	{
	    count_pairs[c].c1_count = counts[c1];
	    count_pairs[c].c2_count = counts[c2];
	    count_pairs[c + half_pair_count].c1_count = counts[c2];
	    count_pairs[c + half_pair_count].c2_count = counts[c1];
	    ++c;
	}
    }

    // Check for program bugs
    if ( c != half_pair_count )
    {
	fprintf(stderr, "%lu != %lu\n", c, half_pair_count);
	exit(EX_SOFTWARE);
    }

    // Down-sampled P-values come up light without shuffling, but average
    // very close to exact with this shuffling.  Why?
    xt_shuffle(count_pairs, pair_count, sizeof(*count_pairs));
    
    if ( Debug )
    {
	for (c = 0; c < pair_count; ++c)
	    printf("%2lu %3.0f, %3.0f   FC = %0.3f\n", c,
		    count_pairs[c].c1_count, count_pairs[c].c2_count,
		    (double)count_pairs[c].c1_count / count_pairs[c].c2_count);
    }
    
    return fc_exact_p_val(count_pairs, pair_count, replicates, observed_fc);
}
