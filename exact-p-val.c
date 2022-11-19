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

double  fc_exact_pval(count_pair_t count_pairs[], size_t pair_count,
			    size_t replicates, double observed_fc)

{
    unsigned long   fc_count, extreme_fcs, actual_fc_count, c, passes;
    double          pval, pval_sum, pval_low, pval_high;

    // fc_count is beyond the range of an unsigned long for > 10 replicates
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
	passes = 5;
    }
    else
	passes = 1;
    
    // Run several reps with the same fold-changes for down-sampled FCs
    // to check stability
    pval_sum = 0.0;
    pval_low = 1.0;
    pval_high = 0.0;
    for (c = 0; c < (replicates >= 5 ? passes : 1); ++c)
    {
	// fprintf(stderr, "c = %lu replicates = %lu\n", c, replicates);
	extreme_fcs = extreme_fcs_count(count_pairs, pair_count, replicates,
			    observed_fc, &actual_fc_count);
	
	pval = (double)extreme_fcs / actual_fc_count;
	// printf("\nLower FC, higher stddev, and outlier counts cause higher P-values.\n");
	if ( Debug )
	    printf("Pass %lu: FCs sampled = %lu  P-value = %lu / %lu = %0.3f\n\n",
		    c, actual_fc_count, extreme_fcs, actual_fc_count, pval);
	pval_sum += pval;
	if ( pval < pval_low ) pval_low = pval;
	if ( pval > pval_high ) pval_high = pval;
    }
    
    if ( Debug )
	printf("P-value span = %0.3f\n\n", pval_high - pval_low);
    return pval_sum / c;
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

double  near_exact_pval(double counts1[], double counts2[], size_t replicates)

{
    unsigned long   c, c1, c2, pair_count, samples, extreme_fcs,
		    half_pair_count;
    const unsigned long max_reps = 12;
    double  c1_sum, c2_sum, observed_fc, *counts;
    count_pair_t    *count_pairs;

    if ( replicates > max_reps )
    {
	fprintf(stderr, "near_exact_pval(): Use mann_whitney_pval() for reps > %lu.\n",
		max_reps);
	return 1.0;
    }

    // Replace 0 counts with 1/1000 the mean of the other condition
    // to avoid FCs of 0 or infinity
    adjust_counts(counts1, counts2, replicates);
    
    c1_sum = c2_sum = 0;
    for (c = 0; c < replicates; ++c)
    {
	c1_sum += counts1[c];
	c2_sum += counts2[c];
    }
    
    // Don't waste time doing combinatorics when both counts are 0
    if ( (c1_sum == 0.0) && (c2_sum == 0.0) )
	return 1.0;
    
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
    
    return fc_exact_pval(count_pairs, pair_count, replicates, observed_fc);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *      Counts of produce fold changes of 0 or infinity, causing problems
 *      for exact P-value computation.  Adjust counts of 0 to at least
 *      1/1000 the mean of the other condition.
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
 *  2022-11-07  Jason Bacon Begin
 ***************************************************************************/

void    adjust_counts(double counts1[], double counts2[], unsigned long replicates)

{
    unsigned long   c;
    double          c1_sum, c2_sum, c1_mean, c2_mean;
    
    for (c1_sum = c2_sum = 0.0, c = 0; c < replicates; ++c)
    {
	c1_sum += counts1[c];
	c2_sum += counts2[c];
    }
    
    if ( (c1_sum == 0.0) && (c2_sum == 0.0) )
    {
	for (c = 0; c < replicates; ++c)
	    counts1[c] = counts2[c] = 1.0;
    }
    else
    {
	c1_mean = c1_sum / replicates;
	c2_mean = c2_sum / replicates;
	for (c = 0; c < replicates; ++c)
	{
	    if ( counts1[c] == 0.0 )
		counts1[c] = c2_mean / 1000.0;
	    if ( counts2[c] == 0.0 )
		counts2[c] = c1_mean / 1000.0;
	}
    }
    
    if ( Debug )
    {
	puts("Adjusted counts:");
	for (c = 0; c < replicates; ++c)
	    printf("%f %f\n", counts1[c], counts2[c]);
    }
}

