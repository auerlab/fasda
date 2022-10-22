/***************************************************************************
 *  Description:
 *      Test driver for computing exact P-values of differential
 *      expression data.
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-09-24  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <xtend/math.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>     // getpid()
#include <sys/time.h>
#include "pval.h"

void    usage(char *argv[])

{
    fprintf(stderr, "Usage:   %s count1-mean count2-mean max-deviation replicates iterations\n", argv[0]);
    fprintf(stderr, "Example: %s 100 200 .2 3 5\n",argv[0]);
    exit(EX_USAGE);
}


int     main(int argc,char *argv[])

{
    unsigned long   i, c, replicates, iterations,
		    count1_mean, count2_mean;
    double          max_deviation, *counts1, *counts2;
    struct timeval  time;

    if ( argc != 6 )
	usage(argv);

    count1_mean = atoi(argv[1]);
    count2_mean = atoi(argv[2]);
    max_deviation = atof(argv[3]);
    replicates = atoi(argv[4]);
    iterations = atoi(argv[5]);
    
    counts1 = malloc(replicates * sizeof(*counts1));
    counts2 = malloc(replicates * sizeof(*counts2));
    
    /*
     *  Generate samples random counts with FC around count2 / count1
     *  These are the "observed" counts
     */
    printf("\ncount1 = %lu +/- to up to %0.0f%%, count2 = %lu +/- up to %0.0f%%\n",
	    count1_mean, max_deviation * count1_mean,
	    count2_mean, max_deviation * count2_mean);
    
    for (i = 0; i < iterations; ++i)
    {
	gettimeofday(&time, NULL);
	srandom(time.tv_usec);
	//if ( i % 100 == 0 )
	//    fprintf(stderr, "%lu\r", i);
	
	puts("Cond1 Cond2 FC      1/FC");
    
	// mean FC = sigma(c2) / sigma(c1), not averages of FCs for each condition
	// P. Auer
	// Comment this out to get the same counts repeatedly
	for (c = 0; c < replicates; ++c)
	{
	    counts1[c] = count1_mean
		     + random() % (unsigned long)(count1_mean * max_deviation * 2)
		     - count1_mean * max_deviation;
	    counts2[c] = count2_mean
		     + random() % (unsigned long)(count2_mean * max_deviation * 2)
		     - count2_mean * max_deviation;
	    printf("%5.0f %5.0f\n", counts1[c], counts2[c]);
	}
	
	near_exact_p_val(counts1, counts2, replicates);
	printf("Mann-Whitney = %0.3f\n",
		mann_whitney_p_val(counts1, counts2, replicates, replicates));
    }
    return EX_OK;
}


/*
 *  Find all possible means of fold-change triplets
 *  Count the number of means >= observed mean
 */

void    fc_exact_p_val(count_pair_t count_pairs[], size_t pair_count,
			    size_t replicates, double observed_fc)

{
    unsigned long   fc_count, extreme_fcs, actual_fc_count, c;

    if ( replicates <= 10 )
    {
	fc_count = xt_n_choose_k(pair_count, replicates);
	printf("\n%zu choose %zu = %lu possible FCs from %lu pairs.\n",
	    pair_count, replicates, fc_count, replicates);
    }
    else
	printf("\nfc_count > 2^64 for replicates > 10.\n");

    printf("P-value = likelihood of FC from %lu pairs >= %0.3f or <= %0.3f\n\n",
	    replicates, observed_fc, 1.0 / observed_fc);

    // Run several reps with the same fold-changes for down-sampled FCs
    // to check stability
    for (c = 0; c < (replicates >= 5 ? 5 : 1); ++c)
    {
	// printf("%lu %lu\n", c, replicates);
	extreme_fcs = extreme_fcs_count(count_pairs, pair_count, replicates,
			    observed_fc, &actual_fc_count);
	
	// Sanity check
	/* Does not apply when down sampling
	if ( actual_fc_count != fc_count )
	{
	    fprintf(stderr, "FC counf mismatch: %lu != %lu\n",
		    actual_fc_count, fc_count);
	    exit(EX_SOFTWARE);
	}
	*/
	
	// printf("\nLower FC, higher stddev, and outlier counts cause higher P-values.\n");
	printf("FCs sampled = %lu  P-value = %lu / %lu = %0.3f\n\n",
		actual_fc_count, extreme_fcs, actual_fc_count,
		(double)extreme_fcs / actual_fc_count);
    }
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
    printf("Observed: FC = %0.3f  1 / Observed FC = %0.3f\n",
	    observed_fc, 1.0 / observed_fc);
	
    /*
     *  Compute fold-change for every possible pairing.
     *  #samples choose 2 * 2 (FC and 1/FC), since this affects
     *  the distribution of means of N samples
     *  Satisfies the null hypothesis P(n1 > n2) = P(n1 < N2)
     */
    
    samples = replicates * 2;
    printf("\n%lu choose %d = %lu combinations of counts\n",
	    samples, 2, xt_n_choose_k(samples, 2));
    puts("2 ordered pairs for each combination:");
    half_pair_count = xt_n_choose_k(samples, 2);
    pair_count = half_pair_count * 2;   // FC and 1/FC
    count_pairs = malloc(pair_count * sizeof(*count_pairs));
    
    /*
     *  Include both FC and 1/FC for each count combination in the set
     *  so that the distribution of FCs covers both cases.
     */
    
    counts = malloc(replicates * 2 * sizeof(*counts));
    for (c = 0; c < replicates; ++c)
    {
	counts[c] = counts1[c];
	counts[c + replicates] = counts2[c];
    }
    for (c = 0; c < replicates * 2; ++c)
	printf("%f\n", counts[c]);

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
	printf("%lu != %lu\n", c, half_pair_count);
	return 1;
    }
    
    // Shuffle
    // Down-sampled P-values come up light without shuffling, but average
    // very close to exact with this shuffling.  Why?
    for (c = 0; c < pair_count - 1; ++c)
    {
	count_pair_t    temp;
	
	c1 = c + 1 + random() % (pair_count - c - 1);
	printf("Swapping %lu with %lu\n", c, c1);
	temp = count_pairs[c];
	count_pairs[c] = count_pairs[c1];
	count_pairs[c1] = temp;
    }
    
    for (c = 0; c < pair_count; ++c)
	printf("%2lu %3.0f, %3.0f   FC = %0.3f\n", c,
		count_pairs[c].c1_count, count_pairs[c].c2_count,
		(double)count_pairs[c].c1_count / count_pairs[c].c2_count);
    
    fc_exact_p_val(count_pairs, pair_count, replicates, observed_fc);
    
    return 1.0;
}
