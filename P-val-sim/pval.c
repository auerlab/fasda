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
#include <time.h>
#include <limits.h>
#include <unistd.h>     // getpid()
#include "pval.h"

void    usage(char *argv[])

{
    fprintf(stderr, "Usage:   %s count1-mean count2-mean max-deviation\n", argv[0]);
    fprintf(stderr, "Example: %s 100 200 .2 3\n",argv[0]);
    exit(EX_USAGE);
}

int     main(int argc,char *argv[])

{
    unsigned long   c, c1, c2, replicates,
		    less, more, equal, *counts,
		    count1_mean, count2_mean, temp, pair_count, samples, fc_ge,
		    half_pair_count;
    double  c1_sum, c2_sum, max_deviation, observed_fc;
    count_pair_t    *count_pairs;

    if ( argc != 5 )
	usage(argv);

    count1_mean = atoi(argv[1]);
    count2_mean = atoi(argv[2]);
    max_deviation = atof(argv[3]);
    replicates = atoi(argv[4]);
    samples = replicates * 2;
    
    counts = malloc(samples * sizeof(*counts));
    
    // Some code below assumes count2_mean > count1_mean
    if ( count1_mean > count2_mean )
    {
	temp = count1_mean;
	count1_mean = count2_mean;
	count2_mean = temp;
    }
    
    /*
     *  Generate samples random counts with FC around count2 / count1
     *  These are the "observed" counts
     */
    printf("\ncount1 = %lu +/- to up to %0.0f%%, count2 = %lu +/- up to %0.0f\n",
	    count1_mean, max_deviation * count1_mean,
	    count2_mean, max_deviation * count2_mean);
    puts("Cond1 Cond2 FC      1/FC");

    // mean FC = sigma(c2) / sigma(c1), not averages of FCs for each condition
    // P. Auer
    // Comment this out to get the same counts repeatedly
    srandom(getpid());
    c1_sum = c2_sum = 0.0;
    for (c = 0; c < replicates; ++c)
    {
	counts[c] = count1_mean
		 + random() % (int)(count1_mean * max_deviation * 2)
		 - count1_mean * max_deviation;
	counts[c + replicates] = count2_mean
		 + random() % (int)(count2_mean * max_deviation * 2)
		 - count2_mean * max_deviation;
	printf("%5lu %5lu\n", counts[c], counts[c + replicates]);
	c1_sum += counts[c];
	c2_sum += counts[c + replicates];
    }
    
    /*
     *  Compute mean and stddev for "observed" values
     */
    
    observed_fc = c2_sum / c1_sum;
    
    /*
     *  Compute fold-change for every possible pairing.
     *  #samples choose 2 * 2 (FC and 1/FC), since this affects
     *  the distribution of means of N samples
     *  Satisfies the null hypothesis P(n1 > n2) = P(n1 < N2)
     */
    
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
    
    c = fc_ge = less = more = equal = 0;
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
    
    for (c = 0; c < pair_count; ++c)
	printf("%2lu %3lu, %3lu   FC = %0.5f\n", c,
		count_pairs[c].c1_count, count_pairs[c].c2_count,
		(double)count_pairs[c].c1_count / count_pairs[c].c2_count);
    
    fc_mean_exact_p_val(count_pairs, pair_count, replicates, observed_fc);
    return EX_OK;
}


/*
 *  Find all possible means of fold-change triplets
 *  Count the number of means >= observed mean
 */

void    fc_mean_exact_p_val(count_pair_t count_pairs[], size_t pair_count,
			    size_t replicates, double observed_fc)

{
    unsigned long   fc_count, fc_ge, c, actual_fc_count;

    if ( replicates <= 10 )
    {
	fc_count = xt_n_choose_k(pair_count, replicates);
	printf("\n%zu choose %zu = %lu possible FCs from %lu pairs.\n",
	    pair_count, replicates, fc_count, replicates);
    }
    else
	printf("\nfc_mean_count > 2^64 for replicates > 10.\n");

    printf("Observed FC = %0.5f  1 / Observed FC = %0.5f\n",
	    observed_fc, 1.0 / observed_fc);
    printf("\nP-value: Likelihood of FC from %lu pairs >= %0.5f or <= %0.5f\n",
	    replicates, observed_fc, 1.0 / observed_fc);

    // Run 10 reps for down-sampled FCs to check stability
    for (c = 0; c < (replicates > 5 ? 5 : 1); ++c)
    {
	fc_ge = fc_ge_count(count_pairs, pair_count, replicates,
			    observed_fc, &actual_fc_count);
	printf("\nLower FC, higher stddev, and outlier counts cause higher P-values.\n");
	printf("Observed: fc = %0.5f  Observed 1/fc = %0.5f\n",
		observed_fc, 1.0 / observed_fc);
	printf("FC count = %lu  FC >= %0.5f = %lu  P(FC >= %0.5f) = %0.5f\n\n",
		fc_count, observed_fc, fc_ge, observed_fc,
		(double)fc_ge / fc_count);
    }
}
