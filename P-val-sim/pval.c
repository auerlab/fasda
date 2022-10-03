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

unsigned long   fc_ge_count(double fc_list[], unsigned long fc_count,
		      unsigned long replicates, double observed_fc_mean,
		      unsigned long *fc_mean_count);
void    fc_mean_exact_p_val(double fc_list[], size_t fc_count,
			    size_t replicates,
			    double observed_fc_mean, double observed_fc_stddev,
			    double dist_fc_mean);

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
		    count1_mean, count2_mean, fc_count, samples, fc_ge,
		    half_fc_count;
    double  fc, observed_fc_mean, *fc_list, max_deviation, fc_sum,
	    dist_fc_mean, fc_var_sum, observed_fc_stddev;

    if ( argc != 5 )
	usage(argv);

    count1_mean = atoi(argv[1]);
    count2_mean = atoi(argv[2]);
    max_deviation = atof(argv[3]);
    replicates = atoi(argv[4]);
    samples = replicates * 2;
    
    counts = malloc(samples * sizeof(*counts));
    
    /*
     *  Generate samples random counts with FC around count2 / count1
     *  These are the "observed" counts
     */
    printf("\ncount1 = %lu +/- to up to %0.0f%%, count2 = %lu +/- same\n",
	    count1_mean, max_deviation * 100, count2_mean);
    puts("Cond1 Cond2");
    srandom(time(NULL));
    for (c = 0, observed_fc_mean = 0.0; c < replicates; ++c)
    {
	counts[c] = count1_mean
		 + random() % (int)(count1_mean * max_deviation * 2)
		 - count1_mean * max_deviation;
	counts[c + replicates] = count2_mean
		 + random() % (int)(count2_mean * max_deviation * 2)
		 - count2_mean * max_deviation;
	fc = (double)counts[c + replicates] / counts[c];
	printf("%5lu %5lu %0.4f\n", counts[c], counts[c + replicates], fc);
	observed_fc_mean += fc;
    }
    
    puts("\nPooled counts:");
    for (c = 0; c < replicates * 2; ++c)
	printf("%3lu\n", counts[c]);
    
    /*
     *  Compute mean and stddev for "observed" values
     */
    
    observed_fc_mean /= replicates;
    fc_var_sum = 0;
    if ( observed_fc_mean < 1.0 )
    {
	observed_fc_mean = 1.0 / observed_fc_mean;
	for (c = 0; c < replicates; ++c)
	{
	    fc = counts[c] / (double)counts[c + replicates];
	    fc_var_sum += (fc - observed_fc_mean) * (fc - observed_fc_mean);
	}
    }
    else
    {
	for (c = 0; c < replicates; ++c)
	{
	    fc = (double)counts[c + replicates] / counts[c];
	    fc_var_sum += (fc - observed_fc_mean) * (fc - observed_fc_mean);
	}
    }
    observed_fc_stddev = sqrt(fc_var_sum / replicates);
    printf("FC mean = %0.4f  FC stddev = %f\n",
	    observed_fc_mean, observed_fc_stddev);
    
    /*
     *  Compute fold-change for every possible pairing.
     *  #samples choose 2 * 2 (FC and 1/FC), since this affects
     *  the distribution of means of N samples
     *  Satisfies the null hypothesis P(n1 > n2) = P(n1 < N2)
     */
    
    puts("\nFold-change of every possible pairing of samples:");
    puts("(Counts are not paired with themselves.)");
    half_fc_count = xt_n_choose_k(samples, 2);
    fc_count = half_fc_count * 2;   // FC and 1/FC
    fc_list = malloc(fc_count * sizeof(*fc_list));
    
    /*
     *  FIXME: Is this needed?  It greatly increases the size of the
     *  FC means set.
     *  Include both FC and 1/FC for each count combination in the set
     *  so that the distribution of FC means covers both cases.
     */
    
    c = fc_sum = fc_ge = less = more = equal = 0;
    for (c1 = 0; c1 < replicates * 2; ++c1)
    {
	for (c2 = c1 + 1; c2 < replicates * 2; ++c2)
	{
	    fc_list[c] = (double)counts[c1] / counts[c2];
	    fc_list[c + half_fc_count] = 1.0 / fc_list[c];
	    printf("%2lu %3lu / %3lu = %0.4f\n", c,
		    counts[c1], counts[c2], fc_list[c]);
	    printf("%2lu %3lu / %3lu = %0.4f\n", c + half_fc_count,
		    counts[c2], counts[c1], fc_list[c + half_fc_count]);
	    if ( fc_list[c] >= observed_fc_mean )
		++fc_ge;
	    
	    if ( fc_list[c] < 1.0 )
		++less, ++more;
	    else if ( fc_list[c] > 1.0 )
		++more, ++less;
	    else
		equal += 2;
	    fc_sum += fc_list[c] + fc_list[c + half_fc_count];
	    ++c;
	}
    }
    
    // Check for program bugs
    if ( c != half_fc_count )
    {
	printf("%lu != %lu\n", c, half_fc_count);
	return 1;
    }
    
    dist_fc_mean = fc_sum / fc_count;
    printf("\nless + more + equal should be %lu.  FC mean should be slightly > 1.\n",
	    fc_count);
    printf("less should equal more to satisfy the null hypothesis.\n");
    printf("Distribution: less = %lu  more = %lu  equal = %lu  FC mean = %0.4f\n",
	    less, more, equal, dist_fc_mean);

    fc_mean_exact_p_val(fc_list, fc_count, replicates,
			observed_fc_mean, observed_fc_stddev, dist_fc_mean);
    return EX_OK;
}


/*
 *  Find all possible means of fold-change triplets
 *  Count the number of means >= observed mean
 */

void    fc_mean_exact_p_val(double fc_list[], size_t fc_count,
			    size_t replicates,
			    double observed_fc_mean, double observed_fc_stddev,
			    double dist_fc_mean)

{
    unsigned long   fc_mean_count, fc_ge;

    // FIXME: This will overflow unsigned long for replicates > 10
    fc_mean_count = xt_n_choose_k(fc_count, replicates);
    printf("%zu choose %zu = %lu possible means of %lu FCs\n",
	    fc_count, replicates, fc_mean_count, replicates);
    
    fc_ge = fc_ge_count(fc_list, fc_count, replicates,
			observed_fc_mean, &fc_mean_count);
    printf("\nObserved FC mean = %0.4f  Observed FC stddev = %0.4f\n",
	    observed_fc_mean, observed_fc_stddev);
    printf("FC mean count = %lu  FC >= %0.4f = %lu  P(FC >= %0.4f) = %0.4f\n",
	    fc_mean_count, observed_fc_mean, fc_ge, observed_fc_mean,
	    (double)fc_ge / fc_mean_count);
}
