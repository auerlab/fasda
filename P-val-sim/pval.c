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
		      unsigned long replicates, double observed_fc_mean);
void    fc_mean_exact_p_val(double fc_list[], size_t fc_count,
			    size_t replicates,
			    double observed_fc_mean, double observed_fc_stddev,
			    double dist_fc_mean);

void    usage(char *argv[])

{
    fprintf(stderr, "Usage:   %s count1-mean count2-mean max_deviation\n", argv[0]);
    fprintf(stderr, "Example: %s 100 200 .2 3\n",argv[0]);
    exit(EX_USAGE);
}

int     main(int argc,char *argv[])

{
    unsigned long   c, c1, c2, pairs, replicates,
		    less, more, equal, *counts,
		    count1_mean, count2_mean, fc_count, samples, fc_ge,
		    half_fc_count;
    double  fc, observed_fc_mean, *fc_list, max_deviation, fc_sum,
	    dist_fc_mean, fc_var_sum, observed_fc_stddev;

    if ( argc != 4 )
	usage(argv);

    count1_mean = atoi(argv[1]);
    count2_mean = atoi(argv[2]);
    max_deviation = atof(argv[3]);
    replicates = 3; // atoi(argv[4]);
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
	printf("%5lu %5lu %0.3f\n", counts[c], counts[c + replicates], fc);
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
    for (c = 0; c < replicates; ++c)
    {
	fc = (double)counts[c + replicates] / counts[c];
	fc_var_sum += (fc - observed_fc_mean) * (fc - observed_fc_mean);
    }
    observed_fc_stddev = sqrt(fc_var_sum / replicates);
    printf("FC mean = %0.3f  FC stddev = %f\n",
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
	    printf("%2lu %3lu / %3lu = %0.3f\n", c,
		    counts[c1], counts[c2], fc_list[c]);
	    printf("%2lu %3lu / %3lu = %0.3f\n", c,
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
    pairs = c * 2;
    
    // Check for program bugs
    if ( c != half_fc_count )
    {
	printf("%lu != %lu\n", fc_count, c);
	return 1;
    }
    
    printf("\nTotal pairings = %lu  FC >= %0.3f = %lu  P(FC >= %0.3f) = %0.3f\n",
	    pairs, observed_fc_mean, fc_ge, observed_fc_mean, (double)fc_ge / pairs);
    dist_fc_mean = fc_sum / fc_count;
    printf("less + more + equal should be %lu.  FC mean should be slightly > 1.\n",
	    fc_count);
    printf("Distribution: less = %lu  more = %lu  equal = %lu  FC mean = %0.3f\n\n",
	    less, more, equal, dist_fc_mean);
    fc_mean_exact_p_val(fc_list, fc_count, replicates,
			observed_fc_mean, observed_fc_stddev, dist_fc_mean);
    return EX_OK;
}


/*
 *  Find all possible means of fold-change triplets
 *  Count the number of means >= observed mean
 *  FIXME: Hard-coded for 3 replicates.  Need to generalize.
 */

void    fc_mean_exact_p_val(double fc_list[], size_t fc_count,
			    size_t replicates,
			    double observed_fc_mean, double observed_fc_stddev,
			    double dist_fc_mean)

{
    size_t  c1, c2, c3;
    unsigned long   fc_mean_count, output_interval, fc_ge,
		    less, more, equal, c;
    double  fc_sum, fc_mean;

    fc_mean_count = xt_n_choose_k(fc_count, replicates);
    printf("%zu choose %zu = %lu\n", fc_count, replicates, fc_mean_count);
    output_interval = fc_mean_count / 10;
    
    puts("Sample of means of all possible combinations of 3 fold-changes:");
    puts("smpl c1 c2 c3   fc1   fc2   fc3 FC mean");
    c = less = more = equal = fc_sum = 0;

    fc_ge = 0;
    for (c1 = 0; c1 < fc_count; ++c1)
    {
	for (c2 = c1 + 1; c2 < fc_count; ++c2)
	{
	    for (c3 = c2 + 1; c3 < fc_count; ++c3)
	    {
		if ( c % output_interval == 0 )
		    printf("%4lu %2zu %2zu %2zu %0.3f %0.3f %0.3f %0.3f\n",
			fc_mean_count, c1, c2, c3,
			fc_list[c1], fc_list[c2], fc_list[c3], fc_mean);
		
		// FIXME: Don't think we need to check 1/FC
		fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3]) / replicates;
		if ( (fc_mean >= observed_fc_mean) )
		    ++fc_ge;
    
		// Debug
		if ( fc_mean < dist_fc_mean )
		    ++less;
		else if ( fc_mean > dist_fc_mean )
		    ++more;
		
		// Debug
		++c;
		fc_sum += fc_mean;
	    }
	}
    }
    
    printf("\nFC >= %0.3f from generated code = %lu\n", observed_fc_mean,
	    fc_ge_count(fc_list, fc_count, replicates, observed_fc_mean));
    
    // Check for program bugs
    if ( c != fc_mean_count )
    {
	fprintf(stderr, "fc_mean_count = %lu  c = %lu\n", fc_mean_count, c);
	exit(EX_SOFTWARE);
    }
    
    printf("Observed FC mean = %0.3f  Observed FC stddev = %0.3f\n",
	    observed_fc_mean, observed_fc_stddev);
    printf("FC mean count = %lu  FC >= %0.3f = %lu  P(FC >= %0.3f) = %0.3f\n",
	    fc_mean_count, observed_fc_mean, fc_ge, observed_fc_mean,
	    (double)fc_ge / fc_mean_count);

    equal = fc_mean_count - less - more;
    printf("Distribution: < %0.3f = %lu  > %0.3f = %lu  equal = %lu\n",
	    dist_fc_mean, less, dist_fc_mean, more, equal);
    
    // Debug: Should equal fc_mean
    printf("Mean of FC means = %0.3f\n", fc_sum / fc_mean_count);
}
