/***************************************************************************
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
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

void    fc_mean_exact_p_val(double fc_list[], size_t fc_count,
			    size_t replicates,
			    double observed_fc_mean, double observed_fc_stddev,
			    double dist_fc_mean);

void    usage(char *argv[])

{
    fprintf(stderr, "Usage:   %s count1-mean count2-mean slop replicates\n", argv[0]);
    fprintf(stderr, "Example: %s 100 200 .2 3\n",argv[0]);
    exit(EX_USAGE);
}

int     main(int argc,char *argv[])

{
    int     c, c1, c2, pairs = 0, replicates = 3,
	    less = 0, more = 0, equal = 0, *counts,
	    count1_mean, count2_mean, fc_count, fc_ge;
    double  fc, observed_fc_mean, *fc_list, slop, fc_sum,
	    dist_fc_mean, fc_var_sum, observed_fc_stddev;

    if ( argc != 4 )
	usage(argv);

    for (replicates = 3; replicates <= 10; ++replicates)
    {
	fc_count = xt_n_choose_k(replicates * 2, 2);
	printf("FC count for %d replicates = %d ", replicates, fc_count);
	printf("FC means = %lu\n", xt_n_choose_k(fc_count, replicates));
    }
    
    count1_mean = atoi(argv[1]);
    count2_mean = atoi(argv[2]);
    slop = atof(argv[3]);
    replicates = 3; // atoi(argv[4]);
    
    counts = malloc(replicates * 2 * sizeof(*counts));
    
    /*
     *  Generate random counts with FC around 2
     */
    printf("\ncount1 = %d +/- to up to %0.0f%%, count2 = %d +/- same\n",
	    count1_mean, slop * 100, count2_mean);
    puts("Cond1 Cond2");
    srandom(time(NULL));
    for (c = 0, observed_fc_mean = 0.0; c < replicates; ++c)
    {
	counts[c] = count1_mean
		 + random() % (int)(count1_mean * slop * 2)
		 - count1_mean * slop;
	counts[c + replicates] = count2_mean
		 + random() % (int)(count2_mean * slop * 2)
		 - count2_mean * slop;
	fc = (double)counts[c + replicates] / counts[c];
	printf("%5d %5d %0.3f\n", counts[c], counts[c + replicates], fc);
	observed_fc_mean += fc;
    }
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
    
    puts("\nPooled counts:");
    for (c = 0; c < replicates * 2; ++c)
	printf("%3d\n", counts[c]);
    
    /*
     *  Compute fold-change for every possible pairing.
     *  Satisfies the null hypothesis P(n1 > n2) = P(n1 < N2)
     */
    puts("\nFold-change of every possible pairing of samples:");
    puts("(Counts are not paired with themselves.)");
    
    // FIXME: Slightly over-allocated since we filter out c1 == c2
    fc_list = malloc(replicates * 2 * replicates * 2 * sizeof(*fc_list));
    fc_count = fc_sum = fc_ge = 0;
    for (c1 = 0; c1 < replicates * 2; ++c1)
	for (c2 = c1 + 1; c2 < replicates * 2; ++c2)
	{
	    ++pairs;
	    fc_list[fc_count] = (double)counts[c1] / counts[c2];
	    printf("%2d %3d / %3d = %0.3f\n", fc_count,
		    counts[c1], counts[c2], fc_list[fc_count]);
	    if ( fc_list[fc_count] >= observed_fc_mean )
		++fc_ge;
	    
	    if ( fc_list[fc_count] < 1.0 )
		++less;
	    else if ( fc_list[fc_count] > 1.0 )
		++more;
	    else
		++equal;
	    fc_sum += fc_list[fc_count];
	    ++fc_count;
	}
    
    printf("\nTotal pairings = %d  FC >= %0.3f = %d  P(FC >= %0.3f) = %0.3f\n",
	    pairs, observed_fc_mean, fc_ge, observed_fc_mean, (double)fc_ge / pairs);
    
    dist_fc_mean = fc_sum / fc_count;
    printf("Distribution: less = %d  more = %d  equal = %d  FC mean = %0.3f\n\n",
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
		    less, more, equal;
    double  fc_sum, fc_mean;

    fc_mean_count = xt_n_choose_k(fc_count, replicates);
    printf("%zu choose %zu = %lu\n", fc_count, replicates, fc_mean_count);
    output_interval = fc_mean_count / 10;
    
    puts("Sample of means of all possible combinations of 3 fold-changes:");
    puts("smpl c1 c2 c3   fc1   fc2   fc3 FC mean");
    fc_mean_count = less = more = equal = fc_sum = 0;

    fc_ge = 0;
    for (c1 = 0; c1 < fc_count; ++c1)
    {
	for (c2 = c1 + 1; c2 < fc_count; ++c2)
	{
	    for (c3 = c2 + 1; c3 < fc_count; ++c3)
	    {
		fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3]) / 3.0;
		if ( fc_mean_count % output_interval == 0 )
		    printf("%4lu %2zu %2zu %2zu %0.3f %0.3f %0.3f %0.3f\n",
			fc_mean_count, c1, c2, c3,
			fc_list[c1], fc_list[c2], fc_list[c3], fc_mean);
		if ( (fc_mean >= observed_fc_mean) ||
		     (fc_mean <= 1 / observed_fc_mean) )
		    ++fc_ge;
		
		if ( fc_mean < dist_fc_mean )   // Allow for round-off
		    ++less;
		else if ( fc_mean > dist_fc_mean )
		    ++more;
		else
		    ++equal;
		++fc_mean_count;
		fc_sum += fc_mean;
	    }
	}
    }
    
    printf("\nObserved FC mean = %0.3f  Observed FC stddev = %0.3f\n",
	    observed_fc_mean, observed_fc_stddev);
    printf("FC mean count = %lu  FC >= %0.3f = %lu  P(FC >= %0.3f) = %0.3f\n",
	    fc_mean_count, observed_fc_mean, fc_ge, observed_fc_mean,
	    (double)fc_ge / fc_mean_count);
    
    // Q: Do we need to satisfy a null hypothesis, e.g. Wilcoxon here?
    // FIXME: Why are we getting about twice as many > 1 as < 1
    printf("Distribution: < %0.3f = %lu  > %0.3f = %lu  equal = %lu\n",
	    dist_fc_mean, less, dist_fc_mean, more, equal);
    printf("Mean of FC means = %0.3f\n", fc_sum / fc_mean_count);
}
