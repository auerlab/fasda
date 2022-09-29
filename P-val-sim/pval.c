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
#include <time.h>

void    usage(char *argv[])

{
    fprintf(stderr, "Usage:   %s count1-mean count2-mean slop replicates\n", argv[0]);
    fprintf(stderr, "Example: %s 100 200 .2 3\n",argv[0]);
    exit(EX_USAGE);
}

int     main(int argc,char *argv[])

{
    int     c, c1, c2, c3, pairs = 0, fc_ge = 0, replicates = 3,
	    less = 0, more = 0, equal = 0, *counts,
	    count1_mean, count2_mean, fc_count;
    double  fc, observed_mean_fc, *fc_list, slop, mean_fc;

    if ( argc != 4 )
	usage(argv);
    
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
    for (c = 0, observed_mean_fc = 0.0; c < replicates; ++c)
    {
	counts[c] = count1_mean
		 + random() % (int)(count1_mean * slop * 2) - count1_mean * slop;
	counts[c + replicates] = count2_mean
		 + random() % (int)(count2_mean * slop * 2) - count2_mean * slop;
	fc = (double)counts[c + replicates] / counts[c];
	printf("%5d %5d %f\n", counts[c], counts[c + replicates], fc);
	observed_mean_fc += fc;
    }
    observed_mean_fc /= replicates;
    printf("Mean FC = %f\n", observed_mean_fc);
    
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
    fc_count = 0;
    for (c1 = 0; c1 < replicates * 2; ++c1)
	for (c2 = 0; c2 < replicates * 2; ++c2)
	{
	    if ( c1 != c2 )
	    {
		++pairs;
		fc_list[fc_count] = (double)counts[c1] / counts[c2];
		printf("%3d / %3d = %f\n", counts[c1], counts[c2],
			fc_list[fc_count]);
		if ( fc_list[fc_count] >= observed_mean_fc )
		    ++fc_ge;
		
		if ( counts[c1] < counts[c2] )
		    ++less;
		else if ( counts[c1] > counts[c2] )
		    ++more;
		else
		    ++equal;
		++fc_count;
	    }
	}
    
    printf("\n\nTotal pairings = %d  FC >= %f = %d  P(FC >= %f) = %f\n",
	    pairs, observed_mean_fc, fc_ge, observed_mean_fc, (double)fc_ge / pairs);
    
    printf("\nDistribution: less = %d  more = %d  equal = %d\n",
	    less, more, equal);

    /*
     *  Find all possible means of fold-change triplets
     *  Count the number of means >= observed mean
     *  FIXME: Hard-coded for 3 replicates.  Need to generalize.
     */
    
    pairs = less = more = equal = 0;
    for (c1 = 0; c1 < fc_count; ++c1)
    {
	for (c2 = c1 + 1; c2 < fc_count; ++c2)
	{
	    for (c3 = c2 + 1; c3 < fc_count; ++c3)
	    {
		mean_fc = (fc_list[c1] + fc_list[c2] + fc_list[c3]) / 3.0;
		// printf("%d %d %d %f\n", c1, c2, c3, fc);
			
		if ( mean_fc >= observed_mean_fc )
		    ++fc_ge;
		if ( fc_list[c1] < 1.0 )
		    ++less;
		else if ( fc_list[c1] > 1.0 )
		    ++more;
		else
		    ++equal;
		++pairs;
	    }
	}
    }
    
    printf("\n\nTotal pairings = %d  FC >= %f = %d  P(FC >= %f) = %f\n",
	    pairs, observed_mean_fc, fc_ge, observed_mean_fc,
	    (double)fc_ge / pairs);
    
    // Q: Do we need to satisfy a null hypothesis, e.g. Wilcoxon here?
    printf("\nDistribution: less = %d  more = %d  equal = %d\n",
	    less, more, equal);
    
    return EX_OK;
}
