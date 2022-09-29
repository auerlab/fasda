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
    int     c, c1, c2, pairs = 0, fc_ge = 0, replicates = 3,
	    less = 0, more = 0, equal = 0, *counts,
	    count1_mean, count2_mean;
    double  fc, mean_fc, slop;

    if ( argc != 5 )
	usage(argv);
    
    count1_mean = atoi(argv[1]);
    count2_mean = atoi(argv[2]);
    slop = atof(argv[3]);
    replicates = atoi(argv[4]);
    
    counts = malloc(replicates * 2 * sizeof(*counts));
    
    // Find all possible combinations of 3 out of 6 samples
    /*
    printf("6 choose 3 = %lu\n", xt_n_choose_k(6, 3));
    for (c1 = 0; c1 < 6; ++c1)
	for (c2 = c1 + 1; c2 < 6; ++c2)
	    for (c3 = c2 + 1; c3 < 6; ++c3)
		printf("%d %d %d\n", c1, c2, c3);
    */

    /*
     *  Generate random counts with FC around 2
     */
    printf("\ncount1 = %d +/- to up to %0.0f%%, count2 = %d +/- same\n",
	    count1_mean, slop * 100, count2_mean);
    puts("Cond1 Cond2");
    srandom(time(NULL));
    for (c = 0, mean_fc = 0.0; c < replicates; ++c)
    {
	counts[c] = count1_mean
		 + random() % (int)(count1_mean * slop * 2) - count1_mean * slop;
	counts[c + replicates] = count2_mean
		 + random() % (int)(count2_mean * slop * 2) - count2_mean * slop;
	fc = (double)counts[c + replicates] / counts[c];
	printf("%5d %5d %f\n", counts[c], counts[c + replicates], fc);
	mean_fc += fc;
    }
    mean_fc /= replicates;
    printf("Mean FC = %f\n", mean_fc);
    
    puts("\nPooled counts:");
    for (c1 = 0; c1 < replicates * 2; ++c1)
	printf("%3d\n", counts[c1]);
    
    /*
     *  Compute fold-change for every possible random pairing.
     *  Satisfies the null hypothesis P(n1 > n2) = P(n1 < N2)
     */
    puts("\nFold-change of every possible pairing of pooled samples:");
    puts("(Counts are not paired with themselves.)");
    for (c1 = 0; c1 < replicates * 2; ++c1)
	for (c2 = 0; c2 < replicates * 2; ++c2)
	{
	    if ( c1 != c2 )
	    {
		++pairs;
		fc = (double)counts[c1] / counts[c2];
		printf("%3d / %3d = %0.1f\n", counts[c1], counts[c2], fc);
		if ( fc >= mean_fc )
		    ++fc_ge;
		
		if ( counts[c1] < counts[c2] )
		    ++less;
		else if ( counts[c1] > counts[c2] )
		    ++more;
		else
		    ++equal;
	    }
	}
    
    printf("\n\nTotal pairings = %d  FC >= %f = %d  P(FC >= %f) = %f\n",
	    pairs, mean_fc, fc_ge, mean_fc, (double)fc_ge / pairs);
    
    printf("\nDistribution: less = %d  more = %d  equal = %d\n",
	    less, more, equal);
    return EX_OK;
}
