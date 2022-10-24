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
#include <stdlib.h>
#include <sysexits.h>
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
	
	printf("Exact P-value = %0.3f  Mann-Whitney = %0.3f\n",
		near_exact_p_val(counts1, counts2, replicates),
		mann_whitney_p_val(counts1, counts2, replicates, replicates));
    }
    return EX_OK;
}
