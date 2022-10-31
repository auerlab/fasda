#include <stdio.h>
#include <math.h>
#include <xtend/math.h>


double  normal_cdf(double x, double mean, double stddev)

{
    return 0.5 * (1 + erf((x - mean) / (stddev * M_SQRT2)));
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *      Compute a p-value using the Mann-Whitney u-test, a.k.a. Wilcoxon
 *      rank sum test.  This test has good power and low FDR for sample
 *      sizes of at least 8.  It has low power if either sample sizes < 8 and
 *      this function will return a p-value of 1.0 in that case.
 *
 *      The formula given here is apparently wrong:
 *      https://support.minitab.com/en-us/minitab/18/help-and-how-to/statistics/nonparametrics/how-to/mann-whitney-test/methods-and-formulas/methods-and-formulas/
 *      Test site for p-values: https://www.statziki.com/Mannwhitneyu
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
 *  2022-05-14  Jason Bacon Begin
 ***************************************************************************/

double  mann_whitney_pval(double rep_counts1[], double rep_counts2[],
			   size_t num_reps1, size_t num_reps2)

{
    double  z, p, s12, w1, w2, w;
    size_t  c1, c2, n = num_reps1, m = num_reps2;
    //unsigned long   total = 0, ties = 0;

    // P-value cannot be computed with fewer than 8 samples
    if ( (n < 8) || (m < 8) )
	return 1.0;

    /*
    puts("\n\nCounts1:");
    for (c1 = 0; c1 < n; ++c1)
	printf("%f ", rep_counts1[c1]);
    puts("\n");
    printf("Counts2:\n");
    for (c2 = 0; c2 < n; ++c2)
	printf("%f ", rep_counts2[c2]);
    puts("\n");
    */
    
    for (c1 = 0, w1 = 0.0; c1 < n; ++c1)
    {
	for (c2 = 0; c2 < m; ++c2)
	{
	    /*
	     *  Ties are almost non-existent, except when count == 0.
	     *  In that case, they can be the majority of comparisons.
	     */
	    
	    if ( rep_counts1[c1] > rep_counts2[c2] )
		s12 = 1;
	    else if ( rep_counts1[c1] < rep_counts2[c2] )
		s12 = 0;
	    else
	    {
		//++ties;
		//fprintf(stderr, "%f == %f\n", rep_counts1[c1], rep_counts2[c2]);
		s12 = 0.5;
	    }
	    w1 += s12;
	    //++total;
	}
    }
    // fprintf(stderr, "Total = %lu  Ties = %lu\n", total, ties);
    w2 = m * n - w1;
    w = XT_MIN(w1, w2);
    
    // Paul A.: Minitab formula is wrong.  Numerator is just n*m/2 per
    // R source for wilcox.test().
    // FIXME: This normal approximation might not be stable for small
    // sample sizes, so maybe compute the p-value combinatorically.
    // (Determine total # of possible rank sums and how many are less
    // (or greater) than this one.
    z = (w - n * m / 2) / sqrt(n * m * (n + m + 1.0) / 12.0);
    p = 2.0 * normal_cdf(z, 0.0, 1.0);
    // printf("  z = %f  p = %f  p(-1.96) = %f\n", z, p, normal_cdf(-1.96, 0.0, 1.0));
    return p;
}

