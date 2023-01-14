/***************************************************************************
 *  Description:
 *      Fast and easy differential analysis tools
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-04  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/param.h>
#include <xtend/file.h>
#include <xtend/dsv.h>
#include <xtend/mem.h>
#include <xtend/math.h>     // XT_MAX()
#include "fold-change.h"
#include "exact-p-val.h"

const int   Debug = 0;

int     main(int argc,char *argv[])

{
    char        *condition_files[MAX_CONDITIONS];
    FILE        *condition_streams[MAX_CONDITIONS],
		*diff_stream = stdout;
    int         conditions, arg;
    unsigned    flags = 0;
    
    if ( argc < 3 )
	usage(argv);

    for (arg = 1; *argv[arg] == '-'; ++arg)
    {
	if ( strcmp(argv[arg], "--output") == 0 )
	{
	    if ( (diff_stream = fopen(argv[++arg], "w")) == NULL )
	    {
		fprintf(stderr, "fold-change: Could not open %s for write: %s.\n",
			argv[arg], strerror(errno));
		return EX_CANTCREAT;
	    }
	}
	else if ( strcmp(argv[arg], "--near-exact") == 0 )
	{
	    flags |= FC_FLAG_NEAR_EXACT;
	}
    }
    
    for (conditions = 0; arg < argc; ++arg, ++conditions)
    {
	condition_files[conditions] = argv[arg];

	if ( (condition_streams[conditions] =
	      xt_fopen(condition_files[conditions], "r")) == NULL )
	{
	    fprintf(stderr, ": Could not open %s for read: %s.\n",
		    condition_files[conditions], strerror(errno));
	    return EX_NOINPUT;
	}
    }
    
    return fold_change(condition_streams, conditions, diff_stream, flags);
}


int     fold_change(FILE *condition_streams[], int conditions,
		    FILE *diff_stream, unsigned flags)

{
    size_t      condition, c, num_repls[MAX_REPLICATES];
    double      cond_tot_counts[MAX_CONDITIONS],
		condition_stddevs[MAX_CONDITIONS],
		*rep_counts[MAX_REPLICATES];    // 2D conditions x replicates
    dsv_line_t  dsv_line[MAX_REPLICATES];
    // Silence GCC 7 uninit warning, later versions OK
    char        *id = NULL, *new_id = NULL;
    int         delim;
    unsigned long    count = 0;
    
    if ( conditions < 2 )
    {
	fprintf(stderr, "fold-change: Must have at least 2 conditions.\n");
	return EX_DATAERR;
    }

    // Skip header line if present
    for (condition = 0; condition < conditions; ++condition)
    {
	if ( dsv_line_read(&dsv_line[condition],
	    condition_streams[condition], "\t") != '\n' )
	{
	    fprintf(stderr, "fold-change: Error reading first feature.\n");
	    return EX_DATAERR;
	}
	id = DSV_LINE_FIELDS_AE(&dsv_line[condition], 0);
	if ( strcmp(id, "target_id") == 0 )
	    dsv_skip_rest_of_line(condition_streams[condition]);
	else
	    rewind(condition_streams[condition]);
    
	// Flag need for malloc below
	rep_counts[condition] = NULL;
    }
    
    print_header(diff_stream, conditions);
    
    while ( !feof(condition_streams[0]) )
    {
	/*
	 *  Fold-change is ratio of total (or avg) normalized abundances
	 *  for each feature across two conditions
	 */
	
	// Get abundances for all replicates of one feature for all conditions
	for (condition = 0; condition < conditions; ++condition)
	{
	    // Read counts for all replicates of one feature
	    delim = dsv_line_read(&dsv_line[condition],
				  condition_streams[condition], "\t");
	    if ( delim != EOF )
	    {
		// Sanity check: incomplete input line
		if ( delim != '\n' )
		{
		    fprintf(stderr, "fold-change: Expected newline on condition stream %zu\n",
			    condition);
		    return EX_DATAERR;
		}
		
		/*
		 *  Sanity check: We should see the same transcript/gene ID
		 *  on corresponding lines from each abundances file.
		 */
		
		new_id = DSV_LINE_FIELDS_AE(&dsv_line[condition], 0);
		if ( (condition > 0) && (strcmp(new_id, id) != 0) )
		{
		    fprintf(stderr, "fold-change: Abundances files out of sync: %s %s\n",
			    id, DSV_LINE_FIELDS_AE(&dsv_line[condition], 0));
		    return EX_DATAERR;
		}
		
		// Save for comparison with next condition
		id = new_id;

		/*
		 *  Allocate array of counts for replicates in this condition
		 */
		
		if ( rep_counts[condition] == NULL )
		{
		    num_repls[condition] = DSV_LINE_NUM_FIELDS(&dsv_line[0]) - 1;
		    rep_counts[condition] =
			xt_malloc(num_repls[condition],
				  sizeof(*rep_counts[condition]));
		    if ( rep_counts[condition] == NULL )
		    {
			fprintf(stderr, "fold-change: Could not allocate rep_counts.\n");
			return EX_UNAVAILABLE;
		    }
		}
    
		/*
		 *  Get sum of normalized counts across all replicates.
		 *  Also get counts for each replicate for computing p-value.
		 */
		
		cond_tot_counts[condition] =
		    dsv_total_counts(&dsv_line[condition],
				     rep_counts[condition],
				     &condition_stddevs[condition]);
	    }
	    else
	    {
		// EOF should be reached on condition == 0 first and on
		// all files at the same time
		if ( condition == 0 )
		{
		    for (c = 1; c < conditions; ++c)
			if ( getc(condition_streams[c]) != EOF )
			{
			    fprintf(stderr, "fold-change: Expected EOF on condition %zu after EOF on condition 0.\n",
				    condition);
			    return EX_DATAERR;
			}
		    break;  // Avoid re-reading EOF on condition 1, etc.
		}
		else
		{
		    fprintf(stderr, "fold-change: Found EOF first on condition %zu, should be 0.\n",
			    condition);
		    return EX_DATAERR;
		}
	    }
	}
	
	// Output fold-change and p-value
	print_fold_change(diff_stream, id,
			  cond_tot_counts, condition_stddevs, conditions,
			  rep_counts, num_repls, flags);
	
	// Progress counter
	if ( ++count % 100 == 0 )
	    fprintf(stderr, "%lu\r", count);
    }    
    
    for (condition = 0; condition < conditions; ++condition)
	xt_fclose(condition_streams[condition]);
    xt_fclose(diff_stream);
    
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Print header for genes, count, and fold-change
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-09  Jason Bacon Begin
 ***************************************************************************/

void    print_header(FILE *diff_stream, int conditions)

{
    int     c1, c2;
    
    fprintf(diff_stream, "%-20s", "Feature");
    for (c1 = 0; c1 < conditions; ++c1)
	fprintf(diff_stream, " %6s%d", "MNC", c1 + 1);
    for (c1 = 0; c1 < conditions; ++c1)
	fprintf(diff_stream, " %5s%d", "SD/C", c1 + 1);
    for (c1 = 0; c1 < conditions; ++c1)
    {
	for (c2 = c1 + 1; c2 < conditions; ++c2)
	    fprintf(diff_stream,"  %4s  FC %d-%d  %5s",
		    "%Agr", c1 + 1, c2 + 1, "P-val");
    }
    putc('\n', diff_stream);
}


/***************************************************************************
 *  Description:
 *      Print count and fold-change stats for a given gene
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-09  Jason Bacon Begin
 ***************************************************************************/

void    print_fold_change(FILE *diff_stream, const char *id,
			  double cond_tot_counts[],
			  double condition_stddevs[], int conditions,
			  double *rep_counts[], size_t num_repls[],
			  unsigned flags)

{
    int     c1, c2;
    double  pval, avg_count;
    
    fprintf(diff_stream,"%-20s", id);
    
    // Report mean counts across all replicates for each condition
    for (c1 = 0; c1 < conditions; ++c1)
	fprintf(diff_stream," %7.1f",
		cond_tot_counts[c1] / num_repls[c1]);

    // Report mean counts / standard deviation as an estimate of consistency
    // in the counts across all replicates
    for (c1 = 0; c1 < conditions; ++c1)
    {
	avg_count = cond_tot_counts[c1] / num_repls[c1];
	fprintf(diff_stream," %6.1f",
		avg_count == 0 ? 0 : condition_stddevs[c1] / avg_count);
    }
    
    /*
     *  Report fold-change for each combination of conditions, i.e.
     *  1 vs 2, 1 vs 3, 2 vs 3.
     */
    
    for (c1 = 0; c1 < conditions; ++c1)
    {
	for (c2 = c1 + 1; c2 < conditions; ++c2)
	{
	    // Agreement: Is fold-change up/down in all replicates?
	    fprintf(diff_stream, "  %4u",
		    agreement(c1, c2, num_repls, rep_counts));
	    
	    // Fold-change
	    if ( (cond_tot_counts[c1] == 0.0) && (cond_tot_counts[c2] == 0.0) )
		fprintf(diff_stream," %7s", "*");
	    else
		fprintf(diff_stream," %7.2f", 
			cond_tot_counts[c2] / cond_tot_counts[c1]);

	    // P-value
	    if ( flags & FC_FLAG_NEAR_EXACT )
	    {
		if ( num_repls[c1] <= 12  )
		    pval = near_exact_pval(rep_counts[c1], rep_counts[c2],
					   num_repls[c1]);
		else
		{
		    fprintf(stderr, "Current limit for near-exact P-values is 12 replicates.\n");
		    exit(EX_USAGE);
		}
	    }
	    else if ( (num_repls[c1] >= 8) && (num_repls[c2] >= 8) )
		pval = mann_whitney_pval(rep_counts[c1], rep_counts[c2],
					 num_repls[c1], num_repls[c2]);
	    else
		pval = near_exact_pval(rep_counts[c1], rep_counts[c2],
				       num_repls[c1]);
	    
	    fprintf(diff_stream, "  %7.5f", pval);
	}
    }
    putc('\n', diff_stream);
}


/***************************************************************************
 *  Description:
 *  
 *  History: 
 *  Date        Name        Modification
 *  2022-11-22  Jason Bacon Begin
 ***************************************************************************/

unsigned agreement(int c1, int c2, size_t num_repls[], double *rep_counts[])

{
    unsigned    r, c1_higher, c2_higher;
    
    for (r = c1_higher = c2_higher = 0;
	 (r < num_repls[c1]) && (r < num_repls[c2]); ++r)
    {
	if ( rep_counts[c2][r] > rep_counts[c1][r] )
	    ++c2_higher;    // Up-regulated
	else
	    ++c1_higher;
    }
    
    // Add .5 to round rather than truncate when demoting to unsigned
    return (100.0 * XT_MAX(c1_higher, c2_higher) / r) + 0.5;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
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
 *  2022-05-16  Jason Bacon Begin
 ***************************************************************************/

double  dsv_total_counts(dsv_line_t *dsv_line, double rep_counts[],
			 double *condition_stddevs)

{
    size_t  f;
    double  total_counts, mean, sum_sq, variance;
    char    *end;
    
    // All but first field are counts
    for (f = 1, total_counts = 0.0; f < DSV_LINE_NUM_FIELDS(dsv_line); ++f)
    {
	rep_counts[f-1] = strtof(DSV_LINE_FIELDS_AE(dsv_line, f), &end);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "fold-change: Invalid abundance: %s\n",
		    DSV_LINE_FIELDS_AE(dsv_line, f));
	    return EX_DATAERR;
	}
	total_counts += rep_counts[f-1];
	//fprintf(stderr, "total_counts = %f\n", total_counts);
	//getchar();
    }
    
    mean = total_counts / DSV_LINE_NUM_FIELDS(dsv_line);
    for (f = 1, sum_sq = 0; f < DSV_LINE_NUM_FIELDS(dsv_line); ++f)
	sum_sq += (rep_counts[f-1] - mean) * (rep_counts[f-1] - mean);
    variance = sum_sq / DSV_LINE_NUM_FIELDS(dsv_line);
    *condition_stddevs = sqrt(variance);
    return total_counts;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s [--output file.txt]\\\n"
	    "\tabundances1.tsv[.gz|.bz2|.xz] abundances2.tsv[.gz|.bz2|.xz] \\\n"
	    "\t[abundances3.tsv[.gz|.bz2|.xz] ...]\n",
	    argv[0]);
    exit(EX_USAGE);
}
