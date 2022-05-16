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
#include "fold-change.h"

int     main(int argc,char *argv[])

{
    char        *condition_files[MAX_CONDITIONS];
    FILE        *condition_streams[MAX_CONDITIONS];
    int         conditions, c;
    
    if ( argc < 3 )
	usage(argv);
    
    for (c = 1, conditions = 0; c < argc; ++c, ++conditions)
    {
	condition_files[conditions] = argv[c];

	if ( (condition_streams[conditions] =
	      xt_fopen(condition_files[conditions], "r")) == NULL )
	{
	    fprintf(stderr, ": Could not open %s for read: %s.\n",
		    condition_files[conditions], strerror(errno));
	    return EX_NOINPUT;
	}
    }
    
    return fold_change(condition_streams, conditions);
}


int     fold_change(FILE *condition_streams[], int conditions)

{
    size_t      condition;
    double      condition_counts[MAX_CONDITIONS];
    dsv_line_t  dsv_line[MAX_CONDITIONS];
    char        *id;
    int         delim;
    static double   *rep_counts[MAX_CONDITIONS];
    size_t          num_reps[MAX_CONDITIONS];
    
    print_header(conditions);

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
		
		if ( (condition > 0) &&
		     (strcmp(DSV_LINE_FIELDS_AE(&dsv_line[condition], 0),
			    id) != 0) )
		{
		    fprintf(stderr, "fold-change: Abundances files out of sync: %s %s\n",
			    id, DSV_LINE_FIELDS_AE(&dsv_line[condition], 0));
		    return EX_DATAERR;
		}
		// Save for comparison with next condition
		id = DSV_LINE_FIELDS_AE(&dsv_line[condition], 0);

		/*
		 *  Allocate array of counts for replicates in this condition
		 */
		
		if ( rep_counts[condition] == NULL )
		{
		    num_reps[condition] = DSV_LINE_NUM_FIELDS(&dsv_line[0]) - 1;
		    rep_counts[condition] =
			xt_malloc(num_reps[condition],
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
		
		condition_counts[condition] =
		    dsv_total_counts(&dsv_line[condition],
				     rep_counts[condition]);
	    }
	    else
	    {
		// EOF should be reached on all files at the same time
		// condition should be 0, check with getc() on 1 - N
	    }
	}
	
	// Output fold-change and p-value
	print_fold_change(id, condition_counts, conditions,
			  rep_counts, num_reps);
    }    
    
    for (condition = 0; condition < conditions; ++condition)
	xt_fclose(condition_streams[condition]);
    
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

void    print_header(int conditions)

{
    int     c1, c2;
    
    printf("%-30s", "Feature");
    for (c1 = 0; c1 < conditions; ++c1)
	printf(" %5s%d", "Cond", c1 + 1);
    for (c1 = 0; c1 < conditions; ++c1)
    {
	for (c2 = c1 + 1; c2 < conditions; ++c2)
	    printf("  FC %d-%d", c1 + 1, c2 + 1);
    }
    putchar('\n');
}


/***************************************************************************
 *  Description:
 *      Print count and fold-change stats for a given gene
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-09  Jason Bacon Begin
 ***************************************************************************/

void    print_fold_change(const char *id, double condition_counts[],
			  int conditions, double *rep_counts[],
			  size_t num_reps[])

{
    int     c1, c2;
    
    printf("%-30s", id);
    
    // Report average counts across all reps
    for (c1 = 0; c1 < conditions; ++c1)
	printf(" %6.2f", condition_counts[c1] / num_reps[c1]);
    
    for (c1 = 0; c1 < conditions; ++c1)
    {
	for (c2 = c1 + 1; c2 < conditions; ++c2)
	{
	    if ( (condition_counts[c1] != 0.0) || (condition_counts[c2] != 0.0) )
		printf(" %7.2f", condition_counts[c2] / condition_counts[c1]);
	    else
		printf(" %7s", "*");
	    
	    // Compute p-value
	    printf(" %f", mann_whitney_p_val(rep_counts[c1], rep_counts[c2],
					     num_reps[c1], num_reps[c2]));
	}
    }
    putchar('\n');
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *      https://support.minitab.com/en-us/minitab/18/help-and-how-to/statistics/nonparametrics/how-to/mann-whitney-test/methods-and-formulas/methods-and-formulas/
 
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

double  mann_whitney_p_val(double rep_counts1[], double rep_counts2[],
			   size_t num_reps1, size_t num_reps2)

{
    double  z, p, s12, w, k;
    size_t  c1, c2, n = num_reps1, m = num_reps2;
    
    for (c1 = 0, w = 0.0; c1 < n; ++c1)
    {
	for (c2 = 0; c2 < m; ++c2)
	{
	    s12 = rep_counts1[c1] > rep_counts2[c2] ? 1 :
		rep_counts1[c1] < rep_counts2[c2] ? 0 : 0.5;
	    w += s12;
	}
    }
    
    // Map to Minitab doc: n = num_reps1, m = num_reps2
    printf(" u = %f", w);
    z = (w - n * (n + m + 1.0) / 2.0) /
	 sqrt(n * m * (n + m + 1.0) / 12.0);
    
    k = MIN(w, n * (n + m + 1.0) - w);
    p = 0.0;    // FIXME
    return p;
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

double  dsv_total_counts(dsv_line_t *dsv_line, double rep_counts[])

{
    size_t  f;
    double  total_counts;
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
    return total_counts;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s \\\n"
	    "\tabundances1.tsv[.gz|.bz2|.xz] abundances2.tsv[.gz|.bz2|.xz] \\\n"
	    "\t[abundances3.tsv[.gz|.bz2|.xz] ...]\n",
	    argv[0]);
    exit(EX_USAGE);
}
