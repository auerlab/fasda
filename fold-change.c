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
#include <xtend/file.h>
#include <xtend/dsv.h>
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
    int         c;
    double      coverage[MAX_CONDITIONS];
    dsv_line_t  dsv_line[MAX_CONDITIONS];
    char        *id, *end;
    
    print_header(conditions);
    
    while ( dsv_line_read(&dsv_line[0], condition_streams[0], "\t") != EOF )
    {
	id = DSV_LINE_FIELDS_AE(&dsv_line[0], 0);
	if ( strcmp(id, "target_id") == 0 )
	{
	    // Skip header line if present
	    for (c = 1; c < conditions; ++c)
		dsv_skip_rest_of_line(condition_streams[c]);
	    continue;
	}
	
	coverage[0] = strtof(DSV_LINE_FIELDS_AE(&dsv_line[0], 1), &end);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "fold-change: Invalid abundance on line %d: %f\n",
		    0, coverage[0]);
	    return EX_DATAERR;
	}
	
	for (c = 1; c < conditions; ++c)
	{
	    coverage[c] = 0;
	    
	    if ( dsv_line_read(&dsv_line[c], condition_streams[c], "\t") != '\n' )
	    {
		fprintf(stderr, "fold-change: Expected newline on condition stream %d\n", c);
		return EX_DATAERR;
	    }
	
	    /*
	     *  We should see the same transcript/gene ID on corresponding
	     *  lines from each abundances file.
	     */
	    
	    if ( strcmp(DSV_LINE_FIELDS_AE(&dsv_line[c], 0), id) != 0 )
	    {
		fprintf(stderr, "fold-change: Abundances files out of sync: %s %s\n",
			id, DSV_LINE_FIELDS_AE(&dsv_line[c], 0));
		return EX_DATAERR;
	    }

	    coverage[c] = strtof(DSV_LINE_FIELDS_AE(&dsv_line[c], 1), &end);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "fold-change: Invalid abundance on line %d: %f\n",
			c, coverage[c]);
		return EX_DATAERR;
	    }
	}
	print_fold_change(id, coverage, conditions);
    }    
    
    for (c = 0; c < conditions; ++c)
	xt_fclose(condition_streams[c]);
    
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Print header for genes, coverage, and fold-change
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
 *      Print coverage and fold-change stats for a given gene
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-09  Jason Bacon Begin
 ***************************************************************************/

void    print_fold_change(const char *id, double coverage[], int conditions)

{
    int     c1, c2;
    
    printf("%-30s", id);
    for (c1 = 0; c1 < conditions; ++c1)
	printf(" %6.2f", coverage[c1]);
    for (c1 = 0; c1 < conditions; ++c1)
    {
	for (c2 = c1 + 1; c2 < conditions; ++c2)
	{
	    if ( (coverage[c1] != 0.0) || (coverage[c2] != 0.0) )
		printf(" %7.2f", coverage[c2] / coverage[c1]);
	    else
		printf(" %7s", "*");
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

/*
double  mann_whitney_p_val()

{
    p = (w - n*(n + m + 1)/2) /
	sqrt(n * m (n + m + 1) / 12);
}
*/

void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s \\\n"
	    "\tabundances1.tsv[.gz|.bz2|.xz] abundances2.tsv[.gz|.bz2|.xz] \\\n"
	    "\t[abundances3.tsv[.gz|.bz2|.xz] ...]\n",
	    argv[0]);
    exit(EX_USAGE);
}
