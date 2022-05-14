/***************************************************************************
 *  Description:
 *      Compute normalized abundances from raw counts provided in a TSV file
 *      similar to Kallisto abundances.tsv.  The TSV can be generated by
 *      the abundance subcommand from alignment data if necessary.
 *
 *  Raw yeast data with 48 biological replicates:
 *      https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/
 *      https://www.ebi.ac.uk/ena/browser/view/PRJEB5348 (FASTQ files)
 *      https://figshare.com/ndownloader/files/2194841 (sample map)
 *
 *  Arguments:

 *  Returns:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-05-04  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sysexits.h>
#include <stdlib.h>
#include <xtend/dsv.h>
#include <xtend/file.h>
#include "normalize.h"

int     main(int argc,char *argv[])

{
    int     arg;
    
    if ( argc < 3 )
	usage(argv);

    for (arg = 1; *argv[arg] == '-'; ++arg)
    {
	if ( strcmp(argv[arg], "--mrn") == 0 )
	    ;
	else
	    usage(argv);
    }
    
    return mrn(argc, argv, arg);
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
 *  2022-05-14  Jason Bacon Begin
 ***************************************************************************/

int     mrn(int argc, char *argv[], int arg)

{
    dsv_line_t  dsv_line;
    int         sample, sample_count, c;
    FILE        *abundance_streams[DIFFANAL_MAX_SAMPLES];
    
    for (sample_count = 0; arg < argc; ++sample_count, ++arg)
    {
	if ( (abundance_streams[sample_count] = xt_fopen(argv[arg], "r")) == NULL )
	{
	    fprintf(stderr, ": Could not open %s for read: %s.\n",
		    argv[arg], strerror(errno));
	    return EX_NOINPUT;
	}
    }
    
    dsv_line_init(&dsv_line);

    // Abundance file format:
    // target_id       length  eff_length      est_counts      tpm
    
    while ( ! feof(abundance_streams[0]) )
    {
	for (sample = 0; sample < sample_count; ++sample)
	{
	    if ( dsv_line_read(&dsv_line, abundance_streams[sample], "\t") == EOF )
	    {
		// Make sure all files reach EOF together
		for (c = 0; c < sample_count; ++c)
		    if ( getc(abundance_streams[c]) != EOF )
			fprintf(stderr, "normalize: EOF reached on %s but not %s\n",
				argv[sample + 1], argv[c + 1]);
		break;
	    }
	    else
	    {
		// Dummy output: Just echo non-normalized counts to test UI
		printf("%s\t%s\n", DSV_LINE_FIELDS_AE(&dsv_line, 0),
			DSV_LINE_FIELDS_AE(&dsv_line, 3));
		/*
		 *  Median of ratios normalization
		 *  https://scienceparkstudygroup.github.io/research-data-management-lesson/median_of_ratios_manual_normalization/index.html
		 *  Similar to TMM but more robust: doi 10.1093/bib/bbx008
		 *
		 *  First sweep:
		 *
		 *  Read raw counts for all genes and all samples
		 *
		 *  1.  Take log of every count (just for filtering in step 3?)
		 *  2.  Average of all samples for the gene (compute pseudo-reference)
		 *  3.  Remove genes witn -inf as average
		 *  4.  Subtract pseudo-reference from each log(expression)
		 *      This is actually a ratio since subtracting a log is dividing
		 *      We'll need to store this value and later sort to find median
		 *
		 *  After first sweep:
		 *
		 *  5.  Take the median of the ratios for each sample
		 *  6.  exp(median) = count scaling factor
		 *  7.  Divide counts by scaling factor
		 */
		
		
	    }
	}
    }
    
    for (sample = 0; sample < sample_count; ++sample)
	fclose(abundance_streams[sample]);
    
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s abundance1.tsv abundance2.tsv ...\n", argv[0]);
    exit(EX_USAGE);
}
