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
#include <sys/param.h>      // MIN(), MAX()
#include <xtend/file.h>
#include <xtend/string.h>   // strlcpy() on Linux
#include <biolibc/gff.h>
#include <biolibc/sam.h>
#include <biolibc/biostring.h>
#include "diffanal.h"

int     main(int argc,char *argv[])

{
    char        *features_file, *condition_files[MAX_CONDITIONS];
    FILE        *feature_stream, *sam_streams[MAX_CONDITIONS];
    int         conditions, c;
    
    if ( argc < 4 )
	usage(argv);
    features_file = argv[1];
    
    if ( (feature_stream = xt_fopen(features_file, "r")) == NULL )
    {
	fprintf(stderr, "diffanal: Could not open %s for read: %s.\n",
		features_file, strerror(errno));
	return EX_NOINPUT;
    }
    bl_gff_skip_header(feature_stream);

    for (c = 2, conditions = 0; c < argc; ++c, ++conditions)
    {
	condition_files[conditions] = argv[c];

	if ( (sam_streams[conditions] =
	      bl_sam_fopen(condition_files[conditions], "r",
	      SAMTOOLS_ARGS)) == NULL )
	{
	    fprintf(stderr, ": Could not open %s for read: %s.\n",
		    condition_files[conditions], strerror(errno));
	    return EX_NOINPUT;
	}
	
	bl_sam_skip_header(sam_streams[conditions]);
    }
    
    return diffanal(feature_stream, sam_streams, conditions);
}


int     diffanal(FILE *feature_stream, FILE *sam_streams[], int conditions)

{
    char        previous_feature_chrom[BL_CHROM_MAX_CHARS + 1],
		previous_alignment_chrom[MAX_CONDITIONS][BL_CHROM_MAX_CHARS + 1];
    bl_gff_t    feature;
    bl_sam_t    alignment;
    int         c, cmp;
    double      coverage[MAX_CONDITIONS];
    FILE        *buffer_streams[MAX_CONDITIONS];
    bl_alignment_stats_t    alignment_stats = BL_ALIGNMENT_STATS_INIT;
    
    /*
     *  Start simple: Count reads overlapping each position
     *  in the gene and compute the average
     *  depth = total bases / gene length.
     *  Thoroughly test and optimize, then explore more
     *  sophisticated depth algorithms.
     */

    bl_gff_init(&feature);
    bl_sam_init(&alignment);
    strlcpy(previous_feature_chrom, "0", BL_CHROM_MAX_CHARS + 1);
    for (c = 0; c < conditions; ++c)
    {
	strlcpy(previous_alignment_chrom[c], "0", BL_CHROM_MAX_CHARS + 1);
	if ( (buffer_streams[c] = tmpfile()) == NULL )
	{
	    fprintf(stderr, "diffanal: Cannot create temp file for SAM buffer.\n");
	    return EX_CANTCREAT;
	}
    }
    print_header(conditions);
    
    while ( bl_gff_read(&feature, feature_stream, GFF_MASK) == BL_READ_OK )
    {
	if ( strcmp(BL_GFF_TYPE(&feature), "gene") == 0 )
	{
	    for (c = 0; c < conditions; ++c)
		coverage[c] = 0;
	    
	    // Verify that features are properly sorted
	    cmp = bl_chrom_name_cmp(BL_GFF_SEQID(&feature), previous_feature_chrom);
	    if ( cmp < 0 )
	    {
		fprintf(stderr, "diffanal: Error: GFF3 chromosomes out of order: %s %s\n",
			previous_feature_chrom, BL_GFF_SEQID(&feature));
		return EX_DATAERR;
	    }
	    else if ( cmp > 0 )
		strlcpy(previous_feature_chrom, BL_GFF_SEQID(&feature), BL_CHROM_MAX_CHARS + 1);

	    for (c = 0; c < conditions; ++c)
	    {
		/*
		 *  Find first overlapping read for this gene.
		 */

		if ( (bl_gff_find_overlapping_alignment(&feature,
			sam_streams[c], buffer_streams[c],
			previous_alignment_chrom[c],
			&alignment, &alignment_stats) == BL_READ_OK) &&
		     (bl_sam_gff_cmp(&alignment, &feature) == 0) )
		{
		    /*
		     *  Now count coverage of all reads overlapping the gene.
		     *  Buffer reads in case some overlap the next gene as well.
		     */
		    
		    coverage[c] = count_coverage(&feature, &alignment,
					sam_streams[c],
					buffer_streams[c],
					previous_alignment_chrom[c],
					&alignment_stats);
		}
	    }
	    print_fold_change(&feature, coverage, conditions);
	}
    }
    
    for (c = 2; c < conditions; ++c)
    {
	bl_sam_fclose(sam_streams[c]);
	fclose(buffer_streams[c]);
    }
    xt_fclose(feature_stream);
    
    printf("\nTotal alignments processed:        %lu\n",
	    BL_ALIGNMENT_STATS_TOTAL(&alignment_stats));
    printf("Alignments overlapping a feature:  %lu\n",
	    BL_ALIGNMENT_STATS_OVERLAPPING(&alignment_stats));
    
    return EX_OK;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *      Should this be part of biolibc?  Not sure if it's generally useful
 *      enough and sufficiently difficult to just rewrite where needed.
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
 *  2022-04-06  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_find_overlapping_alignment(bl_gff_t *feature,
	    FILE *sam_stream, FILE *buffer_stream,
	    char *previous_alignment_chrom, bl_sam_t *alignment,
	    bl_alignment_stats_t *alignment_stats)

{
    int     status, cmp;

    // First search alignments buffered from previous gene
    // fprintf(stderr, "Checking buffered alignments...\n");
    while ( ((status = bl_sam_read(alignment, buffer_stream, SAM_MASK)) == BL_READ_OK) &&
	    (bl_sam_gff_cmp(alignment, feature) < 0) )
    {
	// fprintf(stderr, "buffered: %s %lu\n", BL_SAM_RNAME(alignment), BL_SAM_POS(alignment));
	// Verify that alignments are properly sorted
	cmp = bl_chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom);
	if ( cmp > 0 )
	    strlcpy(previous_alignment_chrom, BL_SAM_RNAME(alignment), BL_CHROM_MAX_CHARS + 1);
	else if ( cmp < 0 )
	{
	    fprintf(stderr, "diffanal, %s(): "
		    "Chromosomes out of order in SAM stream: %s %s\n",
		    __FUNCTION__, previous_alignment_chrom, BL_SAM_RNAME(alignment));
	    exit(EX_DATAERR);
	}
    }
    // fprintf(stderr, "Done with buffered.\n");
    
    // No overlapping alignments in temp buffer
    if ( status == BL_READ_EOF )
    {
	// fprintf(stderr, "Checking new alignments...\n");
	// Discard all buffered alignments.  If they didn't overlap this
	// gene, they won't overlap the next either, since they're sorted.
	rewind(buffer_stream);
	ftruncate(fileno(buffer_stream), 0);
	
	while ( ((status = bl_sam_read(alignment, sam_stream, SAM_MASK)) == BL_READ_OK) &&
		(bl_sam_gff_cmp(alignment, feature) < 0) )
	{
	    ++BL_ALIGNMENT_STATS_TOTAL(alignment_stats);
	    
	    // Verify that alignments are properly sorted
	    cmp = bl_chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom);
	    if ( cmp > 0 )
		strlcpy(previous_alignment_chrom, BL_SAM_RNAME(alignment), BL_CHROM_MAX_CHARS + 1);
	    else if ( cmp < 0 )
	    {
		fprintf(stderr, "diffanal, %s(): "
			"Chromosomes out of order in SAM stream: %s %s\n",
			__FUNCTION__, previous_alignment_chrom, BL_SAM_RNAME(alignment));
		exit(EX_DATAERR);
	    }
	}
    }
    return status;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *      Count coverage of all alignments overlapping the given feature.
 *      Currently using a simple average depth over all positions without
 *      regard for paired-end inner distance between the forward and
 *      reverse reads.  More sophisticated algorithms can be substituted
 *      as time goes by.
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
 *  2022-04-08  Jason Bacon Begin
 ***************************************************************************/

double  count_coverage(bl_gff_t *feature, bl_sam_t *alignment,
		       FILE *sam_stream, FILE *buffer_stream,
		       char *previous_alignment_chrom,
		       bl_alignment_stats_t *alignment_stats)

{
    int64_t overlapping_bases = 0;
    double  coverage;
    int     cmp, status;
    long    buffer_pos;
    
    buffer_pos = ftell(buffer_stream);
    
    /*
     *  "alignment" arg already contains first overlapping read
     *  Loop through all reads overlapping the feature
     */
    
    // Check alignments buffered from previous gene first
    do
    {
	// fprintf(stderr, "Counting buffered alignemnts...\n");
	cmp = bl_chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom);
	if ( cmp > 0 )
	    strlcpy(previous_alignment_chrom, BL_SAM_RNAME(alignment), BL_CHROM_MAX_CHARS + 1);
	else if ( cmp < 0 )
	{
	    fprintf(stderr, "diffanal: %s(): "
		    "Chromosomes out of order in SAM stream: %s %s\n",
		    __FUNCTION__, previous_alignment_chrom, BL_SAM_RNAME(alignment));
	    exit(EX_DATAERR);
	}
	overlapping_bases += bl_gff_sam_overlap(feature, alignment);
    }   while ( ((status = bl_sam_read(alignment, buffer_stream, SAM_MASK))
			    == BL_READ_OK)
		&& (bl_sam_gff_cmp(alignment, feature) == 0) );
    
    // If we ran out of buffered alignments, discard the old ones and
    // continue in primary SAM stream
    if ( status == BL_READ_EOF )
    {
	// fprintf(stderr, "Counting new alignemnts...\n");
	// Discard all buffered alignments.  If they didn't overlap this
	// gene, they won't overlap the next either, since they're sorted.
	rewind(buffer_stream);
	ftruncate(fileno(buffer_stream), 0);
	buffer_pos = 0;     // rewind() obsoletes ftell() above
	
	while ( (bl_sam_read(alignment, sam_stream, SAM_MASK) == BL_READ_OK)
		&& (bl_sam_gff_cmp(alignment, feature) == 0) )
	{
	    ++BL_ALIGNMENT_STATS_TOTAL(alignment_stats);
	    ++BL_ALIGNMENT_STATS_OVERLAPPING(alignment_stats);
	    cmp = bl_chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom);
	    if ( cmp > 0 )
		strlcpy(previous_alignment_chrom, BL_SAM_RNAME(alignment), BL_CHROM_MAX_CHARS + 1);
	    else if ( cmp < 0 )
	    {
		fprintf(stderr, "diffanal: %s(): "
			"Chromosomes out of order in SAM stream: %s %s\n",
			__FUNCTION__, previous_alignment_chrom, BL_SAM_RNAME(alignment));
		exit(EX_DATAERR);
	    }
	    overlapping_bases += bl_gff_sam_overlap(feature, alignment);
    
	    /*
	     *  Buffer new overlapping alignments so they can also be
	     *  checked against the next gene.  The number of overlapping
	     *  alignments is typically in the tens or hundreds, but
	     *  occasionally in the tens of thousands.  Buffering them in
	     *  an array therefore leads to highly variable memory use,
	     *  requiring wasteful memory allocations for batch jobs.
	     *  Using a temporary file trades a modest amount of
	     *  performance for consistently low memory use and simplicity.
	     */
	    bl_sam_write(alignment, buffer_stream, SAM_MASK);
	}
    }
    
    fseek(buffer_stream, buffer_pos, SEEK_SET);
    
    // Divide by length of gene
    coverage = (double)overlapping_bases /
		(BL_GFF_END(feature) - BL_GFF_START(feature) + 1);
    return coverage;
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
    
    printf("%2s %-15s", "Ch", "Gene");
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

void    print_fold_change(bl_gff_t *feature, double coverage[], int conditions)

{
    int     c1, c2;
    
    printf("%2s %-15s",
	   BL_GFF_SEQID(feature), BL_GFF_FEATURE_NAME(feature));
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


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s features.gff3 \\\n", argv[0]);
    fprintf(stderr, "\tcondition1.[sam|bam|cram][.gz|.bz2|.xz] \\\n");
    fprintf(stderr, "\tcondition2.[sam|bam|cram][.gz|.bz2|.xz]\n");
    exit(EX_USAGE);
}
