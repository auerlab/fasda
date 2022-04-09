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
    FILE        *feature_stream, *condition_stream[MAX_CONDITIONS];
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

	if ( (condition_stream[conditions] =
	      bl_sam_fopen(condition_files[conditions], "r",
	      SAMTOOLS_ARGS)) == NULL )
	{
	    fprintf(stderr, ": Could not open %s for read: %s.\n",
		    condition_files[conditions], strerror(errno));
	    return EX_NOINPUT;
	}
	
	bl_sam_skip_header(condition_stream[conditions]);
    }
    
    return diffanal(feature_stream, condition_stream, conditions);
}


int     diffanal(FILE *feature_stream, FILE *condition_stream[], int conditions)

{
    char        previous_feature_chrom[BL_CHROM_MAX_CHARS + 1],
		previous_alignment_chrom[MAX_CONDITIONS][BL_CHROM_MAX_CHARS + 1];
    bl_gff_t    feature;
    bl_sam_t    alignment;
    int         c, cmp;
    double      coverage[MAX_CONDITIONS];
    
    // FIXME: Currently only support two conditions
    printf("%d conditions.\n", conditions);
    
    /*
     *  Start simple: Count reads overlapping each position
     *  in the gene and compute the average
     *  depth = total bases / gene length.
     *  Thoroughly test and optimize, then explore more
     *  sophisticated depth algorithms.
     */

    for (c = 0; c < conditions; ++c)
	strlcpy(previous_alignment_chrom[c], "0", BL_CHROM_MAX_CHARS + 1);

    // FIXME: discard unnecessary fields to improve performance
    bl_gff_init(&feature);
    bl_sam_init(&alignment);
    strlcpy(previous_feature_chrom, "0", BL_CHROM_MAX_CHARS + 1);
    printf("%2s %-20s %-10s %-10s %-10s\n", "Ch", "Gene", "Condition1", "Condition2", "Fold-change");
    
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
			condition_stream[c], previous_alignment_chrom[c],
			&alignment) == BL_READ_OK) &&
		     (bl_sam_gff_cmp(&alignment, &feature) == 0) )
		{
		    /*
		     *  Now count coverage of all reads overlapping the gene.
		     *  Buffer reads in case some overlap the next gene as well.
		     */
		    
		    coverage[c] = count_coverage(&feature, &alignment,
					condition_stream[c],
					previous_alignment_chrom[c]);
		}
	    }
	    if ( (coverage[0] != 0.0) || (coverage[1] != 0.0) )
		printf("%2s %-20s %10.2f %10.2f %10.2f\n",
		       BL_GFF_SEQID(&feature), BL_GFF_FEATURE_NAME(&feature),
		       coverage[0], coverage[1], coverage[1] / coverage[0]);
	}
    }
    
    for (c = 2; c < conditions; ++c)
	bl_sam_fclose(condition_stream[c]);
    xt_fclose(feature_stream);
    
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

int     bl_gff_find_overlapping_alignment(
	    bl_gff_t *feature, FILE *alignment_stream,
	    char *previous_alignment_chrom,
	    bl_sam_t *alignment)

{
    int     status, cmp;
    
    while ( ((status = bl_sam_read(alignment, alignment_stream, SAM_MASK)) == BL_READ_OK) &&
	    (bl_sam_gff_cmp(alignment, feature) < 0) )
    {
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
		       FILE *sam_stream, char *previous_alignment_chrom)

{
    int64_t overlapping_bases = 0;
    double  coverage;
    int     cmp;
    
    // FIXME: Buffer overlapping alignments so they can also
    // be checked against the next gene.
    
    // alignment arg already contains first overlapping read
    // Loop through all reads overlapping the feature
    do
    {
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
    }   while ( (bl_sam_read(alignment, sam_stream, SAM_MASK) == BL_READ_OK)
		&& (bl_sam_gff_cmp(alignment, feature) == 0) );
    
    // Divide by length of gene
    coverage = (double)overlapping_bases /
		(BL_GFF_END(feature) - BL_GFF_START(feature) + 1);
    return coverage;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s features.gff3 \\\n", argv[0]);
    fprintf(stderr, "\tcondition1.[sam|bam|cram][.gz|.bz2|.xz] \\\n");
    fprintf(stderr, "\tcondition2.[sam|bam|cram][.gz|.bz2|.xz]\n");
    exit(EX_USAGE);
}
