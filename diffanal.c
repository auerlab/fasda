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
#include <biolibc/gff.h>
#include <biolibc/sam.h>
#include <biolibc/biostring.h>
#include "diffanal.h"

int     main(int argc,char *argv[])

{
    char        *features_file, *condition_files[MAX_CONDITIONS],
		last_feature_chrom[BL_CHROM_MAX_CHARS + 1];
    FILE        *feature_stream, *condition_stream[MAX_CONDITIONS];
    bl_gff_t    feature;
    bl_sam_t    alignment;
    int         conditions, c, ch, ch2, cmp;
    double      coverage[MAX_CONDITIONS];
    
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
		bl_sam_fopen(condition_files[conditions], "r")) == NULL )
	{
	    fprintf(stderr, ": Could not open %s for read: %s.\n",
		    condition_files[conditions], strerror(errno));
	    return EX_NOINPUT;
	}
	// FIXME: Add this function
	// bl_sam_skip_header(&condition_stream[conditions]);
	while ( (ch = getc(condition_stream[conditions])) == '@' )
	{
	    while ( (ch2 = getc(condition_stream[conditions])) != '\n' )
		;
	}
	ungetc(ch, condition_stream[conditions]);
    }
    printf("%d conditions.\n", conditions);
    
    /*
     *  Start simple: Count reads overlapping each position
     *  in the gene and compute the average
     *  depth = total bases / gene length.
     *  Thoroughly test and optimize, then explore more
     *  sophisticated depth algorithms.
     */
		
    // FIXME: discard unnecessary fields to improve performance
    bl_gff_init(&feature);
    bl_sam_init(&alignment);
    strlcpy(last_feature_chrom, "0", BL_CHROM_MAX_CHARS + 1);
    printf("%-20s %-10s %-10s %-10s\n", "Gene", "Condition1", "Condition2", "Fold-change");
    while ( bl_gff_read(&feature, feature_stream, GFF_MASK) == BL_READ_OK )
    {
	if ( strcmp(BL_GFF_TYPE(&feature), "gene") == 0 )
	{
	    coverage[0] = coverage[1] = 0;
	    cmp = bl_chrom_name_cmp(BL_GFF_SEQID(&feature), last_feature_chrom);
	    // FIXME: Make bl_chrom_name_cmp() return exactly 1 or -1 for
	    // adjacent chromosomes?
	    if ( cmp < 0 )
	    {
		fprintf(stderr, "diffanal: Error: GFF3 chromosomes out of order: %s %s\n",
			last_feature_chrom, BL_GFF_SEQID(&feature));
		return EX_DATAERR;
	    }
	    else if ( cmp > 0 )
	    {
		//fprintf(stderr, "New chrom: %s\n", BL_GFF_SEQID(&feature));
		strlcpy(last_feature_chrom, BL_GFF_SEQID(&feature), BL_CHROM_MAX_CHARS + 1);
	    }
	    /*
	    printf("%s %s %" PRId64 " %" PRId64 " %s\n",
		    BL_GFF_SEQID(&feature), BL_GFF_TYPE(&feature),
		    BL_GFF_START(&feature), BL_GFF_END(&feature),
		    BL_GFF_FEATURE_NAME(&feature));
	    */
	    
	    for (c = 0; c < conditions; ++c)
	    {
		/*
		 *  Find first overlapping read for this gene.
		 */
		
		// FIXME: Verify sort order of both genes and alignments
		if ( (gff_find_overlapping_alignment(&feature,
			condition_stream[c], &alignment) == BL_READ_OK) &&
		     (bl_sam_gff_cmp(&alignment, &feature) == 0) )
		{
		    /*
		    printf("condition %d alignment %s %lu - %lu overlaps feature %s %lu - %lu\n",
			    c, BL_SAM_RNAME(&alignment), BL_SAM_POS(&alignment),
			    BL_SAM_POS(&alignment) + BL_SAM_SEQ_LEN(&alignment),
			    BL_GFF_SEQID(&feature), BL_GFF_START(&feature),
			    BL_GFF_END(&feature));
		    */
		    //getchar();
		    
		    /*
		     *  Now count coverage of all reads overlapping the gene.
		     *  Buffer reads in case some overlap the next gene as well.
		     */
		    
		    coverage[c] = count_coverage(&feature, &alignment, condition_stream[c]);
		}
	    }
	    if ( (coverage[0] != 0.0) || (coverage[1] != 0.0) )
		printf("%-20s %10.2f %10.2f %10.2f\n",
		       BL_GFF_FEATURE_NAME(&feature), coverage[0], coverage[1],
		       coverage[1] / coverage[0]);
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

int     bl_sam_gff_cmp(bl_sam_t *alignment, bl_gff_t *feature)

{
    int     status = bl_chrom_name_cmp(BL_SAM_RNAME(alignment),
					    BL_GFF_SEQID(feature));
    
    if ( status != 0 )
	// Different chromosomes
	return status;
    else if ( BL_SAM_POS(alignment) + BL_SAM_SEQ_LEN(alignment) - 1
		< BL_GFF_START(feature) )
	// Alignment ends before the start of feature
	return -1;
    else if ( BL_SAM_POS(alignment) > BL_GFF_END(feature) )
	// Alignment starts after the end of feature
	return 1;
    else
	// Overlap
	return 0;
}


int     bl_gff_sam_cmp(bl_gff_t *feature, bl_sam_t *alignment)

{
    return -bl_sam_gff_cmp(alignment, feature);
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
 *  2022-04-06  Jason Bacon Begin
 ***************************************************************************/

int     gff_find_overlapping_alignment(bl_gff_t *feature,
				       FILE *stream, bl_sam_t *alignment)

{
    int     status;
    static char    last_chrom[BL_CHROM_MAX_CHARS + 1] = "0";
    
    while ( ((status = bl_sam_read(alignment, stream, SAM_MASK)) == BL_READ_OK) &&
	    (bl_sam_gff_cmp(alignment, feature) < 0) )
    {
	if ( bl_chrom_name_cmp(last_chrom, BL_SAM_RNAME(alignment)) != 0 )
	{
	    //fprintf(stderr, "New alignment chrom: %s\n", BL_SAM_RNAME(alignment));
	    strlcpy(last_chrom, BL_SAM_RNAME(alignment), BL_CHROM_MAX_CHARS + 1);
	}
    }
    /*
	printf("%s %lu %s %lu\n",
		BL_SAM_RNAME(alignment), BL_SAM_POS(alignment),
		BL_GFF_SEQID(feature), BL_GFF_START(feature));
    */
    return status;
}


double  count_coverage(bl_gff_t *feature, bl_sam_t *alignment, FILE *sam_stream)

{
    int64_t overlapping_bases = 0;
    double  coverage;
    
    // alignment arg already contains first overlapping read
    // Loop through all reads overlapping the feature
    do
    {
	overlapping_bases += bl_gff_sam_overlap(feature, alignment);
    }   while ( (bl_sam_read(alignment, sam_stream, SAM_MASK) == BL_READ_OK)
		&& (bl_sam_gff_cmp(alignment, feature) == 0) );
    
    // Divide by length of gene
    coverage = (double)overlapping_bases /
		(BL_GFF_END(feature) - BL_GFF_START(feature) + 1);
    return coverage;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *      |-----------------------|
 *          |-----------------------------------|
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
 *  2022-04-07  Jason Bacon Begin
 ***************************************************************************/

int64_t bl_gff_sam_overlap(bl_gff_t *feature, bl_sam_t *alignment)

{
    int64_t alignment_end = BL_SAM_POS(alignment) + BL_SAM_SEQ_LEN(alignment),
	    overlap_start = MAX(BL_GFF_START(feature), BL_SAM_POS(alignment)),
	    overlap_end = MIN(BL_GFF_END(feature), alignment_end);
    
    //fprintf(stderr, "%" PRId64 " %" PRId64 "\n", overlap_start, overlap_end);
    //fprintf(stderr, "Coverage = %" PRId64 "\n", overlap_end - overlap_start + 1);
    return overlap_end - overlap_start + 1;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s features.gff3 \\\n", argv[0]);
    fprintf(stderr, "\tcondition1.[sam|bam|cram][.gz|.bz2|.xz] \\\n");
    fprintf(stderr, "\tcondition2.[sam|bam|cram][.gz|.bz2|.xz]\n");
    exit(EX_USAGE);
}
