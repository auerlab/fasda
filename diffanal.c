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
#include <biolibc/gff.h>
#include <biolibc/sam.h>
#include <biolibc/biostring.h>

#define MAX_CONDITIONS  128
#define MAX_SEQ_LEN     1024
#define GFF_MASK        BL_GFF_FIELD_SEQID|BL_GFF_FIELD_TYPE|\
			BL_GFF_FIELD_START|BL_GFF_FIELD_END|\
			BL_GFF_FIELD_ATTRIBUTES
#define SAM_MASK        BL_SAM_FIELD_RNAME|BL_SAM_FIELD_POS

int     bl_gff_sam_cmp(bl_gff_t *feature, bl_sam_t *alignment);
int     bl_sam_gff_cmp(bl_sam_t *alignment, bl_gff_t *feature);
int     gff_find_overlapping_alignment(bl_gff_t *feature,
				     FILE *stream, bl_sam_t *alignment);
void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    char        *features_file, *condition_files[MAX_CONDITIONS],
		last_feature_chrom[BL_CHROM_MAX_CHARS + 1];
    FILE        *features_stream, *condition_streams[MAX_CONDITIONS];
    bl_gff_t    feature;
    bl_sam_t    alignment;
    int         conditions, c, ch, ch2, cmp;
    
    if ( argc < 4 )
	usage(argv);
    features_file = argv[1];
    
    if ( (features_stream = xt_fopen(features_file, "r")) == NULL )
    {
	fprintf(stderr, "diffanal: Could not open %s for read: %s.\n",
		features_file, strerror(errno));
	return EX_NOINPUT;
    }
    bl_gff_skip_header(features_stream);

    for (c = 2, conditions = 0; c < argc; ++c, ++conditions)
    {
	condition_files[conditions] = argv[c];

	if ( (condition_streams[conditions] =
		bl_sam_fopen(condition_files[conditions], "r")) == NULL )
	{
	    fprintf(stderr, ": Could not open %s for read: %s.\n",
		    condition_files[conditions], strerror(errno));
	    return EX_NOINPUT;
	}
	// FIXME: Add this function
	// bl_sam_skip_header(&condition_streams[conditions]);
	while ( (ch = getc(condition_streams[conditions])) == '@' )
	{
	    putchar(ch);
	    while ( (ch2 = getc(condition_streams[conditions])) != '\n' )
		putchar(ch2);
	    putchar('\n');
	}
	ungetc(ch, condition_streams[conditions]);
    }
    printf("%d conditions.\n", conditions);
    
    // FIXME: discard unnecessary fields to improve performance
    bl_gff_init(&feature);
    bl_sam_init(&alignment);
    strlcpy(last_feature_chrom, "0", BL_CHROM_MAX_CHARS + 1);
    while ( bl_gff_read(&feature, features_stream, GFF_MASK) == BL_READ_OK )
    {
	if ( strcmp(BL_GFF_TYPE(&feature), "gene") == 0 )
	{
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
		fprintf(stderr, "New chrom: %s\n", BL_GFF_SEQID(&feature));
		strlcpy(last_feature_chrom, BL_GFF_SEQID(&feature), BL_CHROM_MAX_CHARS + 1);
	    }
	    printf("%s %s %" PRId64 " %" PRId64 " %s\n",
		    BL_GFF_SEQID(&feature), BL_GFF_TYPE(&feature),
		    BL_GFF_START(&feature), BL_GFF_END(&feature),
		    BL_GFF_FEATURE_NAME(&feature));
	    
	    for (c = 0; c < conditions; ++c)
	    {
		/*
		 *  Start simple: Count reads overlapping each position
		 *  in the gene and compute the average
		 *  depth = total bases / gene length.
		 *  Thoroughly test and optimize, then explore more
		 *  sophisticated depth algorithms.
		 */
		
		// All we need is a count of overlapping bases
		// Buffer reads or features to capture all overlaps,
		// like ad2vcf
		// FIXME: Verify sort order of both genes and alignments
		if ( (gff_find_overlapping_alignment(&feature,
			condition_streams[c], &alignment) == BL_READ_OK) &&
		     (bl_sam_gff_cmp(&alignment, &feature) == 0) )
		{
		    printf("condition %d alignment %s %lu - %lu overlaps feature %s %lu - %lu\n",
			    c, BL_SAM_RNAME(&alignment), BL_SAM_POS(&alignment),
			    BL_SAM_POS(&alignment) + BL_SAM_SEQ_LEN(&alignment),
			    BL_GFF_SEQID(&feature), BL_GFF_START(&feature),
			    BL_GFF_END(&feature));
		    //getchar();
		}
	    }
	}
    }
    
    for (c = 2; c < conditions; ++c)
	bl_sam_fclose(condition_streams[c]);
    xt_fclose(features_stream);
    
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
	    fprintf(stderr, "New alignment chrom: %s\n", BL_SAM_RNAME(alignment));
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


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s features.gff3 \\\n", argv[0]);
    fprintf(stderr, "\tcondition1.[sam|bam|cram][.gz|.bz2|.xz] \\\n");
    fprintf(stderr, "\tcondition2.[sam|bam|cram][.gz|.bz2|.xz]\n");
    exit(EX_USAGE);
}
