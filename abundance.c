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
#include <unistd.h>         // ftruncate()
#include <sys/param.h>      // MIN(), MAX()
#include <xtend/file.h>
#include <xtend/string.h>   // strlcpy() on Linux
#include <biolibc/gff.h>
#include <biolibc/sam.h>
#include <biolibc/biostring.h>
#include "diffanal.h"

int     main(int argc,char *argv[])

{
    char        *features_file, *sam_files[MAX_FILE_COUNT],
		*abundance_files[MAX_FILE_COUNT], *p;
    FILE        *feature_stream, *sam_streams[MAX_FILE_COUNT],
		*abundance_streams[MAX_FILE_COUNT];
    int         file_count, c, flags;
    
    if ( argc < 3 )
	usage(argv);
    
    for (c = 1, flags = 0x0; (c < argc) && (argv[c][0] == '-'); ++c)
    {
	if ( strcmp(argv[c], "--show-gene-name") == 0 )
	    flags |= DIFFANAL_FLAG_SHOW_GENE;
	else if ( strcmp(argv[c], "--map-to-gene") == 0 )
	    flags |= DIFFANAL_FLAG_MAP_GENE;
	else
	    usage(argv);
    }
    features_file = argv[c++];
    
    if ( (feature_stream = xt_fopen(features_file, "r")) == NULL )
    {
	fprintf(stderr, "diffanal: Could not open %s for read: %s.\n",
		features_file, strerror(errno));
	return EX_NOINPUT;
    }
    bl_gff_skip_header(feature_stream);

    for (file_count = 0; c < argc; ++c, ++file_count)
    {
	sam_files[file_count] = argv[c];

	if ( (sam_streams[file_count] =
	      bl_sam_fopen(sam_files[file_count], "r",
	      SAMTOOLS_ARGS)) == NULL )
	{
	    fprintf(stderr, "abundance: Could not open %s for read: %s.\n",
		    sam_files[file_count], strerror(errno));
	    return EX_NOINPUT;
	}
	bl_sam_skip_header(sam_streams[file_count]);
	
	if ( (p = strrchr(sam_files[file_count], '/')) != NULL )
	    ++p;
	else
	    p = sam_files[file_count];
	abundance_files[file_count] = strdup(p);
	if ( abundance_files[file_count] == NULL )
	{
	    fprintf(stderr, "abundance: Could not strdup() abudance_files[%d]\n",
		    file_count);
	    return EX_UNAVAILABLE;
	}
	if ( (p = strstr(abundance_files[file_count], ".bam")) == NULL )
	    if ( (p = strstr(abundance_files[file_count], ".sam")) == NULL )
		if ( (p = strstr(abundance_files[file_count], ".cram")) == NULL )
		{
		    fprintf(stderr, "abundance: Expected sam|bam|cram: %s\n",
			    abundance_files[file_count]);
		    return EX_NOINPUT;
		}
	
	*p = '\0';
	strlcat(p, "-abundance.tsv", PATH_MAX);
	//fprintf(stderr, "output = %s\n", abundance_files[file_count]);
	if ( (abundance_streams[file_count] =
	      fopen(abundance_files[file_count], "w")) == NULL )
	{
	    fprintf(stderr, "abundance: Could not open %s.\n",
		    abundance_files[file_count]);
	    return EX_CANTCREAT;
	}
	fprintf(abundance_streams[file_count],
		"target_id\tlength\teff_length\test_counts\ttpm\n");
	fflush(abundance_streams[file_count]);
    }
    
    return abundance(feature_stream, sam_streams, abundance_streams,
		     file_count, flags);
}


int     abundance(FILE *feature_stream, FILE *sam_streams[],
		 FILE *abundance_streams[], int file_count, int flags)

{
    char        previous_feature_chrom[BL_CHROM_MAX_CHARS + 1],
		previous_alignment_chrom[MAX_FILE_COUNT][BL_CHROM_MAX_CHARS + 1],
		*feature_type = "mRNA";
    bl_gff_t    feature;
    bl_sam_t    alignment;
    int         c, cmp;
    double      counts = 0.0;
    FILE        *buffer_streams[MAX_FILE_COUNT];
    bl_alignment_stats_t    alignment_stats = BL_ALIGNMENT_STATS_INIT;
    
    /*
     *  Start simple: Count reads overlapping each position
     *  in the gene and compute the average
     *  depth = total bases / gene length.
     *  Thoroughly test and optimize, then explore more
     *  sophisticated depth algorithms.
     */

    if ( flags & DIFFANAL_FLAG_MAP_GENE )
	feature_type = "gene";
    bl_gff_init(&feature);
    bl_sam_init(&alignment);
    strlcpy(previous_feature_chrom, "0", BL_CHROM_MAX_CHARS + 1);
    for (c = 0; c < file_count; ++c)
    {
	strlcpy(previous_alignment_chrom[c], "0", BL_CHROM_MAX_CHARS + 1);
	if ( (buffer_streams[c] = tmpfile()) == NULL )
	{
	    fprintf(stderr, "diffanal: Cannot create temp file for SAM buffer.\n");
	    return EX_CANTCREAT;
	}
    }
    
    while ( bl_gff_read(&feature, feature_stream, GFF_MASK) == BL_READ_OK )
    {
	if ( strcmp(BL_GFF_TYPE(&feature), feature_type) == 0 )
	{
	    // Verify that features are properly sorted
	    cmp = bl_chrom_name_cmp(BL_GFF_SEQID(&feature),
				    previous_feature_chrom);
	    if ( cmp < 0 )
	    {
		fprintf(stderr, "diffanal: Error: GFF3 chromosomes out of order: %s %s\n",
			previous_feature_chrom, BL_GFF_SEQID(&feature));
		return EX_DATAERR;
	    }
	    else if ( cmp > 0 )
		strlcpy(previous_feature_chrom, BL_GFF_SEQID(&feature),
			BL_CHROM_MAX_CHARS + 1);

	    for (c = 0; c < file_count; ++c)
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
		     *  Now count counts of all reads overlapping the gene.
		     *  Buffer reads in case some overlap the next gene as well.
		     */
		    
		    counts = count_coverage(&feature, &alignment,
					sam_streams[c],
					buffer_streams[c],
					previous_alignment_chrom[c],
					&alignment_stats);
		}
		print_abundance(abundance_streams[c], &feature, counts, flags);
	    }
	}
    }
    
    for (c = 0; c < file_count; ++c)
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
    double  counts;
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
    counts = (double)overlapping_bases /
		(BL_GFF_END(feature) - BL_GFF_START(feature) + 1);
    return counts;
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
 *  2022-05-04  Jason Bacon Begin
 ***************************************************************************/

int     print_abundance(FILE *abundance_stream, bl_gff_t *feature,
			double counts, int flags)

{
    char            *id;
    unsigned long   len = 0;
    double          eff_len = 0.0, tpm = 0.0;

    if ( flags & DIFFANAL_FLAG_SHOW_GENE )
	id = BL_GFF_FEATURE_NAME(feature);
    else
	id = BL_GFF_FEATURE_ID(feature);
    
    //printf("%2s %-15s %6.2f", BL_GFF_SEQID(feature), id, counts[0]);
    fprintf(abundance_stream, "%-15s\t%lu\t%6.2f\t%6.2f\t%6.2f\n",
	    id, len, eff_len, counts, tpm);
    return 0;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s features.gff3 \\\n"
	    "\t[--show-gene-name] [--map-to-gene]\n"
	    "\tfile.[sam|bam|cram][.gz|.bz2|.xz] \\\n"
	    "\t[file.[sam|bam|cram][.gz|.bz2|.xz] ...]\n", argv[0]);
    exit(EX_USAGE);
}
