/***************************************************************************
 *  Description:
 *      Fast and easy differential analysis tools
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-04  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <stdbool.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>         // ftruncate()
#include <limits.h>         // PATH_MAX OpenIndiana
#include <stdbool.h>
#include <xtend/file.h>
#include <xtend/string.h>   // strlcpy() on Linux
#include <biolibc/gff.h>
#include <biolibc/sam.h>
#include <biolibc/biostring.h>
#include "abundance.h"

bool    Debug = false;

int     main(int argc,char *argv[])

{
    char        *features_file,
		*sam_files[MAX_FILE_COUNT],
		*abundance_files[MAX_FILE_COUNT],
		*p,
		*feature_type = "mRNA",
		*output_dir = NULL;
    FILE        *feature_stream, *sam_streams[MAX_FILE_COUNT],
		*abundance_streams[MAX_FILE_COUNT];
    int         file_count, c, flags;
    
    if ( argc < 3 )
	usage(argv);

    // FIXME: Add --require-flags, --include-flags, --exclude-flags
    // See samtools-view
    for (c = 1, flags = 0x0; (c < argc) && (argv[c][0] == '-'); ++c)
    {
	if ( strcmp(argv[c], "--debug") == 0 )
	    Debug = true;
	else if ( strcmp(argv[c], "--show-gene-name") == 0 )
	    flags |= FASDA_FLAG_SHOW_GENE;
	else if ( strcmp(argv[c], "--ignore-chromosome-order") == 0 )
	    flags |= FASDA_IGNORE_CHR_ORDER;
	else if ( strcmp(argv[c], "--feature-type") == 0 )
	    feature_type = argv[++c];
	else if ( strcmp(argv[c], "--output-dir") == 0 )
	    output_dir = argv[++c];
	else
	    usage(argv);
    }
    features_file = argv[c++];
    
    if ( (feature_stream = xt_fopen(features_file, "r")) == NULL )
    {
	fprintf(stderr, "fasda: Could not open %s for read: %s.\n",
		features_file, strerror(errno));
	return EX_NOINPUT;
    }
    bl_gff_skip_header(feature_stream);

    for (file_count = 0; c < argc; ++c, ++file_count)
    {
	/*
	 *  Open SAM/BAM/CRAM input file
	 */
	
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
	
	/*
	 *  Open abundance output file for this input
	 */
	
	abundance_files[file_count] = malloc(PATH_MAX + 1);
	if ( abundance_files[file_count] == NULL )
	{
	    fprintf(stderr, "abundance: Could not strdup() abudance_files[%d]\n",
		    file_count);
	    return EX_UNAVAILABLE;
	}
	
	if ( output_dir != NULL )
	{
	    // Get base name of SAM input
	    if ( (p = strrchr(sam_files[file_count], '/')) != NULL )
		++p;
	    else
		p = sam_files[file_count];
	    
	    snprintf(abundance_files[file_count], PATH_MAX + 1, "%s/%s",
		    output_dir, p);
	}
	else
	    strlcpy(abundance_files[file_count], sam_files[file_count], PATH_MAX +1);
	
	// Replace .sam/.bam/.cram with -abundance.tsv
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
	fprintf(stderr, "Writing abundances to %s\n",
	    abundance_files[file_count]);
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
		     file_count, feature_type, flags);
}


int     abundance(FILE *feature_stream, FILE *sam_streams[],
		 FILE *abundance_streams[], int file_count,
		 const char *feature_type, int flags)

{
    char        previous_feature_chrom[BL_CHROM_MAX_CHARS + 1],
		previous_alignment_chrom[MAX_FILE_COUNT][BL_CHROM_MAX_CHARS + 1],
		*id, *p;
    bl_gff_t    feature, subfeature;
    bl_sam_t    alignment;
    int         c, cmp, status;
    int64_t     length, previous_start;
    double      est_counts = 0.0;
    bool        alternate_exons;
    FILE        *buffer_streams[MAX_FILE_COUNT];
    bl_alignment_stats_t    alignment_stats = BL_ALIGNMENT_STATS_INIT;
    
    bl_gff_init(&feature);
    bl_gff_init(&subfeature);
    bl_sam_init(&alignment);

    strlcpy(previous_feature_chrom, "0", BL_CHROM_MAX_CHARS + 1);
    for (c = 0; c < file_count; ++c)
    {
	strlcpy(previous_alignment_chrom[c], "0", BL_CHROM_MAX_CHARS + 1);
	if ( (buffer_streams[c] = tmpfile()) == NULL )
	{
	    fprintf(stderr, "fasda: Cannot create temp file for SAM buffer.\n");
	    return EX_CANTCREAT;
	}
    }
    
    while ( bl_gff_read(&feature, feature_stream, GFF_MASK) == BL_READ_OK )
    {
	// Some IDs contain aliases separated by '|'.  The last is usually
	// the unique ID.
	id = BL_GFF_FEATURE_ID(&feature);
	if ( (id != NULL) && (p = strrchr(id, '|')) != NULL )
	    id = p + 1;

	if ( false ) // Debug )
	    fprintf(stderr, "%s %s %s\n",
		    id, BL_GFF_TYPE(&feature), feature_type);
	if ( strcmp(BL_GFF_TYPE(&feature), feature_type) == 0 )
	{
	    if ( Debug )
		fprintf(stderr, "New %s: %s\n",
			feature_type, id);
	    
	    // Verify that features are properly sorted
	    cmp = chrom_name_cmp(BL_GFF_SEQID(&feature),
				    previous_feature_chrom, flags);
	    if ( cmp > 0 )
		strlcpy(previous_feature_chrom, BL_GFF_SEQID(&feature),
			BL_CHROM_MAX_CHARS + 1);
	    else if ( cmp < 0 )
	    {
		fprintf(stderr, "fasda: Error: GFF3 chromosomes out of order: %s %s\n",
			previous_feature_chrom, BL_GFF_SEQID(&feature));
		return EX_DATAERR;
	    }

	    // Compute sum of exon lengths for abundance TSV output
	    length = 0;
	    previous_start = 0;
	    alternate_exons = false;
	    
	    if ( Debug )
		fprintf(stderr, "Finding exons of %s...\n", id);
	    while ( ((status = bl_gff_read(&subfeature, feature_stream, GFF_MASK))
			== BL_READ_OK)
		    && (strcmp(BL_GFF_TYPE(&subfeature), "###") != 0)
		    && (strcmp(BL_GFF_TYPE(&subfeature), feature_type) != 0) )
	    {
		// fprintf(stderr, "type = %s\n", BL_GFF_TYPE(&subfeature));
		// Stop counting exons as soon as one has a lower start
		// position than the previous one.  This means we're into
		// alternate splicing data.  Lengths computed this way match
		// kallisto's abundance.tsv.
		if ( ! alternate_exons && 
		     (strcmp(BL_GFF_TYPE(&subfeature), "exon") == 0) )
		{
		    //fprintf(stderr, "%ld %ld\n",
		    //        previous_start, BL_GFF_START(&subfeature));
		    if ( BL_GFF_START(&subfeature) > previous_start )
		    {
			if ( Debug )
			    fprintf(stderr, "Adding exon %" PRId64 " %" PRId64 " ",
				    BL_GFF_START(&subfeature), BL_GFF_END(&subfeature));
			length += BL_GFF_END(&subfeature) -
			    BL_GFF_START(&subfeature) + 1;
			if ( Debug )
			    fprintf(stderr, "length = %" PRId64 "\n", length);
			previous_start = BL_GFF_START(&subfeature);
		    }
		    else
			alternate_exons = true;
		}
	    }
	    
	    if ( Debug )
	    {
		fprintf(stderr, "Found %s  length = %" PRId64 " status = %d\n",
			feature_type, length, status);
		fprintf(stderr, "Rewinding to %zu (%s)...\n",
			BL_GFF_FILE_POS(&subfeature),
			BL_GFF_TYPE(&subfeature));
	    }
	    
	    // Rewind to beginning of mRNA/transcript/gene so outer loop reads it
	    // FIXME: Maybe don't use BL_GFF_FILE_POS()
	    if ( strcmp(BL_GFF_TYPE(&subfeature), "###") != 0 )
		fseeko(feature_stream, BL_GFF_FILE_POS(&subfeature), SEEK_SET);
	    
	    for (c = 0; c < file_count; ++c)
	    {
		/*
		 *  Find first overlapping read for this gene.
		 */

		if ( (bl_gff_find_overlapping_alignment(&feature,
			sam_streams[c], buffer_streams[c],
			previous_alignment_chrom[c],
			&alignment, &alignment_stats, flags) == BL_READ_OK) &&
		     (bl_sam_gff_cmp(&alignment, &feature) == 0) )
		{
		    /*
		     *  Now get counts of all reads overlapping the gene.
		     *  Buffer reads in case some overlap the next gene as well.
		     */
		    
		    est_counts = count_coverage(&feature, &alignment,
					sam_streams[c],
					buffer_streams[c],
					previous_alignment_chrom[c],
					&alignment_stats, flags);
		}
		print_abundance(abundance_streams[c], &feature, length,
				est_counts, flags);
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
	    bl_alignment_stats_t *alignment_stats, int flags)

{
    int     status, cmp;

    // First search alignments buffered from previous gene
    // fprintf(stderr, "Checking buffered alignments...\n");
    while ( ((status = bl_sam_read(alignment, buffer_stream, SAM_MASK)) == BL_READ_OK) &&
	    (bl_sam_gff_cmp(alignment, feature) < 0) )
    {
	// fprintf(stderr, "buffered: %s %lu\n", BL_SAM_RNAME(alignment), BL_SAM_POS(alignment));
	// Verify that alignments are properly sorted
	cmp = chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom, flags);
	if ( cmp > 0 )
	    strlcpy(previous_alignment_chrom, BL_SAM_RNAME(alignment), BL_CHROM_MAX_CHARS + 1);
	else if ( cmp < 0 )
	{
	    fprintf(stderr, "fasda, %s(): "
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
	    cmp = chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom, flags);
	    if ( cmp > 0 )
		strlcpy(previous_alignment_chrom, BL_SAM_RNAME(alignment), BL_CHROM_MAX_CHARS + 1);
	    else if ( cmp < 0 )
	    {
		fprintf(stderr, "fasda, %s(): "
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
		       bl_alignment_stats_t *alignment_stats, int flags)

{
    int64_t overlapping_reads = 0;
    double  est_counts;
    int     cmp, read_status;
    long    buffer_pos;
    bool    do_not_count;
    
    buffer_pos = ftell(buffer_stream);
    
    /*
     *  "alignment" arg already contains first overlapping read
     *  Loop through all reads overlapping the feature
     */
    
    // Check alignments buffered from previous gene first
    do
    {
	// fprintf(stderr, "Counting buffered alignemnts...\n");
	cmp = chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom, flags);
	if ( cmp > 0 )
	    strlcpy(previous_alignment_chrom, BL_SAM_RNAME(alignment), BL_CHROM_MAX_CHARS + 1);
	else if ( cmp < 0 )
	{
	    fprintf(stderr, "fasda: %s(): "
		    "Chromosomes out of order in SAM stream: %s %s\n",
		    __FUNCTION__, previous_alignment_chrom, BL_SAM_RNAME(alignment));
	    exit(EX_DATAERR);
	}
	
	// FIXME: What's the best combination of alignment flags to
	// include/exclude reads?
	// PROPER_PAIR means both mates mapped together
	// Paired end fragments will be counted about double, but this
	// doesn't matter for differential analysis, where ratios across
	// across conditions will be the same.  If using this for other
	// purposes, we may need to adjust.
	do_not_count = 
	    BL_SAM_FLAG_SECONDARY|BL_SAM_FLAG_QCFAIL|
	    BL_SAM_FLAG_DUP|BL_SAM_FLAG_SUPPLEMENTARY;
	if ( (BL_SAM_FLAG(alignment) & BL_SAM_FLAG_PROPER_PAIR) &&
	     ! (BL_SAM_FLAG(alignment) & do_not_count) )
	    ++overlapping_reads;
	
	// Finding fragment length for paired end requires
	// finding the mate as well.
	// fragment_length = ??;
	// sum_length += fragment_length;
    }   while ( ((read_status = bl_sam_read(alignment, buffer_stream, SAM_MASK))
			    == BL_READ_OK)
		&& (bl_sam_gff_cmp(alignment, feature) == 0) );
    
    // If we ran out of buffered alignments, discard the old ones and
    // continue in primary SAM stream
    if ( read_status == BL_READ_EOF )
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
	    cmp = chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom, flags);
	    if ( cmp > 0 )
		strlcpy(previous_alignment_chrom, BL_SAM_RNAME(alignment), BL_CHROM_MAX_CHARS + 1);
	    else if ( cmp < 0 )
	    {
		fprintf(stderr, "fasda: %s(): "
			"Chromosomes out of order in SAM stream: %s %s\n",
			__FUNCTION__, previous_alignment_chrom, BL_SAM_RNAME(alignment));
		exit(EX_DATAERR);
	    }
	    ++overlapping_reads;
    
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
	    
	    // Finding fragment length for paired end requires
	    // finding the mate as well.
	    // fragment_length = ??;
	    // sum_length += fragment_length;
	}
    }
    
    fseek(buffer_stream, buffer_pos, SEEK_SET);
    
    // http://robpatro.com/blog/?p=235#efflen
    // Count only fragments <= transcript length
    // avg_fragment_length = sum_fragment_length / fragments <= feature len
    // Actually "expected" effective length.  Effective length refers
    // to individual fragments and we want the average of all fragments
    // shorter than the feature.
    // effective_length = feature_length - avg_fragment_length;

    // FIXME: Find out how kallisto computes est_counts
    est_counts = overlapping_reads;
    
    return est_counts;
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
			int64_t length, double est_counts, int flags)

{
    char            *id, *p;
    double          eff_length, tpm;

    if ( flags & FASDA_FLAG_SHOW_GENE )
	id = BL_GFF_FEATURE_NAME(feature);
    else
	id = BL_GFF_FEATURE_ID(feature);
    
    // Drop "gene:" or "transcript:"
    if ( (p = strchr(id, ':')) != NULL )
	id = p + 1;
    
    // Some IDs contain aliases separated by '|'.  The last is usually
    // the unique ID.
    if ( (p = strrchr(id, '|')) != NULL )
	id = p + 1;
    
    // FIXME: http://robpatro.com/blog/?p=235
    // https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
    eff_length = 0.0;
    tpm = 0.0;
    
    fprintf(abundance_stream, "%s\t%" PRId64 "\t%f\t%f\t%f\n",
	    id, length, eff_length, est_counts, tpm);
    return 0;
}


int     chrom_name_cmp(const char *s1, const char *s2, int flags)

{
    // Return "true" if we're supposed to ignore order
    if ( flags & FASDA_IGNORE_CHR_ORDER )
	return 1;
    else
	return bl_chrom_name_cmp(s1, s2);
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s \\\n"
	    "\t[--debug] \\\n"
	    "\t[--show-gene-name] \\\n"
	    "\t[--ignore-chromosome-order] \\\n"
	    "\t[--feature-type mRNA|transcript|gene (default=mRNA)] \\\n"
	    "\t[--output-dir dir (default=same as SAM/BAM/CRAM input)] \\\n"
	    "\tfeatures.gff3 \\\n"
	    "\tfile.[sam|bam|cram]" XT_COMPRESSION_EXTENSIONS " \\\n"
	    "\t[file.[sam|bam|cram]" XT_COMPRESSION_EXTENSIONS " ...]\n", argv[0]);
    exit(EX_USAGE);
}
