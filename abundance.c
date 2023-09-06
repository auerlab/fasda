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
#include <regex.h>          // Pattern match feature types or IDs
#include <sys/types.h>      // off_t
#include <xtend/file.h>
#include <xtend/string.h>   // strlcpy() on Linux
#include <xtend/dsv.h>
#include <biolibc/gff3.h>
#include <biolibc/sam.h>
#include <biolibc/biostring.h>
#include "abundance.h"

bool    Debug = false;

int     main(int argc,char *argv[])

{
    char        *feature_file,
		*sam_files[MAX_FILE_COUNT],
		*abundance_files[MAX_FILE_COUNT],
		*p, *endp,
		*feature_type = "RNA$|transcript$|gene_segment$",
		*output_dir = NULL,
		*gtf_files[MAX_FILE_COUNT];
    int         file_count, c, flags, status;
    unsigned    read_length;
    abundance_method_t  method = STRINGTIE;
    FILE        *abundance_streams[MAX_FILE_COUNT];
    
    if ( argc < 3 )
	usage(argv);

    // FIXME: Add --require-flags, --include-flags, --exclude-flags
    // See samtools-view
    for (c = 1, flags = 0x0; (c < argc) && (argv[c][0] == '-'); ++c)
    {
	if ( strcmp(argv[c], "--debug") == 0 )
	    Debug = true;
	else if ( strcmp(argv[c], "--exact") == 0 )
	    method = EXACT;
	else if ( strcmp(argv[c], "--stringtie") == 0 )
	    flags |= STRINGTIE;
	else if ( strcmp(argv[c], "--show-gene-name") == 0 )
	    flags |= FASDA_FLAG_SHOW_GENE;
	else if ( strcmp(argv[c], "--ignore-chromosome-order") == 0 )
	{
	    flags |= FASDA_IGNORE_CHR_ORDER;
	    fputs("Warning: Output may be incorrect if GFF and SAM input are not in the same order.\n", stderr);
	}
	else if ( strcmp(argv[c], "--feature-type") == 0 )
	    feature_type = argv[++c];
	// FIXME: Add --feature-id to filter by ID in col 9 rather than col 3
	else if ( strcmp(argv[c], "--output-dir") == 0 )
	    output_dir = argv[++c];
	else
	    usage(argv);
    }

    read_length = strtoul(argv[c++], &endp, 10);
    if ( *endp != '\0' )
	usage(argv);
    
    feature_file = argv[c++];
    if ( c == argc )
	usage(argv);

    for (file_count = 0; c < argc; ++c, ++file_count)
    {
	/*
	 *  Open SAM/BAM/CRAM input file
	 */
	
	sam_files[file_count] = argv[c];
	
	/*
	 *  Open abundance output file for this input
	 */
	
	// FIXME: Compute actual length from removing extension and
	// adding "-abundance.tsv".
	abundance_files[file_count] = malloc(PATH_MAX + 1);
	if ( abundance_files[file_count] == NULL )
	{
	    fprintf(stderr, "abundance: Could not malloc() abundance_files[%d]\n",
		    file_count);
	    return EX_UNAVAILABLE;
	}

	// We only need GTFs for stringtie, but build the filenames
	// here along with abundance.tsv files since the process is the same
	gtf_files[file_count] = strdup(sam_files[file_count]);
	if ( gtf_files[file_count] == NULL )
	{
	    fprintf(stderr, "abundance: Could not strdup() gtf_files[%d]\n",
		    file_count);
	    return EX_UNAVAILABLE;
	}
	
	if ( output_dir != NULL )
	{
	    // Get base name of SAM input
	    // FIXME: Use basename()?  Implementation details vary
	    // according to the FreeBSD man page.
	    if ( (p = strrchr(sam_files[file_count], '/')) != NULL )
		++p;
	    else
		p = sam_files[file_count];
	    
	    snprintf(abundance_files[file_count], PATH_MAX + 1, "%s/%s",
		    output_dir, p);
	    snprintf(gtf_files[file_count], PATH_MAX + 1, "%s/%s",
		    output_dir, p);
	}
	else
	{
	    strlcpy(abundance_files[file_count], sam_files[file_count], PATH_MAX + 1);
	    strlcpy(gtf_files[file_count], sam_files[file_count], strlen(gtf_files[file_count]) + 1);
	}
	
	if ( (p = bl_sam_filename_extension(abundance_files[file_count])) == NULL )
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
		"target_id\tlength\teff_length\test_counts\ttpm\tfpkm\n");
	fflush(abundance_streams[file_count]);
	
	if ( (p = bl_sam_filename_extension(gtf_files[file_count])) == NULL )
	{
	    fprintf(stderr, "abundance: Expected sam|bam|cram: %s\n",
		    abundance_files[file_count]);
	    return EX_NOINPUT;
	}
	strlcpy(p, ".gtf", 5);
	if ( method == STRINGTIE )
	    fprintf(stderr, "Writing GTF to %s\n", gtf_files[file_count]);
    }
    
    switch(method)
    {
	case    EXACT:
	    return exact_abundance(feature_file, sam_files, abundance_streams,
			 file_count, feature_type, flags);
	case    STRINGTIE:
	    /*
	     *  stringtie GTFs come out sorted differently, and
	     *  abundance files are sorted the same.  fasda normalize
	     *  needs them in the same sort order.
	     */
	    status = stringtie_abundance(feature_file, sam_files, 
			 gtf_files, abundance_streams,
			 file_count, feature_type, read_length, flags);
	    sort_abundance(abundance_files, file_count);
	    return status;
	
	default:
	    fprintf(stderr, "Invalid abundance method: %d\n"
		    "This is a software bug.\n", method);
	    return EX_SOFTWARE;
    }
}


int     exact_abundance(const char *feature_file, char *sam_files[],
		 FILE *abundance_streams[], int file_count,
		 const char *feature_type, int flags)

{
    char        previous_feature_chrom[BL_CHROM_MAX_CHARS + 1],
		previous_alignment_chrom[MAX_FILE_COUNT][BL_CHROM_MAX_CHARS + 1],
		*transcript_id, *p;
    bl_gff3_t    feature, subfeature;
    bl_sam_t    alignment;
    int         c, cmp, status;
    int64_t     length, previous_start;
    double      est_counts = 0.0;
    bool        alternate_exons;
    FILE        *feature_stream,
		*sam_streams[MAX_FILE_COUNT],
		*buffer_streams[MAX_FILE_COUNT];
    unsigned long   features_processed;
    bl_alignment_stats_t    alignment_stats = BL_ALIGNMENT_STATS_INIT;
    regex_t     feature_re;
    
    fputs("Warning: Exact abundance calculation is incomplete.\n"
	  "Use the default stringtie abundance calculation for more robust results.\n",
	  stderr);
    
    if ( (feature_stream = xt_fopen(feature_file, "r")) == NULL )
    {
	fprintf(stderr, "fasda: Could not open %s for read: %s.\n",
		feature_file, strerror(errno));
	return EX_NOINPUT;
    }

    for (c = 0; c < file_count; ++c)
    {
	if ( (sam_streams[c] =
	      bl_sam_fopen(sam_files[c], "r",
	      SAMTOOLS_ARGS)) == NULL )
	{
	    fprintf(stderr, "abundance: Could not open %s for read: %s.\n",
		    sam_files[file_count], strerror(errno));
	    return EX_NOINPUT;
	}
	bl_sam_skip_header(sam_streams[c]);
    }
    
    bl_gff3_skip_header(feature_stream);
    
    bl_gff3_init(&feature);
    bl_gff3_init(&subfeature);
    bl_sam_init(&alignment);

    regcomp(&feature_re, feature_type, REG_EXTENDED);
    
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
    
    features_processed = 0;
    while ( bl_gff3_read(&feature, feature_stream, GFF3_MASK) == BL_READ_OK )
    {
	// Some IDs contain aliases separated by '|'.  The last is usually
	// the unique ID.
	transcript_id = BL_GFF3_FEATURE_ID(&feature);
	if ( (transcript_id != NULL) && (p = strrchr(transcript_id, '|')) != NULL )
	    transcript_id = p + 1;

	if ( false ) // Debug )
	    fprintf(stderr, "%s %s %s\n",
		    transcript_id, BL_GFF3_TYPE(&feature), feature_type);
	// if ( strcmp(BL_GFF3_TYPE(&feature), feature_type) == 0 )
	if ( regexec(&feature_re, BL_GFF3_TYPE(&feature), 0, NULL, 0) == 0 )
	{
	    // Verify that features are properly sorted
	    cmp = fasda_chrom_name_cmp(BL_GFF3_SEQID(&feature),
				    previous_feature_chrom, flags);
	    if ( cmp > 0 )
		strlcpy(previous_feature_chrom, BL_GFF3_SEQID(&feature),
			BL_CHROM_MAX_CHARS + 1);
	    else if ( cmp < 0 )
	    {
		fprintf(stderr, "fasda: Error: GFF3 chromosomes out of order: %s %s\n",
			previous_feature_chrom, BL_GFF3_SEQID(&feature));
		return EX_DATAERR;
	    }

	    // Compute sum of exon lengths for abundance TSV output
	    length = 0;
	    previous_start = 0;
	    alternate_exons = false;
	    
	    if ( Debug )
		fprintf(stderr, "Finding exons of %s...\n", transcript_id);
	    while ( ((status = bl_gff3_read(&subfeature, feature_stream, GFF3_MASK))
			== BL_READ_OK)
		    && (strcmp(BL_GFF3_TYPE(&subfeature), "###") != 0)
		    && (regexec(&feature_re, BL_GFF3_TYPE(&subfeature),
			0, NULL, 0) != 0) )
	    {
		// Stop counting exons as soon as one has a lower start
		// position than the previous one.  This means we're into
		// alternate splicing data.  Lengths computed this way match
		// kallisto's abundance.tsv.
		if ( ! alternate_exons && 
		     (strcmp(BL_GFF3_TYPE(&subfeature), "exon") == 0) )
		{
		    //fprintf(stderr, "%ld %ld\n",
		    //        previous_start, BL_GFF3_START(&subfeature));
		    if ( BL_GFF3_START(&subfeature) > previous_start )
		    {
			if ( Debug )
			    fprintf(stderr, "Adding exon %" PRId64 " %" PRId64 " ",
				    BL_GFF3_START(&subfeature), BL_GFF3_END(&subfeature));
			length += BL_GFF3_END(&subfeature) -
			    BL_GFF3_START(&subfeature) + 1;
			if ( Debug )
			    fprintf(stderr, "length = %" PRId64 "\n", length);
			previous_start = BL_GFF3_START(&subfeature);
		    }
		    else
			alternate_exons = true;
		}
	    }
	    
	    if ( Debug )
	    {
		fprintf(stderr, "Found %s  length = %" PRId64 " status = %d\n",
			BL_GFF3_TYPE(&feature), length, status);
	    }
	    
	    // Last feature read was not a subfeature, so rewind to
	    // it so outer loop reads it again as a primary feature
	    // FIXME: Maybe don't use BL_GFF3_FILE_POS()
	    if ( strcmp(BL_GFF3_TYPE(&subfeature), "###") != 0 )
	    {
		if ( Debug )
		    fprintf(stderr, "Rewinding to %zu (%s)...\n",
			    BL_GFF3_FILE_POS(&subfeature),
			    BL_GFF3_TYPE(&subfeature));
		fseeko(feature_stream, BL_GFF3_FILE_POS(&subfeature), SEEK_SET);
	    }
	    
	    for (c = 0; c < file_count; ++c)
	    {
		/*
		 *  Find first overlapping read for this gene.
		 */

		if ( (bl_gff3_find_overlapping_alignment(&feature,
			sam_streams[c], buffer_streams[c],
			previous_alignment_chrom[c],
			&alignment, &alignment_stats, flags) == BL_READ_OK) &&
		     (bl_sam_gff3_cmp(&alignment, &feature) == 0) )
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
	    
	    ++features_processed;
	    if ( isatty(fileno(stdout)) && (features_processed % 100 == 0) )
	    {
		printf("Chrom: %3s  Features processed: %lu\r",
			BL_GFF3_SEQID(&feature), features_processed);
		fflush(stdout);
	    }
	}
    }
    
    for (c = 0; c < file_count; ++c)
    {
	bl_sam_fclose(sam_streams[c]);
	fclose(buffer_streams[c]);
	fclose(abundance_streams[c]);
    }
    xt_fclose(feature_stream);
    
    printf("\nTotal features processed:          %lu\n", features_processed);
    printf("Total alignments processed:        %lu\n",
	    BL_ALIGNMENT_STATS_TOTAL(&alignment_stats));
    printf("Alignments overlapping a feature:  %lu\n",
	    BL_ALIGNMENT_STATS_OVERLAPPING(&alignment_stats));
    
    return EX_OK;
}


unsigned long   sum_exons(FILE *gtf_stream)

{
    xt_dsv_line_t  *dsv_line = xt_dsv_line_new();
    int         delim;
    off_t       pos;
    bool        is_exon;
    char        *gtf_feature_type,
		*gtf_start_text,
		*gtf_end_text,
		*endp;
    unsigned long   start, end, length;
    
    length = 0;
    do 
    {
	pos = ftello(gtf_stream);
	if ( (delim = xt_dsv_line_read(dsv_line, gtf_stream, "\t")) != EOF )
	{
	    gtf_feature_type = xt_dsv_line_get_fields_ae(dsv_line, 2);
	    // fprintf(stderr, "Feature = %s\n", gtf_feature_type);
	    is_exon = (strcmp(gtf_feature_type, "exon") == 0);
	    if ( is_exon )
	    {
		gtf_start_text = xt_dsv_line_get_fields_ae(dsv_line, 3);
		gtf_end_text = xt_dsv_line_get_fields_ae(dsv_line, 4);
		
		start = strtol(gtf_start_text, &endp, 10);
		if ( *endp != '\0' )
		{
		    fprintf(stderr, "Invalid start locus: %s\n", gtf_start_text);
		    exit(EX_DATAERR);
		}
		end = strtol(gtf_end_text, &endp, 10);
		if ( *endp != '\0' )
		{
		    fprintf(stderr, "Invalid end locus: %s\n", gtf_end_text);
		    exit(EX_DATAERR);
		}
		
		// fprintf(stderr, "Adding exon length %lu\n", end - start + 1);
		length += end - start + 1;
	    }
	}
    }   while ( (delim != EOF) && is_exon );
    fseeko(gtf_stream, pos, SEEK_SET);
    
    return length;
}


/***************************************************************************
 *  Description:
 *      Spawn a stringtie process to compute abundances and reformat
 *      the output to a kallisto-style abundance.tsv file.
 *  
 *  History: 
 *  Date        Name        Modification
 *  2023-07-17  Jason Bacon Begin
 ***************************************************************************/

int     stringtie_abundance(const char *feature_file, char *sam_files[],
		 char *gtf_files[],
		 FILE *abundance_streams[], int file_count,
		 const char *feature_type, unsigned read_length, int flags)

{
    FILE    *gtf_stream;
    char    cmd[MAX_CMD_LEN + 1],
	    gene_abundance[PATH_MAX + 1],
	    *gtf_feature_type,
	    *gtf_attributes,
	    *field,
	    *p, *p2, *endp,
	    *transcript_id,
	    *coverage,
	    *tpm, *fpkm;
    unsigned long   length;
    double  eff_length, cov, est_count;
    int     c, status;
    xt_dsv_line_t  *dsv_line = xt_dsv_line_new();
    bool    field_ok;
    
    for (c = 0; c < file_count; ++c)
    {
	
	// Note: BAMs must be properly sorted (samtools sort default order)
	// Write stringtie GTF to a file so we can use seek() while reading
	// It's also just nice to have the GTFs to look at
	strlcpy(gene_abundance, gtf_files[c], PATH_MAX + 1);
	if ( (p = strstr(gene_abundance, ".gtf")) == NULL )
	{
	    fprintf(stderr, "stringtie_abundance(): GTF filename does not contain \".gtf\"\n");
	    exit(EX_SOFTWARE);
	}
	strlcpy(p, "-stringtie-genes.tsv", PATH_MAX + 1);
	snprintf(cmd, MAX_CMD_LEN + 1, "stringtie -e -G %s -o %s -A %s %s",
		 feature_file, gtf_files[c], gene_abundance, sam_files[c]);
	
	fprintf(stderr, "Running %s...\n", cmd);
	if ( (status = system(cmd)) != 0 )
	{
	    fprintf(stderr, "abundance: Error: %s\n", cmd);
	    fprintf(stderr, "status = %d\n", status);
	    return EX_SOFTWARE;
	}
	
	if ( (gtf_stream = fopen(gtf_files[c], "r")) == NULL )
	{
	    fprintf(stderr, "abundance: Could not open %s.\n", gtf_files[c]);
	    return EX_NOINPUT;
	}
	
	// FIXME: Create GTF class in biolibc?
	while ( xt_dsv_line_read(dsv_line, gtf_stream, "\t") != EOF )
	{
	    // Skip header and comment lines
	    if ( *xt_dsv_line_get_fields_ae(dsv_line, 0) != '#' )
	    {
		gtf_feature_type = xt_dsv_line_get_fields_ae(dsv_line, 2);
		gtf_attributes = xt_dsv_line_get_fields_ae(dsv_line, 8);
		if ( strcmp(gtf_feature_type, "transcript") == 0 )
		{
		    if ( Debug )
			fprintf(stderr, "Seq: %s  Feature type: %s\n",
			    xt_dsv_line_get_fields_ae(dsv_line, 0),
			    gtf_feature_type);
		    //printf("%s\n", gtf_start_text);
		    //printf("%s\n", gtf_end_text);
		    //printf("Attributes: %s\n", gtf_attributes);
		    // dsv_line_write(dsv_line, stdout);
		    length = sum_exons(gtf_stream);
		    
		    p = gtf_attributes;
		    while ( (field = strsep(&p, ";")) != NULL )
		    {
			// puts(field);
			field_ok = false;
			
			// FIXME: Factor out attribute extraction function
			// Better yet, implement a GTF class in biolibc
			// If transcript_id comes after a ;, it will start
			// with a space
			if ( (memcmp(field, " transcript_id ", 15) == 0) ||
			     (memcmp(field, "transcript_id ", 14) == 0) )
			{
			    if ( *field == ' ' )
				transcript_id = field + 15;
			    else
				transcript_id = field + 14;
			    
			    if ( *transcript_id == '"' )
			    {
				++transcript_id;
				if ( memcmp(transcript_id, "transcript:", 11) == 0 )
				    transcript_id += 11;

				// Some IDs contain aliases separated by '|'.
				// The last is usually the unique ID.
				if ( (p2 = strrchr(transcript_id, '|')) != NULL )
				    transcript_id = p2 + 1;
				// fprintf(stderr, "Clipped alias: %s\n", transcript_id);
			    
				if ( (p2 = strchr(transcript_id, '"')) != NULL )
				    *p2 = '\0';
				field_ok = true;
			    }
			    
			    if ( ! field_ok )
			    {
				fprintf(stderr, "Malformed transcript_id field: %s.  Expected '\"'.\n",
					field);
				exit(EX_DATAERR);
			    }
			    
			    if ( Debug )
				fprintf(stderr, "transcript_id = %s\n", transcript_id);
			}
			else if ( memcmp(field, " cov ", 5) == 0 )
			{
			    p2 = field + 5;
			    if ( *p2 == '"' )
			    {
				coverage = ++p2;
				if ( (p2 = strchr(coverage, '"')) != NULL )
				    *p2 = '\0';
				field_ok = true;
			    }
			    
			    if ( ! field_ok )
			    {
				fprintf(stderr, "Malformed coverage field: %s.  Expected '\"'.\n",
					field);
				exit(EX_DATAERR);
			    }
			    
			    // printf("coverage = %s\n", coverage);
			}
			else if ( memcmp(field, " TPM ", 5) == 0 )
			{
			    p2 = field + 5;
			    if ( *p2 == '"' )
			    {
				tpm = ++p2;
				if ( (p2 = strchr(tpm, '"')) != NULL )
				    *p2 = '\0';
				field_ok = true;
			    }
			    
			    if ( ! field_ok )
			    {
				fprintf(stderr, "Malformed TPM field: %s.  Expected '\"'.\n",
					field);
				exit(EX_DATAERR);
			    }
			    
			    // printf("TPM = %s\n", tpm);
			}
			else if ( memcmp(field, " FPKM ", 6) == 0 )
			{
			    p2 = field + 6;
			    if ( *p2 == '"' )
			    {
				fpkm = ++p2;
				if ( (p2 = strchr(fpkm, '"')) != NULL )
				    *p2 = '\0';
				field_ok = true;
			    }
			    
			    if ( ! field_ok )
			    {
				fprintf(stderr, "Malformed FPKM field: %s.  Expected '\"'.\n",
					field);
				exit(EX_DATAERR);
			    }
			    
			    // printf("TPM = %s\n", tpm);
			}
		    }

		    // FIXME: Replace coverage with count
		    // From http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
		    // regarding prepDE.py3
		    // reads_per_transcript = coverage * transcript_len / read_len
		    cov = strtod(coverage, &endp);
		    if ( *endp != '\0' )
		    {
			fprintf(stderr, "Invalid coverage: %s\n", coverage);
			exit(EX_DATAERR);
		    }
		    est_count = (double)cov * length / read_length;
		    // FIXME: Explore how kallisto computes this
		    // For now, flag it as invalid
		    eff_length = 0;
		    //fprintf(stderr, "Outputting %s...\n", transcript_id);
		    fprintf(abundance_streams[c],
			    "%s\t%lu\t%0.1f\t%0.1f\t%s\t%s\n",
			    transcript_id, length, eff_length, est_count, tpm, fpkm);
		}
	    }
	    xt_dsv_line_free(dsv_line);
	}
	
	fclose(gtf_stream);
	fclose(abundance_streams[c]);
    }
    
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

int     bl_gff3_find_overlapping_alignment(bl_gff3_t *feature,
	    FILE *sam_stream, FILE *buffer_stream,
	    char *previous_alignment_chrom, bl_sam_t *alignment,
	    bl_alignment_stats_t *alignment_stats, int flags)

{
    int     status, cmp;

    // First search alignments buffered from previous gene
    // fprintf(stderr, "Checking buffered alignments...\n");
    while ( ((status = bl_sam_read(alignment, buffer_stream, SAM_MASK)) == BL_READ_OK) &&
	    (bl_sam_gff3_cmp(alignment, feature) < 0) )
    {
	// fprintf(stderr, "buffered: %s %lu\n", BL_SAM_RNAME(alignment), BL_SAM_POS(alignment));
	// Verify that alignments are properly sorted
	cmp = fasda_chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom, flags);
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
		(bl_sam_gff3_cmp(alignment, feature) < 0) )
	{
	    ++BL_ALIGNMENT_STATS_TOTAL(alignment_stats);
	    
	    // Verify that alignments are properly sorted
	    cmp = fasda_chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom, flags);
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

double  count_coverage(bl_gff3_t *feature, bl_sam_t *alignment,
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
	cmp = fasda_chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom, flags);
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
		&& (bl_sam_gff3_cmp(alignment, feature) == 0) );
    
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
		&& (bl_sam_gff3_cmp(alignment, feature) == 0) )
	{
	    ++BL_ALIGNMENT_STATS_TOTAL(alignment_stats);
	    ++BL_ALIGNMENT_STATS_OVERLAPPING(alignment_stats);
	    cmp = fasda_chrom_name_cmp(BL_SAM_RNAME(alignment), previous_alignment_chrom, flags);
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

int     print_abundance(FILE *abundance_stream, bl_gff3_t *feature,
			int64_t length, double est_counts, int flags)

{
    char            *id, *p;
    double          eff_length, tpm;

    if ( flags & FASDA_FLAG_SHOW_GENE )
	id = BL_GFF3_FEATURE_NAME(feature);
    else
	id = BL_GFF3_FEATURE_ID(feature);
    
    // Drop "gene:" or "transcript:"
    if ( (p = strchr(id, ':')) != NULL )
	id = p + 1;
    
    // Some IDs contain aliases separated by '|'.  The last is usually
    // the unique ID.
    if ( (p = strrchr(id, '|')) != NULL )
	id = p + 1;
    
    // FIXME: http://robpatro.com/blog/?p=235
    // https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
    // Explore how kallisto computes eff_length
    eff_length = length;
    tpm = 0.0;
    
    // FIXME: Add FPKM in last field?
    fprintf(abundance_stream, "%s\t%" PRId64 "\t%0.1f\t%0.1f\t%0.2f\t*\n",
	    id, length, eff_length, est_counts, tpm);
    return 0;
}


int     fasda_chrom_name_cmp(const char *s1, const char *s2, int flags)

{
    // Return "true" if we're supposed to ignore order
    if ( flags & FASDA_IGNORE_CHR_ORDER )
	return 1;
    else
	return bl_chrom_name_cmp(s1, s2);
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
 *  2023-07-18  Jason Bacon Begin
 ***************************************************************************/

char    *bl_sam_filename_extension(const char *filename)

{
    char  *p;
    
    // Replace .sam/.bam/.cram with -abundance.tsv
    if ( (p = strstr(filename, ".bam")) == NULL )
	if ( (p = strstr(filename, ".sam")) == NULL )
	    p = strstr(filename, ".cram");
    
    return p;
}


void    sort_abundance(char *abundance_files[], int file_count)

{
    int     c;
    char    cmd[MAX_CMD_LEN + 1];
    
    // FIXME: Check system() return vals
    for (c = 0; c < file_count; ++c)
    {
	snprintf(cmd, MAX_CMD_LEN + 1, "head -n 1 %s > %s.sorted",
		 abundance_files[c], abundance_files[c]);
	system(cmd);
	snprintf(cmd, MAX_CMD_LEN + 1, "fgrep -v eff_length %s | sort >> %s.sorted",
		 abundance_files[c], abundance_files[c]);
	system(cmd);
	snprintf(cmd, MAX_CMD_LEN + 1, "mv -f %s.sorted %s",
		 abundance_files[c], abundance_files[c]);
	system(cmd);
    }
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s \\\n"
	    "\t[--debug] \\\n"
	    "\t[--stringtie] \\\n"
	    "\t[--exact] \\\n"
	    "\t[--show-gene-name] \\\n"
	    "\t[--ignore-chromosome-order] \\\n"
	    "\t[--feature-type mRNA|transcript|gene (default=mRNA)] \\\n"
	    "\t[--output-dir dir (default=same as SAM/BAM/CRAM input)] \\\n"
	    "\tread-length \\\n"
	    "\tfeatures.gff3 (gtf also works with --stringtie)\\\n"
	    "\tfile.[sam|bam|cram]" XT_COMPRESSION_EXTENSIONS " \\\n"
	    "\t[file.[sam|bam|cram]" XT_COMPRESSION_EXTENSIONS " ...]\n", argv[0]);
    exit(EX_USAGE);
}
