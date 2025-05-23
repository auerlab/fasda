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
#include <math.h>
#include <limits.h>         // PATH_MAX OpenIndiana
#include <stdbool.h>
#include <sys/param.h>      // PATH_MAX
#include <xtend/dsv.h>
#include <xtend/file.h>
#include <xtend/mem.h>
#include <xtend/math.h>     // xt_double_cmp()
#include <xtend/string.h>   // strlcpy() on Linux
#include "normalize.h"

bool    Debug = false;

int     main(int argc, const char *argv[])

{
    int     arg;
    FILE    *norm_all_stream = stdin;
    
    if ( argc < 3 )
	usage(argv);

    for (arg = 1; *argv[arg] == '-'; ++arg)
    {
	// Only median ratios is supported for now, so ignore
	if ( strcmp(argv[arg], "--mrn") == 0 )
	    ;
	
	if ( strcmp(argv[arg], "--debug") == 0 )
	    Debug = true;
	
	else if ( strcmp(argv[arg], "--output") == 0 )
	{
	    if ( (norm_all_stream = fopen(argv[++arg], "w")) == NULL )
	    {
		fprintf(stderr, "normalize: Could not open %s for write: %s.\n",
			argv[arg], strerror(errno));
		return EX_CANTCREAT;
	    }
	}
	else
	    usage(argv);
    }
    
    // Remaining arguments are input files
    if ( argc - arg < 2 )
    {
	fprintf(stderr, "%s requires a minimum of two abundance files.\n",
		argv[0]);
	usage(argv);
    }
    return mrn(argv + arg, norm_all_stream);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *      Perform median ratio normalization (MRN) on a set of abundance
 *      files, writing output to norm_all_stream.  MRN involves
 *      computing the average abundance of all replicates for each
 *      feature (the pseudo-reference), dividing each count for the
 *      feature by this reference, and finally taking the median of
 *      these ratios for each replicate as the scaling factor.
 *      Computations are done on log(count) values here for logistical
 *      reasons, subtracting log(counts) to simulate dividing counts.
 *      This does not produce identical results, but does produce
 *      consistent ratios.
 *
 *      abundance_files is an argv-style pointer array terminated by
 *      a NULL pointer.  Each file should follow the format of
 *      abundance.tsv output by kallisto.
 *
 *      Output written to norm_all_stream is a TSV file containing the
 *      feature name in column 1 followed by a normalized count for each
 *      replicate (input file) in subsequent columns.
 *  
 *  Arguments:
 *      abundance_files     argv-style array of abundance file names
 *      norm_all_output     FILE pointer to TSV output stream
 *
 *  Returns:
 *      Exit status appropriate for passing back from main()
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-05-14  Jason Bacon Begin
 ***************************************************************************/

int     mrn(const char *abundance_files[], FILE *norm_all_stream)

{
    xt_dsv_line_t  *dsv_lines[FASDA_MAX_SAMPLES];
    size_t      sample, sample_count, feature_count, non_zero_feature_count, c;
    FILE        *abundance_streams[FASDA_MAX_SAMPLES],
		*tmp_streams[FASDA_MAX_SAMPLES],
		*norm_sample_streams[FASDA_MAX_SAMPLES];
    char        *end, *target_id,
		norm_sample_file[PATH_MAX + 1], *p;
    double      est_count, sum_log_count, log_count[FASDA_MAX_SAMPLES],
		mean_log_count, *ratios, median_ratio[FASDA_MAX_SAMPLES],
		scaling_factor[FASDA_MAX_SAMPLES];
    
    for (sample_count = 0; abundance_files[sample_count] != NULL; ++sample_count)
    {
	if ( (abundance_streams[sample_count] =
		xt_fopen(abundance_files[sample_count], "r")) == NULL )
	{
	    fprintf(stderr, "normalize: Could not open %s for read: %s.\n",
		    abundance_files[sample_count], strerror(errno));
	    return EX_NOINPUT;
	}
	
	if ( (tmp_streams[sample_count] = tmpfile()) == NULL )
	{
	    fprintf(stderr, "normalize: Could not open temp file: %s\n",
		    strerror(errno));
	    return EX_NOINPUT;
	}
	dsv_lines[sample_count] = xt_dsv_line_new();
    }
    
    skip_headers(abundance_files, abundance_streams, dsv_lines, sample_count);

    // Abundance file format:
    // target_id       length  eff_length      est_counts      tpm
    
    feature_count = non_zero_feature_count = 0;
    while ( ! feof(abundance_streams[0]) )
    {
	for (sample = 0, sum_log_count = 0; sample < sample_count; ++sample)
	{
	    // Read next entry from every abundance file
	    if ( xt_dsv_line_read(dsv_lines[sample], abundance_streams[sample],
				"\t") == EOF )
	    {
		// Check for EOF on every abundance file.  They should all
		// come in the same iteration.
		check_all_eof(abundance_files, abundance_streams, sample, sample_count);
		break;  // Back to outer while loop
	    }
	    else
	    {
		// Feature name
		target_id = xt_dsv_line_get_fields_ae(dsv_lines[sample], 0);
		
		// Corresponding lines in each abundance file should
		// contain the same feature
		// FIXME: Skip header line before loop rather than
		// check sample > 0 in every iteration
		if ( (sample > 0) && (strcmp(target_id,
			xt_dsv_line_get_fields_ae(dsv_lines[sample - 1], 0)) != 0) )
		{
		    fprintf(stderr,
			    "normalize: %s, %s: Different feature IDs on line %zu\n.\n",
			    abundance_files[sample - 1],
			    abundance_files[sample], feature_count + 1);
		    close_all_streams(abundance_streams, norm_sample_streams,
				      norm_all_stream, sample_count);
		    return EX_DATAERR;
		}

		/*
		 *  Median of ratios normalization
		 *  https://scienceparkstudygroup.github.io/research-data-management-lesson/median_of_ratios_manual_normalization/index.html
		 *  Similar to TMM but more robust: doi 10.1093/bib/bbx008
		 *
		 *  First pass:
		 *
		 *  Read raw counts for all genes and all samples
		 *
		 *  1.  Take log of every count (just for filtering in step 3?)
		 *  2.  Mean of all log(counts) for feature (pseudo-reference)
		 */
		
		// Field 3 is est_counts
		est_count =
		    strtof(xt_dsv_line_get_fields_ae(dsv_lines[sample], 3),
		    &end);
		if ( *end != '\0' )
		{
		    fprintf(stderr, "normalize: Invalid count: %s\n",
			    xt_dsv_line_get_fields_ae(dsv_lines[sample_count],0));
		    return EX_DATAERR;
		}
		log_count[sample] = log(est_count);
		/*
		if ( Debug )
		    fprintf(stderr, "%0.1f ", log_count[sample]);
		*/
		sum_log_count += log_count[sample];
	    }
	}

	/*
	 *  3.  Remove genes with -inf as pseudo-reference
	 *      (mean of log(counts) for feature)
	 *  4.  Subtract pseudo-reference from each log(expression)
	 *      This is actually computing a ratio since subtracting from
	 *      log(v) is dividing v. We'll need to store this value and
	 *      later sort to find median.
	 */   
	
	// FIXME: Is this check redundant?
	if ( ! feof(abundance_streams[0]) )
	{
	    mean_log_count = sum_log_count / sample_count;
	    /*
	    if ( Debug )
		fprintf(stderr, "feature mean [avg log(count)] = %f\n",
			feature_mean);
	    */
	    
	    // Filter out genes with -inf mean
	    // FIXME: This check does not seem to work on all CPUs
	    if ( mean_log_count != -INFINITY )
	    {
		for (sample = 0, sum_log_count = 0; sample < sample_count; ++sample)
		    fprintf(tmp_streams[sample], "%f\n",
			    log_count[sample] - mean_log_count);
		// FIXME: This was outside the if, which caused
		// problems reading back the tmp files
		// How was this working on some systems?? (barracuda)
		++non_zero_feature_count;
	    }
	    ++feature_count;
	}
    }
    
    if ( Debug )
    {
	fprintf(stderr, "%s(): After removing features with mean log count = -INFINITY:\n",
		__FUNCTION__);
	fprintf(stderr, "%s(): Info: feature_count = %zu\n",
		__FUNCTION__, feature_count);
	fprintf(stderr, "%s(): Info: non_zero_feature_count = %zu\n",
		__FUNCTION__, non_zero_feature_count);
    }
    
    /*
     *  Second pass:
     *
     *  5.  Find the median of the log ratios for each sample
     *  6.  exp(median) = scaling factor for the sample
     */

    // Reuse same ratios array for all samples
    ratios = xt_malloc(feature_count, sizeof(*ratios));
    if ( ratios == NULL )
    {
	fprintf(stderr, "normalize: Could not allocate ratios[%zu]\n",
		feature_count);
	return EX_UNAVAILABLE;
    }
    
    for (sample = 0; sample < sample_count; ++sample)
    {
	rewind(tmp_streams[sample]);
	
	for (c = 0; c < non_zero_feature_count; ++c)
	{
	    if ( fscanf(tmp_streams[sample], "%lf", &ratios[c]) != 1 )
	    {
		fprintf(stderr, "%s(): fscanf() failed c = %zu in tmp_streams[%zu]\n.\n",
			__FUNCTION__, c, sample);
		exit(EX_DATAERR);
	    }
	}
	
	qsort(ratios, feature_count, sizeof(*ratios),
	      (int (*)(const void *,const void *))xt_double_cmp);

	if ( Debug )
	{
	    // fprintf(stderr, "Sorted %zu ratios:\n", feature_count);
	    for (c = 0; c < feature_count; ++c)
	    {
		if ( c % 10000 == 0 )
		    fprintf(stderr, "%7zu %f\n", c, ratios[c]);
	    }
	}
	
	// If odd count, media is middle, otherwise mean of two middles
	if ( feature_count % 2 == 1 )
	    median_ratio[sample] = ratios[feature_count / 2];
	else
	    median_ratio[sample] = (ratios[feature_count / 2] +
			    ratios[feature_count / 2 + 1]) / 2;
	scaling_factor[sample] = exp(median_ratio[sample]);
	
	if ( Debug )
	{
	    fprintf(stderr, "Median ratio = %f\t", median_ratio[sample]);
	    fprintf(stderr, "Scaling factor[%zu] = %f\n",
		    sample + 1, scaling_factor[sample]);
	}
    }
    
    /*
     *  Third pass:
     *
     *  7.  Divide counts by scaling factor to normalize
     */
    
    // Rewind input abundance.tsv files and create matching output files
    for (sample = 0; sample < sample_count; ++sample)
    {
	rewind(abundance_streams[sample]);
	// FIXME: libxtend strreplace()?
	strlcpy(norm_sample_file, abundance_files[sample], PATH_MAX + 1);
	if ( (p = strrchr(norm_sample_file, '/')) == NULL )
	    p = norm_sample_file;
	else
	    ++p;    // First after '/'
	*p = '\0';
	strlcat(norm_sample_file, "normalized.tsv", PATH_MAX + 1);
	//fprintf(stderr, "sample = %zu  file = %s\n", sample, norm_sample_file);
	if ( (norm_sample_streams[sample] = xt_fopen(norm_sample_file, "w")) == NULL )
	{
	    fprintf(stderr, "normalize: Could not open %s for write: %s.\n",
		    norm_sample_file, strerror(errno));
	    exit(EX_CANTCREAT);
	}
	// Add header to each output file
	fprintf(norm_sample_streams[sample], "target_id\tlength\teff_length\test_counts\ttpm\tnorm_counts\n");
    }
    
    skip_headers(abundance_files, abundance_streams, dsv_lines, sample_count);
    feature_count = 0;
    while ( ! feof(abundance_streams[0]) )
    {
	for (sample = 0, sum_log_count = 0; sample < sample_count; ++sample)
	{
	    if ( xt_dsv_line_read(dsv_lines[sample], abundance_streams[sample],
				"\t") == EOF )
	    {
		check_all_eof(abundance_files, abundance_streams, sample, sample_count);
		break;
	    }
	    else
	    {
		/*
		 *  Add normalized count to original abundance.tsv for
		 *  viewing alongside raw count and TPM.
		 */
		
		// Copy abundance.tsv fields
		for (c = 0; c < xt_dsv_line_get_num_fields(dsv_lines[sample]); ++c)
		    fprintf(norm_sample_streams[sample], "%s\t",
			    xt_dsv_line_get_fields_ae(dsv_lines[sample], c));
		est_count = strtof(xt_dsv_line_get_fields_ae(dsv_lines[sample], 3), &end);
		if ( *end != '\0' )
		{
		    fprintf(stderr, "normalize: Invalid count: %s\n",
			    xt_dsv_line_get_fields_ae(dsv_lines[sample],0));
		    return EX_DATAERR;
		}
		
		// Add normalized count
		// Divide by scaling factor, don't multiply
		fprintf(norm_sample_streams[sample], "%f\n",
			est_count / scaling_factor[sample]);
		//fprintf(stderr, "sample %zu  nc = %f\n",
		//        sample, count * scaling_factor[sample]);
		//getchar();
		
		/*
		 *  Put all normalized counts together in one file for
		 *  easy reading by fold-change.
		 */

		// Feature ID (target_id) just once
		if ( sample == 0 )
		    fprintf(norm_all_stream, "%s",
			    xt_dsv_line_get_fields_ae(dsv_lines[0],0));
		
		fprintf(norm_all_stream, "\t%f",
			est_count / scaling_factor[sample]);
	    }
	}
	
	if ( ! feof(abundance_streams[0]) )
	{
	    putc('\n', norm_all_stream);
	    ++feature_count;
	}
    }

    close_all_streams(abundance_streams, norm_sample_streams,
		      norm_all_stream, sample_count);
    return EX_OK;
}




/***************************************************************************
 *  Description:
 *      Close all input and output streams
 *  
 *  History: 
 *  Date        Name        Modification
 *  2022-05-19  Jason Bacon Begin
 ***************************************************************************/

void    close_all_streams(FILE *abundance_streams[], FILE *norm_sample_streams[],
			  FILE *norm_all_stream, size_t sample_count)

{
    size_t  sample;
    
    for (sample = 0; sample < sample_count; ++sample)
    {
	fclose(abundance_streams[sample]);
	fclose(norm_sample_streams[sample]);
    }
    fclose(norm_all_stream);
}


/***************************************************************************
 *  Description:
 *      Skip header lines if present in all abundance files
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-05-14  Jason Bacon Begin
 ***************************************************************************/

void    skip_headers(const char *abundance_files[], FILE *abundance_streams[],
		     xt_dsv_line_t *dsv_lines[], size_t sample_count)

{
    size_t  c;
    
    for (c = 0; c < sample_count; ++c)
    {
	// Every abundance file should have a 1-line header
	xt_dsv_line_read(dsv_lines[c], abundance_streams[c], "\t");
	//puts(xt_dsv_line_get_fields_ae(dsv_lines[sample_count], 0));
	if ( strcmp(xt_dsv_line_get_fields_ae(dsv_lines[c], 0), "target_id") != 0 )
	{
	    fprintf(stderr, "normalize: %s: Expected header starting with \"target_id\".\n",
		    abundance_files[c]);
	    fprintf(stderr, "Got %s\n", xt_dsv_line_get_fields_ae(dsv_lines[c], 0));
	    exit(EX_DATAERR);
	}
    }
}


/***************************************************************************
 *  Description:
 *      Make sure all files have reached EOF at the same time.  All files
 *      should have exactly the same set of features.
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-05-14  Jason Bacon Begin
 ***************************************************************************/

void    check_all_eof(const char *abundance_files[], FILE *abundance_streams[],
	      size_t sample, size_t sample_count)

{
    size_t  c;
    
    for (c = 0; c < sample_count; ++c)
	if ( getc(abundance_streams[c]) != EOF )
	{
	    fprintf(stderr, "normalize: EOF reached on %s but not %s\n",
		    abundance_files[sample], abundance_files[c]);
	    exit(EX_DATAERR);
	}
}


void    usage(const char *argv[])

{
    fprintf(stderr, "Usage: %s [--mrn] [--output file.tsv] \\\n"
		    "       abundance1.tsv abundance2.tsv ...\n", argv[0]);
    exit(EX_USAGE);
}
