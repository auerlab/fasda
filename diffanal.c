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
#include <biolibc/gff.h>
#include <xtend/file.h>

#define MAX_CONDITIONS  128

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    char        *features_file, *condition_files[MAX_CONDITIONS];
    FILE        *features_stream, *condition_streams[MAX_CONDITIONS];
    bl_gff_t    feature;
    int         mask, conditions, c;
    
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

	// FIXME: Create bl_sam_open() to read SAM/BAM/CRAM
	if ( (condition_streams[conditions] =
		xt_fopen(condition_files[conditions], "r")) == NULL )
	{
	    fprintf(stderr, ": Could not open %s for read: %s.\n",
		    condition_files[conditions], strerror(errno));
	    return EX_NOINPUT;
	}
    }
    
    // FIXME: discard unnecessary fields to improve performance
    mask = BL_GFF_FIELD_ALL;
    bl_gff_init(&feature);
    while ( bl_gff_read(&feature, features_stream, mask) == BL_READ_OK )
    {
	if ( strcmp(BL_GFF_TYPE(&feature), "gene") == 0 )
	{
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
	    }
	}
    }
    
    for (c = 2; c < conditions; ++c)
	fclose(condition_streams[c]);
    fclose(features_stream);
    
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s features.gff3 \\\n", argv[0]);
    fprintf(stderr, "\tcondition1.[sam|bam|cram][.gz|.bz2|.xz] \\\n");
    fprintf(stderr, "\tcondition2.[sam|bam|cram][.gz|.bz2|.xz]\n");
    exit(EX_USAGE);
}
