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

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    char        *features_file; //, *condition1_file, *condition2_file;
    FILE        *features_stream; //, condition1_stream, condition2_stream;
    bl_gff_t    feature;
    int         mask;
    
    if ( argc < 4 )
	usage(argv);
    features_file = argv[1];
    //condition1_file = argv[2];
    //condition2_file = argv[3];
    
    if ( (features_stream = fopen(features_file, "r")) == NULL )
    {
	fprintf(stderr, "diffanal: Could not open %s for read: %s.\n",
		features_file, strerror(errno));
	return EX_NOINPUT;
    }
    bl_gff_skip_header(features_stream);
    
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
	}
    }
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
