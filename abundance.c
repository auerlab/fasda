/***************************************************************************
 *  Description:
 *      Calculate abundances from alignment data.
 *      Output to a TSV file similar to Kallisto abundances.tsv so we
 *      have a common starting point for normalization and differential
 *      analysis.
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-05-04  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    switch(argc)
    {
	case 1:
	    break;
	
	default:
	    usage(argv);
    }
    
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s\n", argv[0]);
    exit(EX_USAGE);
}
