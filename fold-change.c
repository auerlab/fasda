/***************************************************************************
 *  Description:
 *      Calculate fold-change from normalized abundances.
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
