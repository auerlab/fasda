#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <xtend/math.h>

/***************************************************************************
 *  Description:
 *      Adapted from
 *      https://medium.com/enjoy-algorithm/
 *          find-all-possible-combinations-of-k-numbers-from-1-to-n-88f8e3fad33c
 *      Order of magnitude slower than k nested loops.
 *      This program takes 40 seconds on an i5 for 4000 choose 3
 *      A nested loop is instantaneous
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-09-30  Jason Bacon Begin
 ***************************************************************************/

void    combination_list(size_t list[],
			 size_t n, size_t k, size_t index, size_t start)

{
    int     i;
    
    //printf("index = %zu  start = %zu\n", index, start);
    if ( index == k )
    {
	/*
	for (i = 0; i < k; ++i)
	    printf("%zu ", list[i]);
	putchar('\n');
	*/
    }
    else
    {
	for (i = start; (i < n) && (n - i + 1 >= k - index); ++i)
	{
	    list[index] = i + 1;
	    combination_list(list, n, k, index + 1, i + 1);
	}
    }
}


int     main(int argc,char *argv[])

{
    size_t  list[100];
    size_t  n = 4060, k = 3, index = 0, start = 0, i;
    
    combination_list(list, n, k, index, start);
    for (i = 0; i < k; ++i)
	printf("%zu ", list[i]);
    putchar('\n');
    return EX_OK;
}
