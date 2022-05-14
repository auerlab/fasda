/***************************************************************************
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-05-14  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX 10

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    long    v[MAX], c, sum_v;
    double  lv[MAX], sum_lv, avg_v, avg_lv;
    
    srandom(time(NULL));
    for (c = 0, sum_v = 0, sum_lv = 0; c < MAX; ++c)
    {
	v[c] = random();
	lv[c] = log(v[c]);
	printf("%ld %f\n", v[c], lv[c]);
	sum_v += v[c];
	sum_lv += lv[c];
    }
    printf("avg v = %f  avg lv = %f\n",
	    avg_v = (double)sum_v / MAX, avg_lv = sum_lv / MAX);
    
    for (c = 0; c < MAX; ++c)
	printf("%f %f %f\n", v[c] / avg_v, exp(lv[c] - avg_lv),
		(exp(lv[c] - avg_lv)) / (v[c] / avg_v));
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s\n", argv[0]);
    exit(EX_USAGE);
}
