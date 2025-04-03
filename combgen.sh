#!/bin/sh -e

##########################################################################
#   Synopsis:
#       
#   Description:
#       
#   Arguments:
#       
#   Returns:
#
#   Examples:
#
#   Files:
#
#   Environment:
#
#   See also:
#       
#   History:
#   Date        Name        Modification
#   2022-09-30  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 k\n"
    exit 1
}


##########################################################################
#   Function description:
#       
#   Arguments:
#       
#   Returns:
#       
#   History:
#   Date        Name        Modification
#   2022-09-30  Jason Bacon Begin
##########################################################################

print_indent()
{
    if [ $# != 1 ]; then
	printf "Usage: indent n\n"
	exit 1
    fi
    c=$1
    
    for space in $(seq 1 $((c - 1))); do
	printf " "
    done
}


##########################################################################
#   Function description:
#       
#   Arguments:
#       
#   Returns:
#       
#   History:
#   Date        Name        Modification
#   2022-09-30  Jason Bacon Begin
##########################################################################

gen_loop()
{
    local k=$1
    
    # Downsample set of all possible means for replicates > 5 to minimize
    # run time.  Use the highest increment possible to minimize run time
    # while getting close to stability to 2 decimal places.
    # Gradually increase passes until results are stable to 2 decimal places.
    # Inputs chosen to generate P-values that are often > 0.05.  E.g.
    # with means counts 100, 200, we need a max deviation of about .8.
    # FIXME: What is a good max deviation representative of real data?
    # Modify fold-change to report stats on counts
    # Span correlates with P-value.  Span should be < 0.01 when P-values
    # are near or below 0.05.  (Ballpark guestimate for stable averages
    # of 5 P-values)
    case $k in
    5)
	increment=3 # 3, 8: p-values mostly stable to 2 decimal places
	passes=8    # Exact P-value (inc=1, no srandom) = 0.393 ~1 sec
		    # 100 200 .8 5 1
		    # Big difference in P-values between increment 2 and 3
	;;
    6)
	increment=5 # 4, 3: p-values mostly stable to 2 decimal places
	passes=8    # Exact P-value (inc=1, no srandom) = 0.117 ~2 min
		    # 100 200 .8 6 1
	;;
    7)
	increment=7 # 6, 1: p-values mostly stable to 2 decimal places
	passes=4    # 
		    # 100 200 .8 7 1
	;;
    
    8)
	increment=12    # 12, 6: p-values mostly stable to 2 decimal places
	passes=6        # 
			# 100 200 .8 8 1
	;;
    
    9)
	increment=17    # 17, 12: p-values mostly stable to 2 decimal places
	passes=12
			# 100 200 .8 9 1
	;;
    
    10)
	increment=22    # 22, 12: p-values mostly stable to 2 decimal places
	passes=12
			# 100 200 .8 10 1
	;;
    
    11)
	increment=26    # 22, 12: p-values mostly stable to 2 decimal places
	passes=8
			# 100 200 1 11 1
	;;
    
    12)
	increment=30    # 30, 4: p-values mostly stable to 2 decimal places
	passes=4
			# 100 200 1 12 1
	;;

    *)
	increment=1
	passes=1
	;;
    esac
    
    cat << EOM

extern bool Debug;

/*
 *  Generate all combinations n choose $k and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   extreme_fcs$k(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0,
		    increment = $increment, pass, count = 0;
    double          c2_sum, c1_sum, fc;
    
EOM
    # Variable defs
    printf "    unsigned long  "
    for c in $(seq 1 $((k - 1))); do
	printf "c%d, " $c
    done
    printf "c%d;\n\n" $k
    
    # Nested loop
    printf "    for (pass = 0; pass < $passes; ++pass)\n"
    # printf "     #pragma omp parallel for shared(fc_ge, fc_le, count) private(c1_sum, c2_sum, fc)\n"
    printf "     for (c1 = 0; c1 < pair_count; c1 += increment)\n"
    for c in $(seq 2 $k); do
	print_indent $c
	if [ $increment -gt 1 ]; then
	    printf "     for (c$c = c$((c - 1)) + 1 + random() %% increment; c$c < pair_count; c$c += increment)\n"
	else
	    printf "     for (c$c = c$((c - 1)) + 1; c$c < pair_count; c$c += increment)\n"
	fi
    done
    
    # Body
    print_indent $c
    printf "     {\n"
    print_indent $c
    printf "         c2_sum =\n"
    for c2 in $(seq 1 $((k - 1))); do
	print_indent $c
	printf "             count_pairs[c$c2].c2_count +\n"
    done
    print_indent $c
    printf "             count_pairs[c$k].c2_count;\n"
    print_indent $c
    printf "         c1_sum =\n"
    for c2 in $(seq 1 $((k - 1))); do
	print_indent $c
	printf "             count_pairs[c$c2].c1_count +\n"
    done
    print_indent $c
    printf "             count_pairs[c$k].c1_count;\n"
    print_indent $c
    printf "         fc = c2_sum / c1_sum;\n"
    # print_indent $c
    # printf '             if ( count %% 100000000 == 0 ) fprintf(stderr, "%%lu\\r", count);\n'

    print_indent $c
    printf "         if ( fc >= observed_fc ) ++fc_ge;\n"
    print_indent $c
    printf "         else if ( fc <= 1.0 / observed_fc ) ++fc_le;\n"

    print_indent $c
    printf "         ++count;\n"
    print_indent $c
    printf "     }\n"
    
    # Closing braces
    cat << EOM
    if ( Debug )
	printf("FCs > %0.3f = %-5lu     FCs < %0.3f = %-5lu\n",
		observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    return fc_ge + fc_le;
}
EOM
}


##########################################################################
#   Main
##########################################################################

cat << EOM
/*
 *  Generated by combgen.sh.  DO NOT EDIT.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "exact-p-val.h"

EOM
# FC-count choose 11 overflows a 64-bit integer
max_reps=12
for c in $(seq 2 $max_reps); do
    gen_loop $c
done

cat << EOM


unsigned long   extreme_fcs_count(count_pair_t count_pairs[], unsigned long pair_count,
		      unsigned long replicates, double observed_fc,
		      unsigned long *fc_count)

{
    static unsigned long (*extreme_fcs_funcs[])(count_pair_t count_pairs[],
				    unsigned long pair_count,
				    double observed_fc,
				    unsigned long *fc_count) =
    {
EOM

for c in $(seq 2 $((max_reps - 1))); do
    printf "        extreme_fcs$c,\n"
done
printf "        extreme_fcs$max_reps\n    };\n"

cat << EOM
    unsigned long  func_index = replicates - 2;
    
    return extreme_fcs_funcs[func_index](count_pairs, pair_count,
				   observed_fc, fc_count);
}
EOM
