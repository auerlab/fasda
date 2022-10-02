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
    
    cat << EOM

/*
 *  Generate all combinations n choose $k and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge$k(double fc_list[], unsigned long fc_count,
			double observed_fc_mean)

{
    unsigned long   fc_ge = 0;
    double          fc_mean;
    
EOM
    # Variable defs
    printf "    unsigned long  "
    for c in $(seq 1 $((k - 1))); do
	printf "c%d, " $c
    done
    printf "c%d;\n\n" $k
    
    # Nested loop
    printf "    for (c1 = 0; c1 < fc_count; ++c1)\n"
    for c in $(seq 2 $k); do
	print_indent $c
	printf "    for (c$c = c$((c - 1)) + 1; c$c < fc_count; ++c$c)\n"
    done
    
    # Body
    print_indent $c
    printf "    {\n"
    print_indent $c
    printf "        fc_mean = ("
    for c2 in $(seq 1 $((k - 1))); do
	printf "fc_list[c$c2] + "
    done
    printf "fc_list[c$k]) / $k;\n"
    print_indent $c
    printf "        if ( fc_mean >= observed_fc_mean ) ++fc_ge;\n"
    print_indent $c
    printf "    }\n"
    
    # Closing braces
    printf "    return fc_ge;\n}\n"
}


##########################################################################
#   Main
##########################################################################

cat << EOM
// FIXME: Just a skeleton
void    fc_count(void);

EOM
max_reps=8
for c in $(seq 2 $max_reps); do
    gen_loop $c
done

cat << EOM


unsigned long   fc_ge_count(double fc_list[], unsigned long fc_count,
		      unsigned long replicates, double observed_fc_mean)

{
    static unsigned long (*fc_ge_funcs[])(double fc_list[],
				    unsigned long fc_count,
				    double observed_fc_mean) =
    {
EOM

for c in $(seq 2 $((max_reps - 1))); do
    printf "        fc_ge$c,\n"
done
printf "        fc_ge$max_reps\n    };\n"

cat << EOM
    unsigned long  func_index = replicates - 2;
    
    return fc_ge_funcs[func_index](fc_list, fc_count, observed_fc_mean);
}
EOM
