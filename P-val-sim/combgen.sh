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

// FIXME: Just a skeleton
void    fc_count(void);

/*
 *  Generate all combinations n choose $k.  This is much faster than
 *  generic algorithms for generating n choose k lists for any n and k.
 */

void    combinations$k(unsigned long n)

{
EOM
    
    # Variable defs
    printf "    unsigned long  "
    for c in $(seq 1 $((k - 1))); do
	printf "c%d, " $c
    done
    printf "c%d;\n\n" $k
    
    # Nested loop
    printf "    for (c1 = 0; c1 < n; ++c1)\n"
    for c in $(seq 2 $k); do
	print_indent $c
	printf "    for (c$c = c$((c - 1)) + 1; c$c < n; ++c$c)\n"
    done
    
    # Body
    print_indent $c
    printf "    {\n"
    print_indent $c
    printf "        fc_count();\n"
    print_indent $c
    printf "    }\n"
    
    # Closing braces
    printf "}\n"
}


##########################################################################
#   Main
##########################################################################

max_reps=20
for c in $(seq 2 $max_reps); do
    gen_loop $c
done

cat << EOM


void    combinations(unsigned long n, unsigned long k)

{
    static void    (*comb_funcs[])(unsigned long n) =
    {
EOM

for c in $(seq 2 $((max_reps - 1))); do
    printf "        combinations$c,\n"
done
printf "        combinations$max_reps\n    };\n"

cat << EOM
    unsigned long  func_index = k - 2;
    
    comb_funcs[func_index](n);
}
EOM
