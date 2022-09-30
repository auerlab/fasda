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
 *  Generate all combinations n choose $k.  This is much faster than
 *  generic algorithms for generating n choose k lists for any n and k.
 */

void    combinations$k(size_t n)

{
EOM
    
    # Variable defs
    printf "    size_t  "
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

for c in $(seq 2 10); do
    gen_loop $c
done
