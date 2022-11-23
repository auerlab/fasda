#!/bin/sh -e


##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}


while true; do
    clear
    # https://www.ascii-art-generator.org/
    cat << EOM

			_________   _____ ____  ___ 
		       / ____/   | / ___// __ \/   |
		      / /_  / /| | \__ \/ / / / /| |
		     / __/ / ___ |___/ / /_/ / ___ |
		    /_/   /_/  |_/____/_____/_/  |_|

1.. View alignment counts for one feature
2.. View alignment counts based on DA results
Q.. Quit

EOM
    printf "Selection? "
    read selection
    case $selection in
    1)
	printf "Feature name? "
	read feature
	for replicates in $(seq 3 48); do
	    pr=$(printf "%02s" $replicates)
	    if [ -e Data/09-fasda/WT-all-norm-$pr.tsv ]; then
		./show-counts.sh $feature $replicates | more
	    fi
	done
	;;
    
    2)
	printf "Replicates? "
	read replicates
	printf "Maximum P-value? "
	read pval
	printf "Minimum fold-change? "
	read fc
	
	pr=$(printf "%02s" $replicates)
	da_file=Data/09-fasda/WT-SNF2-FC-NE-$pr.txt
	if [ -e $da_file ]; then
	    awk -v pval=$pval -v fc=$fc \
		'$8 <= pval && $7 >= fc { print $0 }' \
		$da_file | more
	else
	    printf "$da_file: No such file.\n"
	    pause
	fi
	;;
    
    Q|q)
	exit 0
	;;
    
    *)
	printf "Invalid selection.\n" >> /dev/stderr
	;;
    esac
done
