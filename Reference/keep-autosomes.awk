#############################################################################
#   Description:
#       Remove non-autosomal info from mouse cDNA fasta file
#
#   History: 
#   Date        Name        Modification
#   2019-10-05  Jason Bacon Begin
#############################################################################

BEGIN {
    FS=":";
}
{
    # As long as the description line does not indicate a numeric chromosome,
    # discard it and everything to the next sequence line.
    while ( ($0 ~ "^>") && ($3 !~ /^[0-9]+$/) )
    {
	do
	{
	    status=getline
	}   while ( (status == 1) && ($0 !~ "^>") );
    }
    print $0
}
