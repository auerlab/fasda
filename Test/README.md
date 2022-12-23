# Running the tests

## Required software

All tools required for this test can be installed via FreeBSD ports
or via pkgsrc on any other POSIX platform.

FreeBSD:

    pkg install rna-seq

pkgsrc:

    cd (pkgsrc-root)/biology/rna-seq
    bmake install

## Test procedure

Run the yeast-test.sh script with the number of replicates as an
argument between 3 and 48 (48 is the number of replicates in the test data):

    ./yeast-test.sh 3
    ./yeast-test.sh 48

See the main FASDA README.md for more information.
