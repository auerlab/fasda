#!/bin/sh -e

for abundance in Results/10-fasda-hisat2/*-abundance.tsv; do
    echo $abundance
    # Compute gene-level abundances
    transcript_list=${abundance%-abundance.tsv}-transcripts.txt
    cut -f 1 $abundance > $transcript_list

    transcript2gene_list=${abundance%-abundance.tsv}-transcript2gene.txt
    echo $transcript2gene_list

    blt ensemblid2gene \
	Results/04-reference/$(Reference/gff-filename.sh) \
	$transcript_list > $transcript2gene_list
    more $transcript2gene_list
    
    # Replace transcripts with genes and combine counts from the
    # same gene in the abundance file
done
