# FASDA - Fast And Simple Differential Analysis

## Description

There are mature, easy-to-use, and efficient tools for all steps in a typical
differential analysis pipeline up to and including alignment (read mapping)
and peak calling.  Well-maintained tools such as FastQC, cutadapt, fastq-trim,
BWA, Bowtie2, HISAT2, samtools, bedtools, and MACS2 make the early stages of
an RNA-Seq, ATAC-Seq, or ChIP-Seq analysis fairly straightforward.

Data-massaging before and after differential analysis is also relatively
simple using standard Unix tools such as awk, cut, grep, sed, and tr, or
simple perl or python scripts.

The differential analysis step itself has been problematic,
however, with few well-developed tools available and many of them requiring
fairly sophisticated R scripting for basic use.  Code maintenance
has also been an issue, with even the most popular tools falling into
disrepair at times and frequently presenting installation issues due to
incompatibility with the latest version of R or other dependencies.

FASDA is a fast and simple differential analysis tool that
just works and does not require any knowledge beyond basic Unix command-line
skills.  It uses a simple command-line interface (CLI) analogous
to popular tools such as bedtools, BWA, and samtools.
The code is written entirely in C and built on
[biolibc](https://github.com/auerlab/biolibc) to maximize efficiency,
portability, and interoperability with other tools.

Starting with kallisto output, computing fold-change and associated P-values
for two conditions can typically be completed in a few seconds with three
simple commands:

```
fasda normalize --output c1-all-norm.tsv c1-*/abundance.tsv
fasda normalize --output c2-all-norm.tsv c2-*/abundance.tsv
fasda fold-change --output FC.txt c1-all-norm.tsv c2-all-norm.tsv
```

Another issue addressed by FASDA is the fact that most popular
differential analysis tools suffer from a high
false discovery rate (FDR) (Li, et al:
[https://doi.org/10.1186/s13059-022-02648-4](https://doi.org/10.1186/s13059-022-02648-4))

FASDA computes exact P-values for experiments with
fewer than 5 replicates and near-exact P-values for 5 to 7 replicates
(possible count pairs are down-sampled to control run time, but resulting
P-values are generally stable to 2 decimal places).  For 8 or more replicates,
we use the Mann-Whitney U-test (A.K.A. Wilcoxon rank-sum test),
a non-parametric test that provides high stability and low FDR.  The main
limitation of Mann-Whitney is that it requires a minimum of 8 samples
to achieve reasonable statistical power.  As a result, it is not
useful for typical RNA-Seq or ATAC-Seq experiments with only 3 replicates.

Parametric tests used by popular tools can provide reasonable power at very
low sample sizes, in exchange for high FDR when the data do not fit the
parametric assumptions well.  However, their high FDR and
poor performance for large sample sizes illustrates a need for
a new approach.  Hence, a major goal with FASDA is to fill an
under-served niche of high-sample studies with a tool that is fast and
produces stable results.

FASDA currently uses Mean Ratios Normalization (MRN) to normalize counts
prior to computing fold-changes and P-values.

## Status

FASDA is currently undergoing alpha-testing.

The kallisto abundance.tsv file format is used as input for
normalization and computing fold-change and P-values.  The "fasda abundance"
command computes abundances from SAM/BAM/CRAM files and
produces a kallisto-style abundance TSV file, so that data from other
aligners can be used.  This is currently used to support output from hisat2,
and will support ChIP-Seq and ATAC-Seq peak data in the future.

FASDA's abundance estimates are currently about 40% higher than kallisto's,
but they are proportional.  I would assume that kallisto's estimates are
more accurate, but the resulting fold-changes are the same.  This
will be explored further when time permits.

The sample output below is from 14 biological replicates of yeast RNA-Seq
data with wild-type and SNF2 mutant conditions.  The run times included
below show that FASDA is extremely fast, normalizing over 92,000 estimated
counts directly from kallisto abundance.tsv files in about 1/4 second
(on a 2.6 GHz i5 laptop) and computing fold-change and P-values in 1/20
second.
Memory use is also very low, with "fasda normalize" peaking around
66 megabytes and "fasda fold-change" around 9 megabytes for the yeast
example.

An automated pipeline script for reproducing these results is provided in
Test/yeast-test.sh.

```
Raw data:

File                    Reads
SNF2-1.fastq.gz       7541320
SNF2-10.fastq.gz      6189284
SNF2-11.fastq.gz      7442520
SNF2-12.fastq.gz      5549584
SNF2-13.fastq.gz      5294592
SNF2-14.fastq.gz      8147372
SNF2-2.fastq.gz       6376640
SNF2-3.fastq.gz       6378780
SNF2-4.fastq.gz       7909012
SNF2-5.fastq.gz       5487300
SNF2-6.fastq.gz       6351828
SNF2-7.fastq.gz       9388052
SNF2-8.fastq.gz       6077820
SNF2-9.fastq.gz       7298096
WT-1.fastq.gz         4375828
WT-10.fastq.gz        5857000
WT-11.fastq.gz        5069724
WT-12.fastq.gz        5895732
WT-13.fastq.gz        6717716
WT-14.fastq.gz        7145472
WT-2.fastq.gz         5870276
WT-3.fastq.gz         4912828
WT-4.fastq.gz         7707420
WT-5.fastq.gz         6053972
WT-6.fastq.gz         10851336
WT-7.fastq.gz         6421324
WT-8.fastq.gz         7078852
WT-9.fastq.gz         6623692

Normalizing WT...
	0.27 real         0.21 user         0.06 sys
Normalizing SNF2...
	0.27 real         0.25 user         0.02 sys
Computing fold-change...
	0.05 real         0.05 user         0.00 sys

Feature                 MNC1    MNC2  SD/C1  SD/C2  %Agr  FC 1-2  P-val
YPL071C_mRNA            36.1    65.2    0.6    0.7    85    1.80  0.035
YLL050C_mRNA           486.9   705.4    0.5    0.4    78    1.45  0.048
YMR172W_mRNA            80.3   157.4    0.6    0.8    78    1.96  0.008
YOR185C_mRNA            61.9    70.2    0.6    0.5    71    1.13  0.198
YLL032C_mRNA            42.4    32.8    0.7    0.3    64    0.77  0.335
YBR225W_mRNA            67.0   122.1    0.6    0.6    85    1.82  0.022
YEL041W_mRNA            22.6    29.5    0.6    0.4    64    1.30  0.081
YOR237W_mRNA            17.1    46.4    0.6    0.8    92    2.72  0.000
YMR027W_mRNA           203.8   294.6    0.5    0.4    71    1.45  0.022
YBR182C-A_mRNA           0.1     0.2    3.5    2.4    85    3.04  0.713
YKL103C_mRNA           195.6   388.8    0.5    0.4    92    1.99  0.001
YOL048C_mRNA            64.0   131.8    0.7    0.4    92    2.06  0.001
YIR015W_mRNA            13.6    25.7    0.8    0.7    85    1.88  0.009
YNR017W_mRNA           258.0   257.3    0.5    0.4    57    1.00  0.854
YBL055C_mRNA           106.7    91.3    0.6    0.3    57    0.86  0.783
...
```

## Interpreting Results

The fasda fold-change output shows the mean normalized counts across
replicates, the standard deviation from that mean (divided by the mean
to show % deviation rather than absolute), the % agreement among replicates
as to the direction of the change (up-regulated vs down or no change),
the fold-change computed from total counts across replicates,
and the exact, near-exact, or Mann-Whitney P-value for this set of counts.

All of these values, along with any available knowledge of the biology,
should be considered when deciding whether a given change is significant.

P-values from any differential analysis tool should never be taken too
seriously. There are countless uncontrollable biological variables that
can affect the RNA abundance in a cell.  There are also numerous sources
of experimental error in sample prep and sequencing that can lead to large
variations in read counts.  Technical replicates (replicates from
the same biological sample) and spike-in controls can reveal some of these
technical issues, but do not address biological variations.

Another problem is that many biology experiments use only 3 replicates.
We simply cannot draw much confidence from any statistics based on 3
samples.  Schurch, et al recommend a minimum of 6 biological replicates,
12 in order to reliably identify DEGs with fold change less than 4
([https://doi.org/10.1261%2Frna.053959.115](https://doi.org/10.1261%2Frna.053959.115).

Note that higher fold-changes and higher counts
lead to higher statistical significance,
but this does not necessarily correlate to higher biological significance.
P-value calculations typically make the same assumptions about all genes.
A fold-change of 1.5 may be highly significant for one gene while meaningless
for another for entirely biochemical reasons.  E.g. increasing the abundance
of a scarce enzyme by 50% could have a profound effect on a cell, while a 50%
increase in many other proteins has little effect.
Statistical routines have no knowledge of the biology that determines this.

There is huge variability on the computational side as well.
Well-established differential analysis tools commonly report very different
sets of genes as differentially expressed.  Li, et al
([https://doi.org/10.1186/s13059-022-02648-4](https://doi.org/10.1186/s13059-022-02648-4))
reported that 23.71% to 75% of
the DEGs identified by DESeq2 were missed by edgeR.  In one data set tested,
DESeq2 and edgeR had only an 8% overlap in the DEGs they identified.

Hence, simply assuming that P-values < 0.05 represent significant
changes while others do not would be foolish.  Rather than try too hard
to produce adjusted P-values that you can take on faith, we provide simple,
honest statistics and leave it to you to consider them carefully.

From our preliminary experiments, we have found that P-values
between 0.05 and 0.20 are often questionable and may warrant a closer look
at the raw data.  A quick look at the individual read counts and
variance for each replicate in these cases
can be enlightening.  You will usually see a high variance in counts across
individuals, sometimes with up-regulation in some individuals and
down-regulation in others.

For example, consider the following gene, for which Sleuth reported a
P-value of 0.0132, while the exact P-value computed by FASDA is 0.116.
The kallisto estimated counts show that
one replicate was up-regulated almost 4-fold, another almost 2-fold, and
the third was slightly down-regulated.  A P-value of 0.01 would not
likely make anyone suspect this situation.  0.116, on the other hand,
tells us that there is a good chance this is significant, but maybe we
should take a minute or two to look at the read counts and consider the
biology behind them.  This is a small investment that will help us better
decide whether a costly experimental verification is warranted.

```
		    FASDA                     Sleuth
	   Feature   Cond1   Cond2  FC   1-2    SC1    SC2  SFC    SPV
ENSMUST00000017610  7620.5 16006.7 2.1 0.116  191.4  433.5  2.3 0.0132

Kallisto estimated counts:

	     R1      R2      R3
Cond1   5382.16  8567.4 6986.43
Cond2   21519.9 16196.7 6307.52
```
Conversely, Sleuth produced a P-value of 0.12 for the following, which
looks like a slam-dunk given the high counts and strong agreement
among all replicates.

```
		    FASDA                     Sleuth
	   Feature   Cond1   Cond2  FC   1-2    SC1    SC2  SFC    SPV
ENSMUST00000036928  1165.0  4064.3 3.5 0.024   70.5  300.4  4.3 0.1203

	     R1      R2      R3
Cond1   1045.41 942.707 1220.02
Cond2   2835.7  4167.81 3718.39
```

The bottom line is, while the 0.05 rule is a good one mathematically,
it is only as reliable as the input data and the statistical methods
used to analyze it..
We cannot count on experimental and computational results reflecting the
biology with that much accuracy.  Give the results of any differential
analysis a generous margin of error, and examine the data more closely for
anything near that margin.

The question then becomes how to narrow down the list of differential
features and identify those of interest.  Our intent is to provide
additional tools that facilitate examination of the results using
multiple criteria rather than filtering first by P-value alone.  For
instance, we might select for genes with a P-value < 0.2, low variance
in the alignment counts, and known relation to the phenotype being
studied.  We might also tolerate lower fold-changes for genes whose
effect is known to be biologically amplified.

For ChIP-Seq or ATAC-Seq , we might select for peaks with a
P-value < 0.2, low variance in the counts, and proximity to a gene of
interest.  Multiple facets must be considered in order to avoid missing
important features that have a P-value above 0.05 due to experimental
imprecision rather than a true lack of biological significance.  Likewise,
we need a multifaceted approach to eliminate false positives.

## Design and Implementation

The code is organized following basic object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with accessor macros and mutator functions
provided, so dependent applications and libraries need not access
structure members directly.  Since the C language cannot enforce this, it's
up to application programmers to exercise self-discipline.

For detailed coding standards, see
https://github.com/outpaddling/Coding-Standards/.

## Building and installing

FASDA is intended to build cleanly in any POSIX environment on any CPU
architecture.  Please don't hesitate to open an issue if you encounter
problems on any Unix-like system.

Primary development is done on FreeBSD with clang, but the code is frequently
tested on Linux, MacOS, NetBSD, and OpenIndiana as well.  MS Windows is not supported,
unless using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

The Makefile is designed to be friendly to package managers, such as
[Debian packages](https://www.debian.org/distrib/packages),
[FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [pkgsrc](https://pkgsrc.org/), etc.
End users should install using a package manager.  Note that pkgsrc can be used by anyone, on virtually any POSIX operating system, with or without administrator privileges..

I maintain a FreeBSD port and a pkgsrc package, which is sufficient to install
cleanly on virtually any POSIX platform.  If you would like to see a
package in another package manager, please consider creating a package
yourself.  This will be one of the easiest packages in the collection and
hence a good vehicle to learn how to create packages.

For an overview of available package managers, see the
[Repology website](https://repology.org/).

### Installing FASDA on FreeBSD:

FreeBSD is a highly underrated platform for scientific computing, with over
2,000 scientific libraries and applications in the FreeBSD ports collection
(of more than 30,000 total), modern clang compiler, fully-integrated ZFS
file system, and renowned security, performance, and reliability.
FreeBSD has a somewhat well-earned reputation for being difficult to set up
and manage compared to user-friendly systems like [Ubuntu](https://ubuntu.com/).
However, if you're a little bit Unix-savvy, you can very quickly set up a
workstation, laptop, or VM using
[desktop-installer](http://www.acadix.biz/desktop-installer.php).
[GhostBSD](https://ghostbsd.org/) offers an experience very similar
to Ubuntu, but is built on FreeBSD rather than Debian Linux.  GhostBSD
packages lag behind FreeBSD ports slightly, but this is not generally
an issue and there are workarounds.

To install the binary package on FreeBSD:

```
pkg install fasda
```

You can just as easily build and install from source.  This is useful for
FreeBSD ports with special build options, for building with non-portable
optimizations such as -march=native, building with debugging info, and for 
[work-in-progress ports](https://github.com/outpaddling/freebsd-ports-wip),
for which binary packages are not yet maintained.

```
cd /usr/ports/biology/fasda && env CFLAGS='-march=native -O2' make install
``` 

### Installing via pkgsrc

pkgsrc is a cross-platform package manager that works on any Unix-like
platform. It is native to [NetBSD](https://www.netbsd.org/) and well-supported
on [Illumos](https://illumos.org/), [MacOS](https://www.apple.com/macos/),
[RHEL](https://www.redhat.com)/[CentOS](https://www.centos.org/), and
many other Linux distributions.
Using pkgsrc does not require admin privileges.  You can install a pkgsrc
tree in any directory to which you have write access and easily install any
of the nearly 20,000 packages in the collection.

The
[auto-pkgsrc-setup](https://github.com/outpaddling/auto-admin/blob/master/User-scripts/auto-pkgsrc-setup)
script will help you install pkgsrc in about 10 minutes.  Just download it
and run

```
sh auto-pkgsrc-setup
```

Then, assuming you selected current packages and the default prefix

```
source ~/Pkgsrc/pkg/etc/pkgsrc.sh   # Or pkgsrc.csh for csh or tcsh
cd ~/Pkgsrc/pkgsrc/biology/fasda
sbmake install clean clean-depends
```

See the pkgsrc documentation for more information.

Community support for pkgsrc is available through the
[pkgsrc-users](http://netbsd.org/mailinglists) mailing list.
