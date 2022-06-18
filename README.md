# Diffanal

## Description

There are mature, easy-to-use, and efficient tools for all steps in a typical
differential analysis pipeline up to and including alignment (read mapping)
and peak calling.  Well-maintained tools such as FastQC, cutadapt, fastq-trim,
BWA, Bowtie2, HISAT2, samtools, bedtools, and MACS2 make the early stages of
an RNA-Seq, ATAC-Seq, or ChIP-Seq analysis fairly straightforward.

Data massaging before and after differential analysis is also relatively
easy using standard Unix tools such as awk, cut, grep, sed, and tr, or
simple perl or python scripts.

The differential analysis step itself has until now has been problematic,
however, with few well-developed tools available, and many of them
requiring fairly sophisticated R scripting for basic use.  Code maintenance
has also been an issue, with even the most popular tools falling into
disrepair at times and frequently presenting installation issues due to
incompatibility with the latest version of R or other dependencies.

Diffanal aims to provide a fast and simple differential analysis tool that
just works and does not require any knowledge beyond basic Unix command-line
skills.  The code is written entirely in C and built on
[biolibc](https://github.com/auerlab/biolibc) to maximize efficiency and
portability, and to provide a simple command-line user interface.

Starting with kallisto output, computing fold-change and associated P-values
for two conditions can be done in seconds with three simple commands:

```
diffanal normalize --output c1-all-norm.tsv c1-*/abundance.tsv
diffanal normalize --output c2-all-norm.tsv c2-*/abundance.tsv
diffanal fold-change --output FC.txt c1-all-norm.tsv c2-all-norm.tsv
```

Another issue is that most popular differential analysis tools show a high
false discovery rate (FDR) regardless of sample size (biological replicates):
[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4)

Diffanal utilizes the Mann-Whitney U-test (A.K.A. Wilcoxon rank-sum test), a
non-parametric test that provides high stability and low FDR.  The main
limitation is that Mann-Whitney
requires a minimum sample size of about 8 to achieve reasonable statistical
power and hence is not useful for typical RNA-Seq or ATAC-Seq experiments
with only 3 replicates.

Parametric tests used by popular tools can provide reasonable power at very
low sample sizes in exchange for high FDR. However, their high FDR and
otherwise poor performance for large sample sizes illustrates a need for
a new approach.  Our first goal with diffanal is to address this need
order to fill an under-served niche.  Additional use cases including
data with fewer replicates may be addressed at a later date.

## Status

We're still in the fairly early stages of development.  We are able
to normalize counts using Mean Ratios Normalization (MRN) and compute
fold-change and Mann-Whitney P-values for an arbitrary number of conditions.

Currently only kallisto abundance.tsv files can be used as input.  A
tool to compute abundances from SAM/BAM/CRAM files and produce a kallisto
style abundance file is in the works so that RNA-Seq data from other aligners
can eventually be used.  We also plan to eventually support ChIP and ATAC
peak data.

The sample output below is from 14 biological replicates of yeast RNA-Seq
data with wild-type and SNF2 mutant conditions.  The run times included
below show that diffanal is extremely fast, normalizing over 92,000 estimated
counts directly from kallisto abundance.tsv files in about 1/4 second
(on a 2.6 GHz i5 laptop) and computing fold-change and P-values in 1/20
second.
Memory use is also very low, with "diffanal normalize" peaking around
66 megabytes and "diffanal fold-change" around 9 megabytes for the yeast
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

Feature                           Cond1    Cond2  FC 1-2  P-val
YPL071C_mRNA                      43.60    81.96    1.88 0.0169
YLL050C_mRNA                     603.23   865.14    1.43 0.0432
YMR172W_mRNA                      81.78   157.15    1.92 0.0131
YOR185C_mRNA                      76.50    88.00    1.15 0.3346
YLL032C_mRNA                      45.14    35.93    0.80 0.5814
YBR225W_mRNA                      67.91   121.57    1.79 0.0274
YEL041W_mRNA                      24.57    33.46    1.36 0.0731
YOR237W_mRNA                      18.87    48.77    2.58 0.0003
YMR027W_mRNA                     217.29   312.18    1.44 0.0244
YBR182C-A_mRNA                     0.22     0.70    3.19 0.1753
YKL103C_mRNA                     207.99   414.52    1.99 0.0006
YOL048C_mRNA                      72.43   147.76    2.04 0.0007
...
```

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

Diffanal is intended to build cleanly in any POSIX environment on any CPU
architecture.  Please don't hesitate to open an issue if you encounter
problems on any Unix-like system.

Primary development is done on FreeBSD with clang, but the code is frequently
tested on Linux, MacOS, and NetBSD as well.  MS Windows is not supported,
unless using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

The Makefile is designed to be friendly to package managers, such as
[Debian packages](https://www.debian.org/distrib/packages),
[FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [pkgsrc](https://pkgsrc.org/), etc.
End users should install via one of these if at all possible.

I maintain a FreeBSD port and a pkgsrc package, which is sufficient to install
cleanly on virtually any POSIX platform.  If you would like to see a
package in another package manager, please consider creating a package
yourself.  This will be one of the easiest packages in the collection and
hence a good vehicle to learn how to create packages.

For an overview of available package managers, see the
[Repology website](https://repology.org/).

### Installing diffanal on FreeBSD:

FreeBSD is a highly underrated platform for scientific computing, with over
2,000 scientific libraries and applications in the FreeBSD ports collection
(of more than 30,000 total), modern clang compiler, fully-integrated ZFS
file system, and renowned security, performance, and reliability.
FreeBSD has a somewhat well-earned reputation for being difficult to set up
and manage compared to user-friendly systems like [Ubuntu](https://ubuntu.com/).
However, if you're a little bit Unix-savvy, you can very quickly set up a
workstation, laptop, or VM using
[desktop-installer](http://www.acadix.biz/desktop-installer.php).

To install the binary package on FreeBSD:

```
pkg install diffanal
```

You can just as easily build and install from source.  This is useful for
FreeBSD ports with special build options, for building with non-portable
optimizations such as -march=native, building with debugging info, and for 
[work-in-progress ports](https://github.com/outpaddling/freebsd-ports-wip),
for which binary packages are not yet maintained.

```
cd /usr/ports/biology/diffanal && env CFLAGS='-march=native -O2' make install
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
[auto-pkgsrc-setup](https://github.com/outpaddling/auto-admin/blob/master/Scripts/auto-pkgsrc-setup)
script will help you install pkgsrc in about 10 minutes.  Just download it
and run

```
sh auto-pkgsrc-setup
```

Then, assuming you selected current packages and the default prefix

```
source ~/Pkgsrc/pkg/etc/pkgsrc.sh   # Or pkgsrc.csh for csh or tcsh
cd ~/Pkgsrc/biology/diffanal
sbmake install clean clean-depends
```

See the pkgsrc documentation for more information.

Community support for pkgsrc is available through the
[pkgsrc-users](http://netbsd.org/mailinglists) mailing list.

### Building diffanal locally

Below are cave man install instructions for development purposes, not
recommended for regular use.
Diffanal depends on [biolibc](https://github.com/auerlab/biolibc).
Install biolibc before attempting to build diffanal.

1. Clone the repository
2. Run "make depend" to update Makefile.depend
3. Run "make install"

The default install prefix is ../local.  Clone diffanal, biolibc and dependent
apps into sibling directories so that ../local represents a common path to all
of them.

To facilitate incorporation into package managers, the Makefile respects
standard make/environment variables such as CC, CFLAGS, PREFIX, etc.  

Add-on libraries required for the build, such as biolibc, should be found
under ${LOCALBASE}, which defaults to ../local.
The library, headers, and man pages are installed under
${DESTDIR}${PREFIX}.  DESTDIR is empty by default and is primarily used by
package managers to stage installations.  PREFIX defaults to ${LOCALBASE}.

To install directly to /myprefix, assuming biolibc is installed there as well,
using a make variable:

```
make LOCALBASE=/myprefix clean depend install
```

Using an environment variable:

```
# C-shell and derivatives
setenv LOCALBASE /myprefix
make clean depend install

# Bourne shell and derivatives
LOCALBASE=/myprefix
export LOCALBASE
make clean depend install
```

View the Makefile for full details.
