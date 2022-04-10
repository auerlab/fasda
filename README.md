# Diffanal

## Description

There are mature, easy-to-use, and efficient tools for all steps in a typical
differential analysis pipeline up to and including alignment (read mapping)
and peak calling.  Well-maintained tools such as FastQC, cutadapt, BWA,
Bowtie2, HISAT2, samtools, and MACS2 make the early stages of an RNA-Seq,
ATAC-Seq, or ChIP-Seq analysis fairly straightforward.

The differential analysis step itself has until now has been problematic,
with few well-developed tools available, and many of them
requiring fairly sophisticated R scripting for basic use.  Code maintenance
has also been an issue, with even the most popular tools falling into
disrepair at times and frequently presenting installation issues due to
incompatibility with the latest version of R or other dependencies.

Diffanal aims to provide a fast and simple differential analysis tool that
just works and does not require any knowledge beyond basic Unix command-line
skills.  The code is written entirely in C to maximize efficiency and
portability, and to provide a simple command-line user interface.

## Status

We're in the early stages of development.  So far we are able to efficiently
calculate fold-change for an arbitrary number of conditions using a
simplistic coverage measurement.  Once this code
is robust we will move onto computing P-values, exploring more sophisticated
coverage algorithms, and adding other features.

The sample output below is from Mus musculus annotations (GRCm39) and 3
hisat2 BAM files with 68, 89, and 76 million reads, representing three time
points during development.  This 3-condition differential analysis runs in
6 minutes on a Core i5 2.9GHz (ThinkCenter M92p-tiny) and uses about
2 megabytes (yes, megabytes - not gigabytes) of RAM. The analysis for just
conditions 1 and 2 runs in 4.5 minutes.

```
diffanal mouse-sorted.gff3 time1.bam time2.bam time3.bam

Ch Gene             Cond1  Cond2  Cond3  FC 1-2  FC 1-3  FC 2-3
 1 4933401J01Rik     0.00   0.00   0.00       *       *       *
 1 Xkr4              0.33   0.45   0.10    1.36    0.29    0.22
 1 Gm37180           0.00   0.00   0.00       *       *       *
 1 Gm37363           0.00   0.00   0.00       *       *       *
 1 Gm37686           0.00   0.00   0.00       *       *       *
 1 Gm37329           0.00   0.00   0.00       *       *       *
 1 Gm38148           0.00   0.00   0.00       *       *       *
 1 Gm10568           0.00   0.00   0.00       *       *       *
 1 Gm38385           0.27   0.10   0.13    0.39    0.47    1.21
 1 Rp1               0.03   0.04   0.02    1.29    0.68    0.52
 1 Gm37483           0.00   0.00   0.00       *       *       *
 1 Sox17             0.36   0.19   0.14    0.52    0.39    0.75
 1 Mrpl15           27.97  36.81  22.93    1.32    0.82    0.62
 1 Gm37144           0.00   0.00   0.00       *       *       *
 1 Lypla1           42.11  34.87  29.35    0.83    0.70    0.84
 1 Gm37988           0.15   0.06   0.09    0.40    0.60    1.51
 1 Tcea1             0.20   0.05   0.11    0.25    0.54    2.17
 1 Gm37277           0.62   0.00   1.56    0.00    2.53     inf
 1 Rgs20             3.02   3.30   3.13    1.09    1.04    0.95
 1 Gm37079           0.00   0.00   0.00       *       *       *
 1 Atp6v1h           1.19   1.47   1.35    1.24    1.13    0.92
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
of the nearly 20,000 packages in the collection.  The
[auto-pkgsrc-setup](http://netbsd.org/~bacon/) script can assist you with
basic setup.

First bootstrap pkgsrc using auto-pkgsrc-setup or any
other method.  Then run the following commands:

```
cd pkgsrc-dir/biology/diffanal
bmake install clean
```

There may also be binary packages available for your platform.  If this is
the case, you can install by running:

```
pkgin install diffanal
```

See the [Joyent Cloud Services Site](https://pkgsrc.joyent.com/) for
available package sets.

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
