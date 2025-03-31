############################################################################
#
#              Another Programmer's Editor Makefile Template
#
# This is a template Makefile for a simple program.
# It is meant to serve as a starting point for creating a portable
# Makefile, suitable for use under ports systems like *BSD ports,
# MacPorts, Gentoo Portage, etc.
#
# The goal is a Makefile that can be used without modifications
# on any Unix-compatible system.
#
# Variables that are conditionally assigned (with ?=) can be overridden
# by the environment or via the command line as follows:
#
#       make VAR=value
#
# For example, MacPorts installs to /opt/local instead of the default
# ../local, and hence might use the following:
# 
#       make PREFIX=/opt/local
#
# Different systems may also use different compilers and keep libraries in
# different locations:
#
#       make CC=gcc CFLAGS=-O2 LDFLAGS="-L/usr/X11R6 -lX11"
#
# Variables can also inheret values from parent Makefiles (as in *BSD ports).
#
# Lastly, they can be overridden by the environment, e.g.
# 
#       setenv CFLAGS "-O -Wall -pipe"  # C-shell and derivatives
#       export CFLAGS="-O -Wall -pipe"  # Bourne-shell and derivatives
#       make
#
# All these override methods allow the Makefile to respect the environment
# in which it is used.
#
# You can append values to variables within this Makefile (with +=).
# However, this should not be used to add compiler-specific flags like
# -Wall, as this would disrespect the environment.
#
#   History: 
#   Date        Name        Modification
#   2022-04-04  Jason Bacon Begin
############################################################################

############################################################################
# Installed targets

BIN     = fasda
LIBEXEC = abundance normalize fold-change pval-sim
SCRIPTS = filter
BINS    = ${BIN} ${LIBEXEC}

############################################################################
# List object files that comprise BIN.

OBJS_FASDA          = fasda.o alignment-stats-mutators.o
OBJS_ABUNDANCE      = abundance.o
OBJS_NORMALIZE      = normalize.o
OBJS_FOLD_CHANGE    = fold-change.o mann-whitney.o exact-p-val.o fc-ge.o
OBJS_PVAL_SIM       = pval-sim.o fc-ge.o mann-whitney.o exact-p-val.o

############################################################################
# Compile, link, and install options

# Install in ../local, unless defined by the parent Makefile, the
# environment, or a command line option such as PREFIX=/opt/local.
# FreeBSD ports sets this to /usr/local, MacPorts to /opt/local, etc.
PREFIX      ?= ../local

# Where to find local libraries and headers.  If you want to use libraries
# from outside ${PREFIX} (not usually recommended), you can set this
# independently.
LOCALBASE   ?= ${PREFIX}

# Allow caller to override either MANPREFIX or MANDIR
MANPREFIX   ?= ${PREFIX}
MANDIR      ?= ${MANPREFIX}/man
# FIXME: Need to realpath this if relative (e.g. ../local) or fasda won't
# find subcommands from arbitrary CWD
# Currently must use cave-man-install.sh for this until a bmake/gmake
# portable method is found
LIBEXECDIR  ?= ${PREFIX}/libexec/fasda

############################################################################
# Build flags
# Override with "make CC=gcc", "make CC=icc", etc.
# Do not add non-portable options (such as -Wall) using +=
# Make sure all compilers are part of the same toolchain.  Do not mix
# compilers from different vendors or different compiler versions unless
# you know what you're doing.

# Defaults that should work with GCC and Clang.
CC          ?= cc
CFLAGS      ?= -Wall -g -O
CFLAGS      += -DLIBEXECDIR=\"${LIBEXECDIR}\"
CFLAGS      += -DVERSION=\"`./version.sh`\"

# Link command:
# Use ${FC} to link when mixing C and Fortran
# Use ${CXX} to link when mixing C and C++
# When mixing C++ and Fortran, use ${FC} and -lstdc++ or ${CXX} and -lgfortran
LD          = ${CC}

CPP         ?= cpp

AR          ?= ar
RANLIB      ?= ranlib

INCLUDES    += -isystem ${PREFIX}/include -isystem ${LOCALBASE}/include
CFLAGS      += ${INCLUDES}
LDFLAGS     += -L${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib \
	       -L${LOCALBASE}/lib -Wl,-rpath,${LOCALBASE}/lib \
	       -lbiolibc -lxtend -lm

############################################################################
# Assume first command in PATH.  Override with full pathnames if necessary.
# E.g. "make INSTALL=/usr/local/bin/ginstall"
# Do not place flags here (e.g. RM = rm -f).  Just provide the command
# and let flags be specified separately.

CP      ?= cp
MV      ?= mv
MKDIR   ?= mkdir
LN      ?= ln
RM      ?= rm

# No full pathnames for these.  Allow PATH to dtermine which one is used
# in case a locally installed version is preferred.
PRINTF  ?= printf
INSTALL ?= install
STRIP   ?= strip

############################################################################
# Standard targets required by package managers

.PHONY: all depend clean realclean install install-strip help

all:    ${BINS}

fasda: ${OBJS_FASDA}
	${LD} -o fasda ${OBJS_FASDA} ${LDFLAGS}

abundance: ${OBJS_ABUNDANCE}
	${LD} -o abundance ${OBJS_ABUNDANCE} ${LDFLAGS}
	
normalize: ${OBJS_NORMALIZE}
	${LD} -o normalize ${OBJS_NORMALIZE} ${LDFLAGS}
	
fold-change: ${OBJS_FOLD_CHANGE}
	${LD} -o fold-change ${OBJS_FOLD_CHANGE} ${LDFLAGS}

pval-sim:   ${OBJS_PVAL_SIM}
	${LD} -o pval-sim ${OBJS_PVAL_SIM} ${LDFLAGS}

fc-ge.c:    combgen.sh
	./combgen.sh > fc-ge.c

############################################################################
# Include dependencies generated by "make depend", if they exist.
# These rules explicitly list dependencies for each object file.
# See "depend" target below.  If Makefile.depend does not exist, use
# generic source compile rules.  These have some limitations, so you
# may prefer to create explicit rules for each target file.  This can
# be done automatically using "cpp -M" or "cpp -MM".  Run "man cpp"
# for more information, or see the "depend" target below.

# Rules generated by "make depend"
# If Makefile.depend does not exist, "touch" it before running "make depend"
include Makefile.depend

############################################################################
# Self-generate dependencies the old-fashioned way
# Edit filespec and compiler command if not using just C source files

depend:
	rm -f Makefile.depend
	for file in *.c; do \
	    ${CC} ${INCLUDES} -MM $${file} >> Makefile.depend; \
	    ${PRINTF} "\t\$${CC} -c \$${CFLAGS} $${file}\n\n" >> Makefile.depend; \
	done

############################################################################
# Remove generated files (objs and nroff output from man pages)

clean:
	rm -f *.o ${BINS} *.nr README.html Utils/log-or-not

# Keep backup files during normal clean, but provide an option to remove them
realclean: clean
	rm -f .*.bak *.bak *.BAK *.gmon core *.core

############################################################################
# Install all target files (binaries, libraries, docs, etc.)

install: all
	${MKDIR} -p ${DESTDIR}${PREFIX}/bin ${DESTDIR}${LIBEXECDIR} \
		    ${DESTDIR}${MANDIR}/man1
	${INSTALL} -m 0755 ${BIN} ${DESTDIR}${PREFIX}/bin
	${INSTALL} -m 0755 ${LIBEXEC} ${SCRIPTS} ${DESTDIR}${LIBEXECDIR}
	${INSTALL} -m 0644 Man/*.1 ${DESTDIR}${MANDIR}/man1

install-strip: install
	${STRIP} ${DESTDIR}${PREFIX}/bin/${BIN}
	for f in ${LIBEXEC}; do \
	    ${STRIP} ${DESTDIR}${LIBEXECDIR}/$${f}; \
	done

help:
	@printf "Usage: make [VARIABLE=value ...] all\n\n"
	@printf "Some common tunable variables:\n\n"
	@printf "\tCC        [currently ${CC}]\n"
	@printf "\tCFLAGS    [currently ${CFLAGS}]\n"
	@printf "\tCXX       [currently ${CXX}]\n"
	@printf "\tCXXFLAGS  [currently ${CXXFLAGS}]\n"
	@printf "\tF77       [currently ${F77}]\n"
	@printf "\tFC        [currently ${FC}]\n"
	@printf "\tFFLAGS    [currently ${FFLAGS}]\n\n"
	@printf "View Makefile for more tunable variables.\n\n"

