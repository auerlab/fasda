/***************************************************************************
 *  Description:
 *      Wrapper to turn FASDA commands into subcommands.  This will help
 *      avoid future conflicts with other programs without sacrificing
 *      desriptive command names.
 *
 *      This wrapper can be installed under the standard
 *      PATH and used to to execute fasda commands installed under a
 *      private prefix, without altering PATH, activating a special
 *      environment, opening a container, etc.  This sub-command paradigm
 *      is already familiar to bioinformaticians thanks to other suites
 *      like samtools, bedtools, etc.
 *
 *      The standard location for executables not meant to be in PATH
 *      is $PREFIX/libexec.
 *
 *      Example:
 *
 *          fasda fastx-derep args
 *
 *      instead of one of the following:
 *
 *          prefix/bin/fastx-derep args
 *
 *          env PATH=prefix/bin:$PATH fastx-derep args
 *
 *  Arguments:
 *      The subcommand and its specific arguments, as if it were run
 *      directly.
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-09-13  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <sysexits.h>
#include <limits.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

int     main(int argc,char *argv[])

{
    char    cmd[PATH_MAX + 1];
    DIR     *dp;
    struct dirent   *dir_entry;
    struct stat     inode;
    
    if ( (argc == 2) && (strcmp(argv[1],"--version") == 0) )
    {
	printf("%s %s\n", argv[0], VERSION);
	return EX_OK;
    }
    else if ( argc < 2 )
    {
	// LIBEXECDIR must be set by Makefile
	fprintf(stderr, "Usage: %s --version\n", argv[0]);
	fprintf(stderr, "       %s subcommand [args]\n", argv[0]);
	fprintf(stderr, "\nSubcommands:\n\n");
	if ( (dp = opendir(LIBEXECDIR)) != NULL )
	{
	    while ( (dir_entry = readdir(dp)) != NULL )
	    {
		if ( (*dir_entry->d_name != '.') &&
		     (strstr(dir_entry->d_name, ".awk") == NULL) )
		fprintf(stderr, "%s\n", dir_entry->d_name);
	    }
	    closedir(dp);
	}
	fprintf(stderr, "\nRun \"fasda subcommand\" or \"man fasda-subcommand\" for more details.\n\n");
	return EX_USAGE;
    }

    fputs("\n********************************* Note *********************************\n"
	  "FASDA is beta-quality software.  Be sure to verify all results.\n"
	  "(This is not to imply that other DE tools should be trusted, either.)\n"
	  "Please contribute by reporting problems and offering suggestions at\n"
	  "https://github.com/auerlab/fasda.\n"
	  "************************************************************************\n\n",
	  stderr);
    snprintf(cmd, PATH_MAX, "%s/%s", LIBEXECDIR, argv[1]);
    if ( stat(cmd, &inode) == 0 )
	execv(cmd, argv + 1);
    else
    {
	fprintf(stderr, "%s: No %s found in %s.\n", argv[0], argv[1], LIBEXECDIR);
	return EX_USAGE;
    }
}
