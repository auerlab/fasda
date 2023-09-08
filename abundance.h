#ifndef __ABUNDANCE_H__

#define MAX_CMD_LEN         8192
#define MAX_FILE_COUNT      1024    // SAM stream args to abundance.c
#define MAX_SEQ_LEN         1024
#define MAX_GTF_LINE_LEN    8192

#define GFF3_MASK        BL_GFF3_FIELD_SEQID|BL_GFF3_FIELD_TYPE|\
			BL_GFF3_FIELD_START|BL_GFF3_FIELD_END|\
			BL_GFF3_FIELD_ATTRIBUTES
#define SAM_MASK        BL_SAM_FIELD_RNAME|BL_SAM_FIELD_POS
// -@ 2 does not help, at least with BAM files.  CRAM files are
// more expensive to decode, so it may help there.
#define SAMTOOLS_ARGS   "--exclude-flags UNMAP"

#define FASDA_FLAG_SHOW_GENE    0x01
#define FASDA_IGNORE_CHR_ORDER  0x02

typedef enum { STRINGTIE, EXACT } abundance_method_t;

#ifndef _BIOLIBC_ALIGNMENT_STATS_H_
#include "alignment-stats.h"
#endif

#include "abundance-protos.h"

#endif
