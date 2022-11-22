#ifndef __ABUNDANCE_H__

#define MAX_FILE_COUNT  1024    // SAM stream args to abundance.c
#define MAX_SEQ_LEN     1024
#define GFF_MASK        BL_GFF_FIELD_SEQID|BL_GFF_FIELD_TYPE|\
			BL_GFF_FIELD_START|BL_GFF_FIELD_END|\
			BL_GFF_FIELD_ATTRIBUTES
#define SAM_MASK        BL_SAM_FIELD_RNAME|BL_SAM_FIELD_POS
// -@ 2 does not help, at least with BAM files.  CRAM files are
// more expensive to decode, so it may help there.
#define SAMTOOLS_ARGS   "--exclude-flags UNMAP"

#define FASDA_FLAG_SHOW_GENE 0x01
#define FASDA_FLAG_MAP_GENE  0x02

#ifndef _BIOLIBC_ALIGNMENT_STATS_H_
#include "alignment-stats.h"
#endif

#include "abundance-protos.h"

#endif
