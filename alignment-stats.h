#ifndef _BIOLIBC_ALIGNMENT_STATS_H_
#define _BIOLIBC_ALIGNMENT_STATS_H_

#define BL_ALIGNMENT_STATS_INIT { 0, 0 }

typedef struct
{
    unsigned long   total;
    unsigned long   overlapping;
}   bl_alignment_stats_t;

#include "alignment-stats-rvs.h"
#include "alignment-stats-accessors.h"
#include "alignment-stats-mutators.h"

#endif  // _BIOLIBC_ALIGNMENT_STATS_H_
