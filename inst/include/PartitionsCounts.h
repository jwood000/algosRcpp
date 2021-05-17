#ifndef PARTITIONS_COUNTS_H
#define PARTITIONS_COUNTS_H

#include "CleanConvert.h"

double GetComputedPartsComps(const std::vector<int> &z, PartitionType PartType, 
                             int target, int m, bool IsComb, bool IncludeZero, bool mIsNull);
    
#endif
