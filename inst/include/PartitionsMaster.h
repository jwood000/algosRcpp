#ifndef PARTITIONS_MASTER_H
#define PARTITIONS_MASTER_H

#include "CleanConvert.h"

namespace Partitions {
    template <typename typeRcpp>
    typeRcpp PartitionsMaster(const std::vector<int64_t> &v, std::vector<int> &z, const std::vector<int> &Reps,
                              PartitionType PartType, int64_t target, int lenV, int m, bool IsRep, bool IsMult,
                              double numParts, bool IsComb, bool xtraCol, bool bUserRows, bool mIsNull, bool getAll);
}

#endif
