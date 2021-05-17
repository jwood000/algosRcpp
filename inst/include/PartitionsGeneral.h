#ifndef PARTITIONS_GENERAL_H
#define PARTITIONS_GENERAL_H

#include "CleanConvert.h"
#include <vector>
#include <cstdint>

int PartitionsMultiSet(int m, const std::vector<std::int64_t> &v, std::int64_t target, 
                       int lastElem, int lastCol, int maxRows, bool isComb,
                       const std::vector<int> &Reps, std::vector<std::int64_t> &partitionsVec);
    
int PartitionsRep(int m, const std::vector<std::int64_t> &v, std::int64_t target, int lastElem,
                  int lastCol, int maxRows, bool isComb, std::vector<std::int64_t> &partitionsVec);

int PartitionsDistinct(int m, const std::vector<std::int64_t> &v, std::int64_t target, int lastElem,
                       int lastCol, int maxRows, bool isComb, std::vector<std::int64_t> &partitionsVec);

#endif
