#include "StandardCount.h"
#include "NextPartitions.h"
#include "Cpp14MakeUnique.h"

void BinaryNextElem(int &uppBnd, int &lowBnd, int &ind, int lastElem,
                    int64_t target, int64_t partial, const std::vector<int64_t> &v) {
    
    int64_t dist = target - (partial + v[ind]);
    
    while ((uppBnd - lowBnd) > 1 && dist != 0) {
        const int mid = (uppBnd - lowBnd) / 2;
        ind = lowBnd + mid;
        dist = target - (partial + v[ind]);
        
        if (dist > 0)
            lowBnd = ind;
        else
            uppBnd = ind;
    }
    
    // Check last index. N.B. There are some cases when ind == lowBnd and dist < 0.
    // This will not matter as we simply reassign ind and recompute dist 
    if (dist < 0) {
        ind = lowBnd;
        dist = target - (partial + v[ind]);
    }
    
    // We must have dist < 0. Below is an informal proof.
    // The sub-sequences are defined as below:
    //                  A_max = a_(i + 1), a_(i + 2), ..., a_m
    //                  A_set = a_1, a_2, ..., a_(i - 1)
    // A_set are those elements that have already been determined by the algorithm.
    // A_max is maximal (i.e. constructed of the the (m - i) largest elements). We
    // seek to determine the i_th element given the following contraints:
    //                      A_sum = A_set + a_i + A_max
    //                       dist = target - A_sum
    // With the goal of finding the minimum lexicographic combination such that the
    // dist = 0 (i.e. target = A_sum). If we have dist > 0 for any i, then it will
    // be impossible to obtain dist = 0. dist > 0 implies that the target > A_sum,
    // and since A_max is already maximal, we are not able to increase A_sum in
    // later iterations, thus we must have dist <= 0 for all i.
    if (dist > 0 && ind < lastElem)
        ++ind;
}

int GetFirstCombo(int m, const std::vector<int64_t> &v, bool IsRep, bool IsMult,
                  std::vector<int> &z, const std::vector<int> &freqs, int64_t target,
                  const std::vector<int> &Reps, int lastCol, int lastElem) {
    
    int64_t testMax = 0;
    constexpr int64_t zero64 = 0;
    
    if (IsRep) {
        testMax = v[lastElem] * m;
    } else if (IsMult) {
        const int lenMinusM = freqs.size() - m;
        
        for (int i = freqs.size() - 1, j = 0; i >= lenMinusM; --i, ++j)
            testMax += v[freqs[i]];
    } else {
        testMax = std::accumulate(v.cend() - m, v.cend(), zero64);
    }
    
    if (testMax < target)  {return 0;}
    int zExpCurrPos = IsMult ? freqs.size() - m : 0;
    int currPos = IsMult ? freqs[zExpCurrPos] : (IsRep ? lastElem : (v.size() - m));
    
    int64_t partial = testMax;
    partial -= v[currPos];
    int64_t testMin = 0;
    
    if (IsRep) {
        testMin = v[0] * m;
    } else if (IsMult) {
        for (int i = 0; i < m; ++i)
            testMin += v[freqs[i]];
    } else {
        testMin = std::accumulate(v.cbegin(), v.cbegin() + m, zero64);
    }
    
    if (testMin > target)  {return 0;}
    int mid = currPos / 2;
    int64_t dist = target - (partial + v[mid]);
    
    int lowBnd = (dist > 0) ? mid : 0;
    int uppBnd = (dist > 0) ? currPos : mid;
    int ind = mid;
    
    std::vector<int> repsCounter;
    
    if (IsMult)
        repsCounter.assign(Reps.cbegin(), Reps.cend());
    
    for (int i = 0; i < lastCol; ++i) {
        BinaryNextElem(uppBnd, lowBnd, ind, lastElem, target, partial, v);
        z[i] = ind;
        partial += v[ind];
        
        if (IsMult) {
            --repsCounter[ind];
            
            if (repsCounter[ind] == 0)
                ++ind;
            
            ++zExpCurrPos;
            currPos = freqs[zExpCurrPos];
        } else if (!IsRep) {
            ++ind;
            ++currPos;
        }
        
        lowBnd = ind;
        uppBnd = currPos;
        mid = (uppBnd - lowBnd) / 2;
        
        ind = lowBnd + mid;
        partial -= v[currPos];
    }
    
    BinaryNextElem(uppBnd, lowBnd, ind, lastElem, target, partial, v);
    z[lastCol] = ind;
    
    // The algorithm above finds the first possible sum that equals
    // target. If there is no combination of elements from v that sum
    // to target, the algo returns the combination such that its sum
    // is closest to target and greater than target
    int64_t finalCheck = 0;
    
    for (int i = 0; i < m; ++i)
        finalCheck += v[z[i]];
    
    if (finalCheck != target)
        return 0;
    
    return 1;
}

inline void PopulateVec(int m, const std::vector<int64_t> &v, std::vector<int> &z,
                        int &count, int maxRows, bool isComb, std::vector<int64_t> &partitionsVec) {
    
    if (isComb) {
        for (int k = 0; k < m; ++k)
            partitionsVec.push_back(v[z[k]]);
        
        ++count;
    } else {
        do {
            for (int k = 0; k < m; ++k)
                partitionsVec.push_back(v[z[k]]);
            
            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
    }
}

inline bool keepGoing(const std::vector<int> &rpsCnt, int lastElem,
                      const std::vector<int> &z, int edge, int boundary) {
    
    if (edge >= 0) {
        const int myDiff = z[boundary] - z[edge];
        
        if (myDiff < 2) {
            return false;
        } else if (myDiff == 2) {
            return (rpsCnt[z[edge] + 1] > 1);
        } else {
            return (rpsCnt[z[edge] + 1] && rpsCnt[z[boundary] - 1]);
        }
    } else {
        return false;
    }
}

int PartitionsMultiSet(int m, const std::vector<int64_t> &v, int64_t target, 
                       int lastElem, int lastCol, int maxRows, bool isComb,
                       const std::vector<int> &Reps, std::vector<int64_t> &partitionsVec) {
    
    std::vector<int> z(m);
    std::vector<int> zExpanded;
    
    for (std::size_t i = 0; i < v.size(); ++i)
        for (int j = 0; j < Reps[i]; ++j)
            zExpanded.push_back(i);
    
    if (!GetFirstCombo(m, v, false, true, z, zExpanded, target, Reps, lastCol, lastElem))
        return 0;
    
    std::vector<int> rpsCnt(Reps.cbegin(), Reps.cend());
    
    for (const auto ind: z)
        --rpsCnt[ind];
    
    // boundary is the greatest index that can be decremented
    // pivot is the greatest index that can be incremented
    // edge is the greatest index smaller than boundary that can be incremented
    int p, e, b = lastCol, count = 0;
    PrepareMultiSetPart(rpsCnt, z, b, p, e, lastCol, lastElem);
    
    while (keepGoing(rpsCnt, lastElem, z, e, b)) {
        PopulateVec(m, v, z, count, maxRows, isComb, partitionsVec);
        
        if (count >= maxRows)
            break;
        
        NextMultiSetGenPart(rpsCnt, z, e, b,  p, lastCol, lastElem);
    }
    
    
    if (count < maxRows)
        PopulateVec(m, v, z, count, maxRows, isComb, partitionsVec);
    
    return 1;
}

int PartitionsRep(int m, const std::vector<int64_t> &v, int64_t target, int lastElem,
                  int lastCol, int maxRows, bool isComb, std::vector<int64_t> &partitionsVec) {
    
    std::vector<int> z(m);
    std::vector<int> trivVec;
    
    if (!GetFirstCombo(m, v, true, false, z, trivVec, target, trivVec, lastCol, lastElem))
        return 0;
    
    // smallest index such that z[boundary] == currMax
    int boundary = lastCol;
    
    while (boundary > 1 && z[boundary - 1] == z[lastCol])
        --boundary;
    
    // pivot is the greatest index that can be incremented.
    // We know that if z[boundary] < lastElem ==>> pivot = lastCol
    int pivot = (z[boundary] < lastElem) ? lastCol : boundary - 1;
    
    // edge is the greatest index such that z[boundary] - z[edge] >= 2
    // This is the index that will be be used as a starting point
    // to determine the next combination that meets the criteria
    int edge = boundary - 1;
    int edgeTest = z[boundary] - 2;
    
    while (edge && edgeTest < z[edge])
        --edge;
    
    int count = 0;
    
    while ((edge >= 0) && (z[boundary] - z[edge] >= 2)) {
        PopulateVec(m, v, z, count, maxRows, isComb, partitionsVec);
        
        if (count >= maxRows)
            break;
        
        NextRepGenPart(z, boundary, edge, pivot, lastCol, lastElem);
    }
    
    if (count < maxRows)
        PopulateVec(m, v, z, count, maxRows, isComb, partitionsVec);
    
    return 1;
}

int PartitionsDistinct(int m, const std::vector<int64_t> &v, int64_t target, int lastElem,
                       int lastCol, int maxRows, bool isComb, std::vector<int64_t> &partitionsVec) {
    
    std::vector<int> z(m);
    std::vector<int> trivVec;
    
    if (!GetFirstCombo(m, v, false, false, z, trivVec, target, trivVec, lastCol, lastElem))
        return 0;
    
    // Largest index such that z[boundary] - z[boundary - 1] > 1
    // This is the index that will be decremented at the same
    // time edge is incremented. See below.
    int boundary = lastCol;
    
    while (boundary > 1 && (z[boundary] - z[boundary - 1]) < 2)
        --boundary;
    
    // pivot is the greatest index that can be incremented...
    // Either z[pivot + 1] - z[pivot] > 1 or if z[lastCol] < lastElem
    // pivot = lastCol since incrementing z[lastCol] is possible
    int pivot = (z[lastCol] < lastElem) ? lastCol : boundary - 1;
    
    // edge is the greatest index such that when incremented
    // the result will be at least one less than its neighbor
    // even if its neighbor is decremented
    int edge = boundary - 1;
    int tarDiff = 3;
    
    while (edge && (z[boundary] - z[edge]) < tarDiff) {
        --edge;
        ++tarDiff;
    }
    
    int count = 0;
    
    while (edge >= 0 && (z[boundary] - z[edge]) >= tarDiff) {
        PopulateVec(m, v, z, count, maxRows, isComb, partitionsVec);
        
        if (count >= maxRows)
            break;
        
        NextDistinctGenPart(z, boundary, edge, pivot, tarDiff, lastCol, lastElem);
    }
    
    if (count < maxRows)
        PopulateVec(m, v, z, count, maxRows, isComb, partitionsVec);
    
    return 1;
}
