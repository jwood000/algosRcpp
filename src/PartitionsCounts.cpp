#include "StandardCount.h"
#include "CleanConvert.h"

// Credit to Robin K. S. Hankin, author of the excellent partitions package.
// From the partitions.c, here are Hankin's comments for c_numbdiffparts:
//      "the recursion on p826 of Abramowitz and Stegun"
double c_numbdiffparts(int n) {
    
    const int n1 = n + 1;
    std::vector<double> qq(n1, 1);
    
    for(int i = 2 ; i < n1; ++i) {
        qq[i] = 0;
        
        for (int s = 1, f = 5, r = 2; i >= r; r += f, f += 3, s *= -1) {
            qq[i] += s * qq[i - r];
            if(i == r * 2) {qq[i] -= s;}
        }
        
        for (int s = 1, f = 4, r = 1; i >= r; r += f, f += 3, s *= -1) {
            qq[i] += s * qq[i - r];
            if(i == r * 2) {qq[i] -= s;}
        }
    }
    
    return qq.back();
}

// Consider finding partitions of 50 of length 5 with no zeros.
// The first few results are:
//                  1 1 1 1 46
//                  1 1 1 2 45
//                  1 1 1 3 44
//
// If we cont. until the num in the 3rd col changes, we have:
//              [22,] 1 1 1 22 25
//              [23,] 1 1 1 23 24
//              [24,] 1 1 2  2 44
//
// That's 23 results, which is expected as we are decrementing
// the right most number and incrementing the number in the
// 4th column. Thus we have:
//                   (46 - n) > (1 + n)
//              ==>        45 > 2n
//              ==>      22.5 > n
//              ==>         n = 22
//
// Including the first result (i.e. 1 46), we have n + 1 = 23.
// In general, in can be show that when we are presented with
// a situaiton of the form: 1 1 ... 1 m, then the number of
// results for the "1st" group (i.e. while the 2nd to the last
// column stays the same), is given by floor((1 + m) / 2).
//
// Continuing on with the next group, we have:
//                   1 1 2 2 44
//
// We need to find the number of result of (2 44). This is
// equivalent to (1 43). The next few groups are (3 42), (4, 40),
// (5, 38), and (6, 36). These all get converted to (1, 40),
// (1, 37), (1, 34), and (1, 31). The number of results for
// each are given by (assume integer division and don't 
// forget to add 1):
//          44 / 2, 41 / 2, 38 / 2, 35 / 2, 32 / 2
//
// The total number of result for all groups is given by:
//          47 / 2 + 44 / 2 + ... + 5 / 2 + 2 / 2
//
// This can be converted to:
// {(47 - 0) + (47 - 3) + ... + (47 - 42) + (47 - 45)} / 2
//
// The 8 below comes from adding 1 to every other term
// that needs to be removed (leave as a proof to the reader)
//   ==>>  (47 * 16 - 3 * (0 + 1 + 2 + .... + 14 + 15) - 8) / 2
//   ==>>    (47 * 16 - 3 * (16 * 15 / 2) - 8) / 2

double SumSection(int n) {
    const int numIter = n / 3;
    const double sumOne = numIter * (n - 2);
    const double sumTwo = 3 * (numIter * (numIter - 1)) / 2;
    return std::floor((sumOne - sumTwo - numIter / 2) / 2);
}

// The algos below are modified versions of that found in:
// (Credit to user @m69) https://stackoverflow.com/a/32918426/4408538
//
// These are faster as they skip many iterations as a result of
// returning a result when m == 3 as opposed to m == 2. There are
// more opportunities for improvements such as deriving a constant
// time calculation when m == 4 (Currently under development).
// *********************************************************************** //
static std::vector<double> memoize;

double pStd(int n, int m, int w) {
    
    if (memoize[(n - m) * w + m - 2])
        return memoize[(n - m) * w + m - 2];
    
    if (m == 3)
        return SumSection(n + 1);
    
    int myLim = n / m;
    double count = 0;
    
    for (; myLim--; n -= m)
        count += (memoize[(n - m) * w + m - 3] = pStd(n - 1, m - 1, w));
    
    return count;
}

double CountStdPartLen(int n, int m) {
    
    if (m < 2)
        return 1;
    
    if (m == 2)
        return std::floor(n / 2);
    
    memoize.assign((n - m + 1) * m, 0.0);
    double count = pStd(n, m, m);
    
    return count;
}

// *********************************************************************** //

double pDist(int n, int m, int w) {
    
    if (memoize[n * w + m]) 
        return memoize[n * w + m];
    
    if (m == 3)
        return SumSection(n - 2);
    
    int myLim = n / m;
    double count = 0;
    
    for (; --myLim; n -= m)
        count += (memoize[(n - m) * w + m - 1] = pDist(n - m, m - 1, w));
    
    return count;
}

double CountDistinctPartLen(int n, int m) {
    
    if (m < 2)
        return 1;
    
    if (m == 2)
        return std::floor((n - 1) / 2);
    
    memoize.assign((n + 1) * (m + 1), 0.0);
    double count = pDist(n, m, m + 1);
    
    return count;
}

double GetComputedPartsComps(const std::vector<int> &z, PartitionType PartType, 
                             int target, int m, bool IsComb, bool IncludeZero, bool mIsNull) {
    
    double partCountTest = 0;
    const bool IsDistinct = PartType > PartitionType::PartTradNoZero;
    const int startLen = std::count_if(z.cbegin(), z.cend(), 
                                       [](int i){return i > 0;});
    
    if (IsDistinct) {
        if (IsComb) {
            if (PartType == PartitionType::PartDstctStdAll) {
                partCountTest = c_numbdiffparts(target);
            } else {
                for (int i = startLen; i <= m; ++i)
                    partCountTest += CountDistinctPartLen(target, i);
            }
        } else {
            if (IncludeZero) {
                // Given z = c(0, 0, 3, 4)
                // When m is Null, here is the output (only permute last
                // 2 indices):
                //                  c(0, 0, 3, 4)
                //                  c(0, 0, 4, 3)
                //
                // And when m is give, we have (permute all indices):
                //                   [,1] [,2] [,3] [,4]
                //              [1,]    0    0    3    4
                //              [2,]    0    0    4    3
                //              [3,]    0    3    0    4
                //                 .    .    .    .    .
                //             [10,]    4    0    0    3
                //             [11,]    4    0    3    0
                //             [12,]    4    3    0    0
                
                if (mIsNull) {
                    
                    for (int i = startLen; i <= m; ++i)
                        partCountTest += (CountDistinctPartLen(target, i) * NumPermsNoRep(i, i));
                    
                } else {
                    std::vector<int> permCountVec(m);
                    std::iota(permCountVec.begin(), permCountVec.begin() + startLen, 1);
                    
                    for (int i = startLen; i <= m; ++i) {
                        permCountVec[i - 1] = i;
                        partCountTest += (CountDistinctPartLen(target, i)
                                              * NumPermsWithRep(permCountVec));
                    }
                }
            } else {
                partCountTest = CountDistinctPartLen(target, m) * NumPermsNoRep(m, m);
            }
        }
    } else {
        if (IsComb) {
            for (int i = startLen; i <= m; ++i)
                partCountTest += CountStdPartLen(target, i);
            
        } else if (mIsNull && IncludeZero) {
            partCountTest = std::pow(2.0, static_cast<double>(target - 1));
        } else {
            partCountTest = (IncludeZero) ? nChooseK(target + m - 1, m - 1)
                                          : nChooseK(target - 1, m - 1);
        }
    }
    
    return partCountTest;
}
