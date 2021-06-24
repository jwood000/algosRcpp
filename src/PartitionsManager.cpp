#include "Partitions/PartitionsMultiset.h"
#include "Partitions/PartitionsDistinct.h"
#include "Partitions/PartitionsTypes.h"
#include "Partitions/PartitionsRep.h"
#include "RMatrix.h"

void PartsStdManager(int* mat, std::vector<int> &z, int width,
                     int lastCol, int boundary, int edge,
                     int strt, int nRows, bool IsComb, bool IsRep) {

    if (IsRep && IsComb) {
        PartsRep(mat, z, width, boundary, edge, lastCol, strt, nRows);
    } else if (IsRep) {
        PartsPermRep(mat, z, width, boundary, edge, lastCol, nRows);
    } else if (IsComb) {
        PartsDistinct(mat, z, width, boundary, lastCol, edge, strt, nRows);
    } else {
        const auto it = std::find(z.rbegin(), z.rend(), 0);

        if (it != z.rend()) {
            PartsPermZeroDistinct(mat, z, width, boundary,
                                  lastCol, edge, nRows);
        } else {
            PartsPermDistinct(mat, z, width, boundary,
                              lastCol, edge, nRows);
        }
    }
}

template <typename T>
void PartsGenManager(T* mat, const PartDesign &part, const std::vector<T> &v,
                     std::vector<int> &z, int width, int lastCol,
                     int lastElem, int boundary, int edge, int strt,
                     int nRows, bool IsComb, bool IsRep) {

    if (IsComb) {
        if (IsRep) {
            PartsGenRep(mat, v, z, width,
                        lastElem, lastCol, strt, nRows);
        } else {
            PartsGenDistinct(mat, v, z, width,
                             lastElem, lastCol, strt, nRows);
        }
    } else {
        if (part.ptype <= PartitionType::RepShort) {
            PartsGenPermRep(mat, v, z, width, lastElem, lastCol, nRows);
        } else if (part.ptype <= PartitionType::DstctCapped) {

            const auto it = std::find(z.rbegin(), z.rend(), 0);

            if (it != z.rend()) {
                PartsGenPermZeroDistinct(mat, v, z, width,
                                         lastElem, lastCol, nRows);
            } else {
                PartsGenPermDistinct(mat, v, z, width,
                                     lastElem, lastCol, nRows);
            }
        }
    }
}

template <typename T>
void PartsGenManager(std::vector<T> &partsVec, const std::vector<T> &v,
                     const std::vector<int> &Reps, std::vector<int> &z,
                     PartitionType ptype, int width, int nRows,
                     bool IsComb) {

    if (IsComb) {
        PartsGenMultiset(partsVec, v, Reps, z, width, nRows);
    } else {
        if (ptype == PartitionType::RepCapped) {
            PartsGenPermRep(partsVec, v, z, width, nRows);
        } else {
            PartsGenPermMultiset(partsVec, v, Reps, z, width, nRows);
        }
    }
}

// template <typename T>
// void ComboParallel(RcppParallel::RMatrix<T> &mat, const std::vector<T> &v,
//                    std::vector<int> &z, int n, int m, int strt, int nRows,
//                    const std::vector<int> &freqs, bool IsMult, bool IsRep) {
//
//     if (IsMult) {
//         MultisetCombination(mat, v, z, n, m, strt, nRows, freqs);
//     } else if (IsRep) {
//         CombinationsRep(mat, v, z, n, m, strt, nRows);
//     } else {
//         CombinationsDistinct(mat, v, z, n, m, strt, nRows);
//     }
// }
//

template void PartsGenManager(int*, const PartDesign&,
                              const std::vector<int>&, std::vector<int>&,
                              int, int, int, int, int, int, int, bool, bool);

template void PartsGenManager(double*, const PartDesign&,
                              const std::vector<double>&, std::vector<int>&,
                              int, int, int, int, int, int, int, bool, bool);

template void PartsGenManager(std::vector<int>&, const std::vector<int>&,
                              const std::vector<int>&, std::vector<int>&,
                              PartitionType, int, int, bool);

template void PartsGenManager(std::vector<double>&, const std::vector<double>&,
                              const std::vector<int>&, std::vector<int>&,
                              PartitionType, int, int, bool);

// template void ComboManager(Rbyte*, const std::vector<Rbyte>&,
//                            std::vector<int>&, int, int, int,
//                            const std::vector<int>&, bool, bool);
//
// template void ComboManager(Rcomplex*, const std::vector<Rcomplex>&,
//                            std::vector<int>&, int, int, int,
//                            const std::vector<int>&, bool, bool);
//
// template void ComboParallel(RcppParallel::RMatrix<int>&, const std::vector<int>&,
//                             std::vector<int>&, int, int, int, int,
//                             const std::vector<int>&, bool, bool);
//
// template void ComboParallel(RcppParallel::RMatrix<double>&, const std::vector<double>&,
//                             std::vector<int>&, int, int, int, int,
//                             const std::vector<int>&, bool, bool);