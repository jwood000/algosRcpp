#ifndef COMB_PERM_PTR_H
#define COMB_PERM_PTR_H

#include "Combinations.h"
#include "Permutations.h"

template <typename T, typename U>
using combPermPtr = void (*const)(T &matRcpp, const U &v, std::vector<int> z,
                          int n, int m, int strt, int nRows, const std::vector<int> &freqs);

template <typename T, typename U>
Rcpp::XPtr<combPermPtr<T, U>> putCombPtrInXPtr(bool IsComb, bool IsMult, bool IsRep, bool IsGen) {
    
    if (IsComb) {
        if (IsMult)
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&MultisetCombination)));
        else if (IsRep)
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&CombinationsRep)));
        else
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&CombinationsNoRep)));
    } else {
        if (IsMult) {
            return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&MultisetPermutation)));
        } else if (IsGen) {
            if (IsRep)
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteGeneralRep)));
            else
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteGeneralNoRep)));
        } else {
            if (IsRep)
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteSerialRep)));
            else
                return(Rcpp::XPtr<combPermPtr<T, U>>(new combPermPtr<T, U>(&PermuteSerialNoRep)));
        }
    }
}

#endif
