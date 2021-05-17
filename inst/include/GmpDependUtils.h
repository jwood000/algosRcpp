#ifndef GMP_DEPEND_UTILS_H
#define GMP_DEPEND_UTILS_H

#include "StandardUtils.h"
#include "StandardCount.h"
#include "BigNumCount.h"
#include "NthResult.h"
#include <cmath>

// Based off the internal limitations of sample, we
// cannot utilize the full range of 53-bit significand
// precision. Here is the condition from `do_dample2`:
//     if (!R_FINITE(dn) || dn < 0 || dn > 4.5e15 || (k > 0 && dn == 0))
// Here is the source (line ~1800): 
//     https://github.com/wch/r-source/blob/trunk/src/main/unique.c
constexpr double sampleLimit = 4500000000000000.0;

SEXP GetCount(bool IsGmp, mpz_t computedRowMpz, double computedRows);

void SetBounds(const SEXP &Rlow, const SEXP &Rhigh,bool IsGmp, bool &bLower, bool &bUpper,
               double &lower, double &upper, mpz_t *const lowerMpz,
               mpz_t *const upperMpz, mpz_t computedRowMpz, double computedRows);

void SetNumResults(bool IsGmp, bool bLower, bool bUpper, bool IsGenCnstrd, mpz_t *const upperMpz,
                   mpz_t *const lowerMpz, double lower, double upper, double computedRows,
                   mpz_t &computedRowMpz, int &nRows, double &userNumRows);

void SetRandomSampleMpz(const SEXP &RindexVec, const SEXP &RmySeed, std::size_t sampSize,
                        bool IsGmp, mpz_t &computedRowMpz, mpz_t *myVec);

void SetStartZ(int n, int m, double &lower, int stepSize, mpz_t &lowerMpz, bool IsRep,
               bool IsComb, bool IsMultiset, bool IsGmp, const std::vector<int> &myReps,
               const std::vector<int> &freqs, std::vector<int> &z, 
               const nthResultPtr nthResFun);

template <typename typeRcpp>
inline void SetSampleNames(bool IsGmp, std::size_t sampSize, typeRcpp &objRcpp,
                           const std::vector<double> &mySample, 
                           mpz_t *const myBigSamp, bool IsMatrix = true) {
    
    Rcpp::CharacterVector myRowNames(sampSize);
    
    if (IsGmp) {
        constexpr int base10 = 10;
        
        for (std::size_t i = 0; i < sampSize; ++i) {
            mpz_add_ui(myBigSamp[i], myBigSamp[i], 1);
            auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(myBigSamp[i], base10) + 2);
            mpz_get_str(buffer.get(), base10, myBigSamp[i]);
            const std::string tempRowName = buffer.get();
            myRowNames[i] = tempRowName;
        }
    } else {
        for (std::size_t i = 0; i < sampSize; ++i)
            myRowNames[i] = std::to_string(static_cast<int64_t>(mySample[i] + 1));
    }
    
    if (IsMatrix)
        Rcpp::rownames(objRcpp) = myRowNames;
    else
        objRcpp.attr("names") = myRowNames;
}

#endif
