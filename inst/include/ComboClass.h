#ifndef COMBO_CLASS_H
#define COMBO_CLASS_H

#include "GmpDependUtils.h"
#include "ComboClassUtils.h"
#include "CombPermPtr.h"
#include "NextCombinatorics.h"
#include "PrevCombinatorics.h"

class Combo {
private:
    const int n;
    const int m;
    const int m1;
    const SEXP sexpVec;
    const bool IsFactor;
    const bool IsComb;
    const bool IsMult;
    const bool IsRep;
    const bool IsGmp;
    const double computedRows;
    mpz_t computedRowsMpz[1];
    
    double dblTemp;
    mpz_t mpzTemp;
    std::vector<int> z;
    
    const std::vector<int> freqs;
    const std::vector<int> myReps;
    
    double dblIndex;
    mpz_t mpzIndex[1];
    
    const Rcpp::IntegerVector testFactor;
    const Rcpp::CharacterVector myClass;
    const Rcpp::CharacterVector myLevels;
    
    const nthResultPtr nthResFun;
    const nextIterPtr nextIter;
    const prevIterPtr prevIter;
    
    // This has to be initialized later becuase it
    // depends on freqs.size, IsMult, and n
    const int n1;
    
    template <int RTYPE>
    Rcpp::Vector<RTYPE> VecRes(const Rcpp::Vector<RTYPE> &v);
    
    template <int RTYPE>
    Rcpp::Matrix<RTYPE> MatCombo(const Rcpp::Vector<RTYPE> &v, int nRows);
    
    template <int RTYPE>
    Rcpp::Matrix<RTYPE> MatComboReverse(const Rcpp::Vector<RTYPE> &v, int nRows);
    
    template <int RTYPE>
    Rcpp::Matrix<RTYPE> SampMatrix(const Rcpp::Vector<RTYPE> &v,
                                   const std::vector<double> &mySample,
                                   mpz_t *const myBigSamp, std::size_t sampSize);
    
public:
    
    Combo(SEXP Rv, int Rm, SEXP RcompRows,
          Rcpp::LogicalVector bVec, std::vector<std::vector<int>> freqInfo);
    
    void startOver();
    SEXP nextComb();
    SEXP prevComb();
    SEXP nextNumCombs(SEXP RNum);
    SEXP prevNumCombs(SEXP RNum);
    SEXP nextGather();
    SEXP prevGather();
    SEXP currComb();
    SEXP combIndex(SEXP RindexVec);
    SEXP front();
    SEXP back();
    SEXP sourceVector() const;
    Rcpp::List summary();
};

#endif
