#include "GmpDependUtils.h"
#include "CombinationApply.h"
#include "PermutationsApply.h"

template <int RTYPE>
Rcpp::List ApplyFunction(const Rcpp::Vector<RTYPE> &v, int n, int m, bool IsRep, int nRows,
                         bool IsComb, const std::vector<int> &freqs, std::vector<int> &z, 
                         bool IsMult, SEXP stdFun, SEXP rho) {
    
    Rcpp::List myList(nRows);
    SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));
    
    if (IsComb) {
        if (IsMult)
            MultisetComboApplyFun(myList, v, z, n, m, nRows, sexpFun, rho, freqs);
        else
            ComboGeneralApplyFun(myList, v, z, n, m, IsRep, nRows, sexpFun, rho);
    } else {
        PermutationApplyFun(myList, v, z, n, m, IsRep, IsMult, nRows, sexpFun, rho);
    }
    
    UNPROTECT(1);
    return myList;
}

// [[Rcpp::export]]
SEXP CombinatoricsApply(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow,
                        SEXP Rhigh, bool IsComb, SEXP stdFun, SEXP myEnv) {

    int n, m = 0, lenFreqs = 0, nRows = 0;
    VecType myType = VecType::Integer;
    bool IsMult = false;

    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqs;
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");

    SetType(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    SetFreqsAndM(RFreqs, IsMult, myReps, IsRep, lenFreqs, freqs, Rm, n, m);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep, n,
                                                m, Rm, lenFreqs, freqs, myReps);

    const bool IsGmp = (computedRows > Significand53);
    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult, 
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }
    
    double lower = 0, upper = 0;
    bool bLower = false, bUpper = false;
    
    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);
    mpz_init(lowerMpz[0]); mpz_init(upperMpz[0]);

    SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
              lowerMpz.get(), upperMpz.get(), computedRowsMpz, computedRows);

    std::vector<int> startZ(m);
    Rcpp::XPtr<nthResultPtr> xpComb = putNthResPtrInXPtr(IsComb, IsMult, IsRep, IsGmp);
    const nthResultPtr nthResFun = *xpComb;
    
    SetStartZ(n, m, lower, 0, lowerMpz[0], IsRep, IsComb,
              IsMult, IsGmp, myReps, freqs, startZ, nthResFun);
    
    double userNumRows = 0;
    SetNumResults(IsGmp, bLower, bUpper, false, upperMpz.get(), lowerMpz.get(),
                  lower, upper, computedRows, computedRowsMpz, nRows, userNumRows);
    
    const SEXP sexpCopy = CopyRv(Rv, vInt, vNum, myType);
    RCPP_RETURN_VECTOR(ApplyFunction, sexpCopy, n, m, IsRep, nRows,
                       IsComb, freqs, startZ, IsMult, stdFun, myEnv);
}

