#include "GmpDependUtils.h"
#include "CombPermPtr.h"
#include "RMatrix.h"

template <int T>
Rcpp::Matrix<T> SerialReturn(const Rcpp::Vector<T> &v, std::vector<int> &z, int n,
                             int m, int nRows, bool IsComb, bool IsRep, bool IsMult,
                             bool generalRet, const std::vector<int> &freqs) {

    Rcpp::Matrix<T> matRcpp = Rcpp::no_init_matrix(nRows, m);
    Rcpp::XPtr<combPermPtr<Rcpp::Matrix<T>, Rcpp::Vector<T>>> xpFunCoPePtr = 
        putCombPtrInXPtr<Rcpp::Matrix<T>, Rcpp::Vector<T>>(IsComb, IsMult, IsRep, generalRet);

    const combPermPtr<Rcpp::Matrix<T>, Rcpp::Vector<T>> myFunCombPerm = *xpFunCoPePtr;
    myFunCombPerm(matRcpp, v, z, n, m, 0, nRows, freqs);
    return matRcpp;
}

template <typename typeRcpp, typename T>
void MasterReturn(typeRcpp &matRcpp, std::vector<T> v, int n, int m, bool IsRep, 
                  bool IsComb, bool IsMult, bool IsGmp, bool generalRet, const std::vector<int> &freqs,
                  std::vector<int> z, const std::vector<int> &myReps, double lower, mpz_t &lowerMpz,
                  int nRows, int nThreads, bool Parallel, int phaseOne) {

    if (Parallel) {
        RcppParallel::RMatrix<T> parMat(matRcpp);

        if (generalRet) {
            RcppThread::ThreadPool pool(nThreads);
            const int stepSize = nRows / nThreads;
            int nextStep = stepSize;
            int step = 0;
            
            Rcpp::XPtr<combPermPtr<RcppParallel::RMatrix<T>, 
                                   std::vector<T>>> xpFunCoPePtr = 
                    putCombPtrInXPtr<RcppParallel::RMatrix<T>, 
                                     std::vector<T>>(IsComb, IsMult, IsRep, generalRet);
            
            const combPermPtr<RcppParallel::RMatrix<T>, 
                              std::vector<T>> myFunParCombPerm = *xpFunCoPePtr;
            
            Rcpp::XPtr<nthResultPtr> xpNthComb = putNthResPtrInXPtr(IsComb, IsMult, IsRep, IsGmp);
            const nthResultPtr nthResFun = *xpNthComb;

            for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
                pool.push(std::cref(myFunParCombPerm), std::ref(parMat),
                          std::cref(v), z, n, m, step, nextStep, std::cref(freqs));
                
                SetStartZ(n, m, lower, stepSize, lowerMpz, IsRep,
                          IsComb, IsMult, IsGmp, myReps, freqs, z, nthResFun);
            }
            
            pool.push(std::cref(myFunParCombPerm), std::ref(parMat),
                      std::cref(v), z, n, m, step, nRows, std::cref(freqs));

            pool.join();
        } else {
            PermuteParallel(parMat, v, z, n, m, nRows, phaseOne, nThreads, IsRep);
        }
    } else {
        Rcpp::XPtr<combPermPtr<typeRcpp, std::vector<T>>> xpFunCoPePtr = putCombPtrInXPtr<typeRcpp, 
                               std::vector<T>>(IsComb, IsMult, IsRep, generalRet);
        const combPermPtr<typeRcpp, std::vector<T>> myFunCombPerm = *xpFunCoPePtr;
        myFunCombPerm(matRcpp, v, z, n, m, 0, nRows, freqs);
    }
}

// [[Rcpp::export]]
SEXP CombinatoricsStndrd(SEXP Rv, SEXP Rm, SEXP RisRep, SEXP RFreqs, SEXP Rlow,
                         SEXP Rhigh, bool IsComb, bool IsFactor, SEXP Rparallel,
                         SEXP RNumThreads, int maxThreads) {

    int n, m = 0, lenFreqs = 0, nRows = 0;
    VecType myType = VecType::Integer;
    bool IsMult = false;

    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqs;

    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
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
    const bool permNonTriv = (!IsComb && bLower);
    
    Rcpp::XPtr<nthResultPtr> xpComb = putNthResPtrInXPtr(IsComb, IsMult, IsRep, IsGmp);
    const nthResultPtr nthResFun = *xpComb;
    
    SetStartZ(n, m, lower, 0, lowerMpz[0], IsRep, IsComb,
              IsMult, IsGmp, myReps, freqs, startZ, nthResFun);
    
    double userNumRows = 0;
    SetNumResults(IsGmp, bLower, bUpper, false, upperMpz.get(), lowerMpz.get(),
                  lower, upper, computedRows, computedRowsMpz, nRows, userNumRows);
    
    int nThreads = 1;
    const int limit = 20000;
    SetThreads(Parallel, maxThreads, nRows, myType, nThreads, RNumThreads, limit);

    const double phaseOneDbl = (!permNonTriv && !IsComb) ?
                               ((IsRep) ? std::pow(static_cast<double>(n),
                                             static_cast<double>(m - 1)) :
                               NumPermsNoRep(n - 1, m - 1)) : 0;
    
    bool generalRet = IsComb || IsMult || permNonTriv || 
                      n == 1 || phaseOneDbl > std::numeric_limits<int>::max();

    if (!generalRet) {
        if ((phaseOneDbl * 2.0) > nRows) {
            generalRet = true;
        } else {
            // Here, we estimate the maximum size of an array of ints
            // by taking advantage of the max_size method of vectors
            std::vector<int> sizeTestVec;
            const double first = (IsRep) ? 1 : 0;
            
            if (phaseOneDbl * (static_cast<double>(m) - first) > sizeTestVec.max_size()) {
                generalRet = true;
            }
        }
    }
    
    const int phaseOne = (generalRet) ? 0 : static_cast<int>(phaseOneDbl);

    if (myType > VecType::Logical) {
        const SEXP sexpCopy(Rcpp::clone(Rv));
        RCPP_RETURN_VECTOR(SerialReturn, sexpCopy, startZ, n, m, 
                           nRows, IsComb, IsRep, IsMult, generalRet, freqs);
    } else if (myType == VecType::Logical) {
        Rcpp::LogicalMatrix matBool = Rcpp::no_init_matrix(nRows, m);
        MasterReturn(matBool, vInt, n, m, IsRep, IsComb, IsMult,
                     IsGmp, generalRet, freqs, startZ, myReps, lower,
                     lowerMpz[0], nRows, nThreads, Parallel, phaseOne);
        return matBool;
    } else if (myType == VecType::Integer) {
        Rcpp::IntegerMatrix matInt = Rcpp::no_init_matrix(nRows, m);
        MasterReturn(matInt, vInt, n, m, IsRep, IsComb, IsMult,
                     IsGmp, generalRet, freqs, startZ, myReps, lower,
                     lowerMpz[0], nRows, nThreads, Parallel, phaseOne);
        if (IsFactor) {SetFactorClass(matInt, Rv);}
        return matInt;
    } else {
        Rcpp::NumericMatrix matNum = Rcpp::no_init_matrix(nRows, m);
        MasterReturn(matNum, vNum, n, m, IsRep, IsComb, IsMult,
                     IsGmp, generalRet, freqs, startZ, myReps, lower,
                     lowerMpz[0], nRows, nThreads, Parallel, phaseOne);
        return matNum;
    }
}

