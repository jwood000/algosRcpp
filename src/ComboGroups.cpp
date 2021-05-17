#include "ComboGroupsUtils.h"
#include "RMatrix.h"
#include <RcppThread/ThreadPool.hpp>

template <typename typeRcpp>
void FinalTouch(typeRcpp &GroupsMat, bool IsArray, int grpSize, int r, int n, 
                int nRows, bool IsNamed, const std::vector<double> &mySample, 
                mpz_t *const myBigSamp, bool IsGmp) {
    
    std::vector<std::string> myColNames(r, "Grp");
    
    for (int j = 0; j < r; ++j)
        myColNames[j] += std::to_string(j + 1);
    
    if (IsArray) {
        Rcpp::CharacterVector rcppCols(r);
        rcppCols = myColNames;
        GroupsMat.attr("dim") = Rcpp::IntegerVector::create(nRows, grpSize, r);
        GroupsMat.attr("dimnames") = Rcpp::List::create(R_NilValue, R_NilValue, rcppCols);
    } else {
        std::vector<std::string> extendedName;
        
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < grpSize; ++j)
                extendedName.push_back(myColNames[i]);
        
        Rcpp::CharacterVector rcppExtended(n);
        rcppExtended = extendedName;
        Rcpp::colnames(GroupsMat) = rcppExtended; 
    }
    
    if (IsNamed)
        SetSampleNames(IsGmp, nRows, GroupsMat, mySample, myBigSamp);
}

template <typename typeRcpp, typename typeVector>
void SampleWorkerGmp(std::size_t n, const typeVector &v, typeRcpp &GroupsMat,
                     int r, int grpSize, std::size_t strtIdx, std::size_t endIdx,
                     mpz_t *const myBigSamp, mpz_t &computedRowMpz) {
    
    for (std::size_t i = strtIdx; i < endIdx; ++i) {
        const std::vector<int> z = nthComboGroupGmp(n, grpSize, r, myBigSamp[i], computedRowMpz);
        
        for (std::size_t j = 0; j < n; ++j)
            GroupsMat(i, j) = v[z[j]];
    }
}

template <typename typeRcpp, typename typeVector>
void SampleWorker(std::size_t n, const typeVector &v, typeRcpp &GroupsMat,
                  int r, int grpSize, std::size_t strtIdx, std::size_t endIdx,
                  const std::vector<double> &mySample, double computedRows) {

    for (std::size_t i = strtIdx; i < endIdx; ++i) {
        const std::vector<int> z = nthComboGroup(n, grpSize, r, mySample[i], computedRows);
        
        for (std::size_t j = 0; j < n; ++j)
            GroupsMat(i, j) = v[z[j]];
    }
}

template <typename typeRcpp, typename typeVector>
void GroupWorker(std::size_t n, const typeVector &v, typeRcpp &GroupsMat, std::vector<int> z,
                 int r, int grpSize, std::size_t strtIdx, std::size_t endIdx) {
    
    const int idx1 = (r - 1) * grpSize - 1;
    const int idx2 = v.size() - 1;
    const int last1 = (r - 2) * grpSize + 1;
    const std::size_t lastRow = endIdx - 1;
    
    for (std::size_t i = strtIdx; i < lastRow; ++i) {
        for (std::size_t j = 0; j < n; ++j)
            GroupsMat(i, j) = v[z[j]];
        
        nextComboGroup(z, r, grpSize, idx1, idx2, last1);
    }
    
    // Get last combo group
    for (std::size_t j = 0; j < n; ++j)
        GroupsMat(lastRow, j) = v[z[j]];
}

template <typename typeRcpp, typename typeVector>
void SerialGlue(std::size_t n, const typeVector &v, typeRcpp &GroupsMat, bool IsGmp,
                const std::vector<int> &z, int r, int grpSize, std::size_t nRows,
                bool IsSample, const std::vector<double> &mySample, mpz_t *const myBigSamp,
                double computedRows, mpz_t &computedRowMpz, bool IsArray, bool IsNamed) {
    
    if (IsSample) {
        if (IsGmp) {
            SampleWorkerGmp(n, v, GroupsMat, r, grpSize, 0, nRows, myBigSamp, computedRowMpz);
        } else {
            SampleWorker(n, v, GroupsMat, r, grpSize, 0, nRows, mySample, computedRows);
        }
    } else {
        GroupWorker(n, v, GroupsMat, z, r, grpSize, 0, nRows);
    }
    
    FinalTouch(GroupsMat, IsArray, grpSize, r, n,
               nRows, IsNamed, mySample, myBigSamp, IsGmp);
}

template <int T>
Rcpp::Matrix<T> NonThreadSafe(const Rcpp::Vector<T> &v, const std::vector<int> &z, int n,
                              bool IsGmp, int r, int grpSize, int nRows, bool IsSample,
                              const std::vector<double> &mySample, mpz_t *const myBigSamp,
                              double computedRows, mpz_t &computedRowMpz, bool IsArray, bool IsNamed) {
    
    Rcpp::Matrix<T> matRcpp = Rcpp::no_init_matrix(nRows, n);
    SerialGlue(n, v, matRcpp, IsGmp, z, r, grpSize, nRows, IsSample, mySample, 
               myBigSamp, computedRows, computedRowMpz, IsArray, IsNamed);
    return matRcpp;
}

template <typename typeRcpp, typename typeElem>
void ParallelGlue(std::size_t n, const std::vector<typeElem> &v, typeRcpp &GroupsMat, bool IsGmp,
                  const std::vector<int> &z, int r, int grpSize, std::size_t strtIdx, std::size_t endIdx,
                  bool IsSample, const std::vector<double> &mySample, mpz_t *const myBigSamp,
                  double computedRows, mpz_t &computedRowMpz) {
    
    if (IsSample) {
        if (IsGmp) {
            SampleWorkerGmp(n, v, GroupsMat, r, grpSize, strtIdx, endIdx, myBigSamp, computedRowMpz);
        } else {
            SampleWorker(n, v, GroupsMat, r, grpSize, strtIdx, endIdx, mySample, computedRows);
        }
    } else {
        GroupWorker(n, v, GroupsMat, z, r, grpSize, strtIdx, endIdx);
    }
}

void GetStartGrp(bool IsGmp, std::vector<int> &z, int n, int grpSize, 
                 int r, double &lower, double computedRows, mpz_t lowerMpz,
                 mpz_t computedRowMpz, std::size_t stepSize) {
    
    if (IsGmp) {
        mpz_add_ui(lowerMpz, lowerMpz, stepSize);
        z = nthComboGroupGmp(n, grpSize, r, lowerMpz, computedRowMpz);
    } else {
        lower += stepSize;
        z = nthComboGroup(n, grpSize, r, lower, computedRows);
    }
}

template <typename typeRcpp, typename typeElem>
void GroupsMaster(std::size_t n, const std::vector<typeElem> &v, typeRcpp &GroupsMat, std::vector<int> z,
                  int r, int grpSize, std::size_t nRows, bool Parallel, int nThreads, bool IsGmp,
                  double lower, mpz_t &lowerMpz, double computedRows, mpz_t &computedRowMpz, 
                  bool IsSample, const std::vector<double> &mySample, mpz_t *const myBigSamp,
                  bool IsArray, bool IsNamed) {
    
    if (Parallel) {
        RcppParallel::RMatrix<typeElem> parMat(GroupsMat);
        RcppThread::ThreadPool pool(nThreads);
        std::size_t step = 0, stepSize = nRows / nThreads;
        std::size_t nextStep = stepSize;
        
        for (int j = 0; j < (nThreads - 1); ++j, step += stepSize, nextStep += stepSize) {
            pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), n,
                      std::cref(v), std::ref(parMat), IsGmp, z, r, grpSize, step, nextStep, IsSample, 
                      std::cref(mySample), myBigSamp, computedRows, std::ref(computedRowMpz));
            
            GetStartGrp(IsGmp, z, n, grpSize, r, lower, 
                        computedRows, lowerMpz, computedRowMpz, stepSize);
        }
        
        pool.push(std::cref(ParallelGlue<RcppParallel::RMatrix<typeElem>, typeElem>), n,
                  std::cref(v), std::ref(parMat), IsGmp, z, r, grpSize, step, nRows, IsSample, 
                  std::cref(mySample), myBigSamp, computedRows, std::ref(computedRowMpz));
        
        pool.join();
        
        FinalTouch(GroupsMat, IsArray, grpSize, r,
                   n, nRows, IsNamed, mySample, myBigSamp, IsGmp);
    } else {
        SerialGlue(n, v, GroupsMat, IsGmp, z, r, grpSize, nRows, IsSample,
                   mySample, myBigSamp, computedRows, computedRowMpz, IsArray, IsNamed);
    }
}

// [[Rcpp::export]]
SEXP ComboGroupsCountCpp(SEXP Rv, SEXP RNumGroups) {
    
    int n, numGroups;
    VecType myType = VecType::Integer;
    CleanConvert::convertPrimitive(RNumGroups, numGroups, "numGroups");
    
    std::vector<double> vNum;
    std::vector<int> vInt;
    
    SetType(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    
    if (n % numGroups != 0) {
        Rcpp::stop("The length of v (if v is a vector) or v (if v is a"
                       " scalar) must be divisible by numGroups");
    }
    
    const int grpSize = n / numGroups;
    const double computedRows = numGroupCombs(n, numGroups, grpSize);
    bool IsGmp = (computedRows > Significand53);
    
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);
    
    if (IsGmp) {
        mpz_set_ui(computedRowMpz, 1);
        numGroupCombsGmp(computedRowMpz, n, numGroups, grpSize);
    }
    
    return GetCount(IsGmp, computedRowMpz, computedRows);
}

// [[Rcpp::export]]
SEXP ComboGroupsRcpp(SEXP Rv, SEXP RNumGroups, SEXP RRetType, SEXP Rlow, SEXP Rhigh,
                     bool IsFactor, SEXP Rparallel, SEXP RNumThreads, int maxThreads,
                     bool IsSample, SEXP RindexVec, SEXP RmySeed, SEXP RNumSamp,
                     Rcpp::Function baseSample, SEXP RNamed) {
    
    int n, numGroups;
    VecType myType = VecType::Integer;
    CleanConvert::convertPrimitive(RNumGroups, numGroups, "numGroups");
    bool Parallel = CleanConvert::convertLogical(Rparallel, "Parallel");
    bool IsNamed = (IsSample) ? CleanConvert::convertLogical(RNamed, "namedSample") : false;
    
    std::vector<double> vNum;
    std::vector<int> vInt;
    
    SetType(myType, Rv);
    SetValues(myType, vInt, vNum, n, Rv);
    
    if (n % numGroups != 0) {
        Rcpp::stop("The length of v (if v is a vector) or v (if v is a"
                       " scalar) must be divisible by numGroups");
    }
    
    const int grpSize = n / numGroups;
    const double computedRows = numGroupCombs(n, numGroups, grpSize);
    bool IsGmp = (computedRows > Significand53);
    
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);
    
    if (IsGmp) {
        mpz_set_ui(computedRowMpz, 1);
        numGroupCombsGmp(computedRowMpz, n, numGroups, grpSize);
    }
    
    double lower = 0, upper = 0;
    bool bLower = false, bUpper = false;
    auto lowerMpz = FromCpp14::make_unique<mpz_t[]>(1);
    auto upperMpz = FromCpp14::make_unique<mpz_t[]>(1);
    mpz_init(lowerMpz[0]); mpz_init(upperMpz[0]);
    
    if (!IsSample) {
        SetBounds(Rlow, Rhigh, IsGmp, bLower, bUpper, lower, upper,
                  lowerMpz.get(), upperMpz.get(), computedRowMpz, computedRows);
    }
    
    std::vector<int> startZ;
    double dblLower = lower;
    if (!IsGmp) mpz_set_d(lowerMpz[0], dblLower);
    
    if (bLower && mpz_cmp_ui(lowerMpz[0], 0) > 0) {
        if (IsGmp)
            startZ = nthComboGroupGmp(n, grpSize, numGroups, lowerMpz[0], computedRowMpz);
        else
            startZ = nthComboGroup(n, grpSize, numGroups, lower, computedRows);
    } else {
        startZ.resize(n);
        std::iota(startZ.begin(), startZ.end(), 0);
    }
    
    int nRows = 0;
    
    if (!IsSample) {
        double userNumRows = 0;
        SetNumResults(IsGmp, bLower, bUpper, false, upperMpz.get(), lowerMpz.get(),
                      lower, upper, computedRows, computedRowMpz, nRows, userNumRows);
    }
    
    const std::string retType = Rcpp::as<std::string>(RRetType);
    bool IsArray = false;

    if (retType != "3Darray" && retType != "matrix") {
        Rcpp::stop("retType must be '3Darray' or 'matrix'");
    } else {
        if (retType == "3Darray")
            IsArray = true;
    }
    
    std::size_t sampSize;
    std::vector<double> mySample;
    
    if (IsSample) {
        // sampleLimit defined in GmpCombPermUtils.cpp.. see comments in GmpCombPermUtils.cpp
        IsGmp = computedRows > sampleLimit;
        SetRandomSample(RindexVec, RNumSamp, sampSize, IsGmp, computedRows, mySample, baseSample);
    }
    
    int nThreads = 1;
    const int limit = (IsSample) ? 2 : 20000;
    const int numResults = (IsSample) ? sampSize : nRows;
    SetThreads(Parallel, maxThreads, numResults, myType, nThreads, RNumThreads, limit);
    
    const std::size_t bigSampSize = (IsSample && IsGmp) ? sampSize : 1;
    auto myVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);
    
    for (std::size_t i = 0; i < bigSampSize; ++i)
        mpz_init(myVec[i]);
        
    if (IsSample) {
        nRows = sampSize;
        SetRandomSampleMpz(RindexVec, RmySeed, sampSize, IsGmp, computedRowMpz, myVec.get());
    }
    
    if (myType > VecType::Logical) {
        const SEXP sexpCopy(Rcpp::clone(Rv));
        RCPP_RETURN_VECTOR(NonThreadSafe, sexpCopy, startZ, n, IsGmp, numGroups, grpSize,
                           nRows, IsSample, mySample, myVec.get(), computedRows, 
                           computedRowMpz, IsArray, IsNamed);
    } else if (myType == VecType::Logical) {
        Rcpp::LogicalMatrix boolGroupsMat(nRows, n);
        GroupsMaster(n, vInt, boolGroupsMat, startZ, numGroups, grpSize, nRows, Parallel,
                     nThreads, IsGmp, lower, lowerMpz[0], computedRows, computedRowMpz,
                     IsSample, mySample, myVec.get(), IsArray, IsNamed);
        return boolGroupsMat;
    } else if (myType == VecType::Integer) {
        Rcpp::IntegerMatrix intGroupsMat = Rcpp::no_init_matrix(nRows, n);
        GroupsMaster(n, vInt, intGroupsMat, startZ, numGroups, grpSize, nRows, Parallel,
                     nThreads, IsGmp, lower, lowerMpz[0], computedRows, computedRowMpz,
                     IsSample, mySample, myVec.get(), IsArray, IsNamed);

        if (IsFactor) {SetFactorClass(intGroupsMat, Rv);}
        return intGroupsMat;
    } else {
        Rcpp::NumericMatrix numGroupsMat = Rcpp::no_init_matrix(nRows, n);
        GroupsMaster(n, vNum, numGroupsMat, startZ, numGroups, grpSize, nRows, Parallel,
                     nThreads, IsGmp, lower, lowerMpz[0], computedRows, computedRowMpz,
                     IsSample, mySample, myVec.get(), IsArray, IsNamed);
        return numGroupsMat;
    }
}
