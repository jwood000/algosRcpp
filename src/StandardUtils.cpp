#include "CleanConvert.h"

void SetType(VecType &myType, const SEXP &Rv) {
    
    switch(TYPEOF(Rv)) {
        case LGLSXP: {
            myType = VecType::Logical;
            break;
        }
        case INTSXP: {
            myType = VecType::Integer;
            break;
        }
        case REALSXP: {
            myType = VecType::Numeric;
            break;
        }
        case STRSXP: {
            myType = VecType::Character;
            break;
        }
        case CPLXSXP: {
            myType = VecType::Complex;
            break;
        }
        case RAWSXP: {
            // Vectors of class bigZ and mpfr cause a lot of headaches, and for this
            // we simply exclude all raw vectors that have any attributes. If you
            // think there is a clean solution for including these cases, please
            // contact me @ jwood000@gmail.com. N.B., see commit 655 which includes
            // a function for returning a matrix of class bigz. I observed terrible
            // performance compared to simply converting to a character vector.
            Rcpp::RawVector tempRaw(Rcpp::clone(Rv));
            
            if (tempRaw.attributeNames().size() == 0) {
                myType = VecType::Raw;
                break;
            }
        }
        default: {
            Rcpp::stop("Only atomic types are supported for v");
        }
    }
}

std::vector<int> zUpdateIndex(SEXP x, SEXP y, int m) {
    
    Rcpp::IntegerVector nonZeroBased(m);
    constexpr double myTolerance = 8 * std::numeric_limits<double>::epsilon();
    const int n1 = Rf_length(x) - 1;
    
    switch (TYPEOF(x)) {
    case LGLSXP: {
        Rcpp::LogicalVector xBool(Rcpp::clone(x));
        Rcpp::LogicalVector yBool(Rcpp::clone(y));
        
        for (int j = 0; j < m; ++j) {
            int ind = 0;
            
            while (ind < n1 && xBool[ind] != yBool[j])
                ++ind;
            
            nonZeroBased[j] = ind + 1;
        }
        
        break;
    }
    case INTSXP: {
        Rcpp::IntegerVector xInt(Rcpp::clone(x));
        Rcpp::IntegerVector yInt(Rcpp::clone(y));
        nonZeroBased = Rcpp::match(yInt, xInt);
        break;
    }
    case REALSXP: {
        Rcpp::NumericVector xNum(Rcpp::clone(x));
        Rcpp::NumericVector yNum(Rcpp::clone(y));
        
        for (int j = 0; j < m; ++j) {
            int ind = 0;
            
            while (ind < n1 && std::abs(xNum[ind] - yNum[j]) > myTolerance)
                ++ind;
            
            nonZeroBased[j] = ind + 1;
        }
        
        break;
    }
    case STRSXP: {
        Rcpp::CharacterVector xBool(Rcpp::clone(x));
        Rcpp::CharacterVector yBool(Rcpp::clone(y));
        nonZeroBased = Rcpp::match(yBool, xBool);
        break;
    }
    case CPLXSXP: {
        Rcpp::ComplexVector xCmplx(Rcpp::clone(x));
        Rcpp::ComplexVector yCmplx(Rcpp::clone(y));
        
        for (int j = 0; j < m; ++j) {
            int ind = 0;
            bool bTestImg = std::abs(xCmplx[ind].i - yCmplx[j].i) > myTolerance;
            bool bTestReal = std::abs(xCmplx[ind].r - yCmplx[j].r) > myTolerance;
            
            while (ind < n1 && (bTestImg || bTestReal)) {
                ++ind;
                bTestImg = std::abs(xCmplx[ind].i - yCmplx[j].i) > myTolerance;
                bTestReal = std::abs(xCmplx[ind].r - yCmplx[j].r) > myTolerance;
            }
            
            nonZeroBased[j] = ind + 1;
        }
        
        break;
    }
    case RAWSXP: {
        Rcpp::RawVector xRaw(Rcpp::clone(x));
        Rcpp::RawVector yRaw(Rcpp::clone(y));
        
        for (int j = 0; j < m; ++j) {
            int ind = 0;
            
            while (ind < n1 && xRaw[ind] != yRaw[j])
                ++ind;
            
            nonZeroBased[j] = ind + 1;
        }
        
        break;
    }
    default:{
        Rcpp::stop("Only atomic types are supported for v");
    }
    }
    
    std::vector<int> res(m, 0);
    
    for (int i = 0; i < m; ++i)
        res[i] = nonZeroBased[i] - 1;
    
    return res;
}

void SetFactorClass(Rcpp::IntegerMatrix &matInt, const SEXP &Rv) {
    Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
    Rcpp::CharacterVector myClass = testFactor.attr("class");
    Rcpp::CharacterVector myLevels = testFactor.attr("levels");
    matInt.attr("class") = myClass;
    matInt.attr("levels") = myLevels;
}

void SetValues(VecType &myType, std::vector<int> &vInt,
               std::vector<double> &vNum, int &n, const SEXP &Rv) {
    
    if (myType > VecType::Logical) {
        n = Rf_length(Rv);
    } else if (myType == VecType::Logical) {
        vInt = Rcpp::as<std::vector<int>>(Rv);
        n = vInt.size();
    } else {
        if (Rf_length(Rv) == 1) {
            int seqEnd, m1, m2;         // numOnly = true, checkWhole = true, negPoss = true
            CleanConvert::convertPrimitive(Rv, seqEnd, "If v is not a character and of length 1, it", true, true, true);
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            Rcpp::IntegerVector vTemp = Rcpp::seq(m1, m2);
            vNum = Rcpp::as<std::vector<double>>(vTemp);
            myType = VecType::Integer;
        } else {
            vNum = Rcpp::as<std::vector<double>>(Rv);
        }
        
        n = vNum.size();
    }
    
    if (myType == VecType::Integer) {
        bool IsInteger = true;
        
        for (int i = 0; i < n && IsInteger; ++i) {
            if (Rcpp::NumericVector::is_na(vNum[i])) {
                IsInteger = false;
                myType = VecType::Numeric;
            }
        }
        
        if (IsInteger)
            vInt.assign(vNum.cbegin(), vNum.cend());
    }
}

void SetFreqsAndM(SEXP RFreqs, bool &IsMultiset, std::vector<int> &Reps, bool &IsRepetition,
                  int &lenFreqs, std::vector<int> &freqsExpanded, const SEXP &Rm, int n, int &m) {
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        Reps.push_back(1);
    } else {
        IsRepetition = false;
        CleanConvert::convertVector(RFreqs, Reps, "freqs");
        lenFreqs = static_cast<int>(Reps.size());
        bool allOne = std::all_of(Reps.cbegin(), Reps.cend(), 
                                  [](int v_i) {return v_i == 1;});
        if (allOne) {
            IsMultiset = false;
            freqsExpanded = Reps;
        } else {
            IsMultiset = true;
            
            for (int i = 0; i < lenFreqs; ++i)
                for (int j = 0; j < Reps[i]; ++j)
                    freqsExpanded.push_back(i);
        }
    }
    
    if (Rf_isNull(Rm)) {
        if (freqsExpanded.empty()) {
            m = n;
        } else {
            m = freqsExpanded.size();
        }
    } else {
        if (Rf_length(Rm) > 1)
            Rcpp::stop("length of m must be 1");
        
        CleanConvert::convertPrimitive(Rm, m, "m");
    }
}

SEXP CopyRv(const SEXP &Rv, const std::vector<int> &vInt,
            const std::vector<double> &vNum, VecType myType, bool IsFactor) {
    
    if (myType > VecType::Numeric || IsFactor) {
        return Rcpp::clone(Rv);
    } else if (myType == VecType::Integer) {
        return Rcpp::wrap(vInt);
    } else {
        return Rcpp::wrap(vNum);
    }
}

void SetThreads(bool &Parallel, int maxThreads, int nRows,
                VecType myType, int &nThreads, SEXP RNumThreads, int limit) {
    
    const int halfLimit = limit / 2;
    
    // Determined empirically. Setting up threads can be expensive,
    // so we set the cutoff below to ensure threads aren't spawned
    // unnecessarily. We also protect users with fewer than 2 threads
    if ((nRows < limit) || (maxThreads < 2) || myType > VecType::Logical) {
        Parallel = false;
    } else if (!Rf_isNull(RNumThreads)) {
        int userThreads = 1;
        
        if (!Rf_isNull(RNumThreads))
            CleanConvert::convertPrimitive(RNumThreads, userThreads, "nThreads");
        
        if (userThreads > maxThreads)
            userThreads = maxThreads;
        
        // Ensure that each thread has at least halfLimit
        if ((nRows / userThreads) < halfLimit)
            userThreads = nRows / halfLimit;
        
        if (userThreads > 1) {
            Parallel = true;
            nThreads = userThreads;
        } else {
            Parallel = false;
        }
    } else if (Parallel) {
        // We have already ruled out cases when the user has fewer than 2 
        // threads. So if user has exactly 2 threads, we enable them both.
        nThreads = (maxThreads > 2) ? (maxThreads - 1) : maxThreads;
        
        // Ensure that each thread has at least halfLimit
        if ((nRows / nThreads) < halfLimit)
            nThreads = nRows / halfLimit;
    }
}

void SetRandomSample(SEXP RindexVec, SEXP RNumSamp, std::size_t &sampSize,
                     bool IsGmp, double computedRows, std::vector<double> &mySample,
                     Rcpp::Function baseSample) {
    
    // We must treat gmp case special. We first have to get the size of our sample
    // vector, as we have to declare a mpz_t array with known size. You will note
    // that in the base case below, we simply populate mySample, otherwise we just
    // get the size. This size var will be used in the next block (If (IsGmp)...)
    if (Rf_isNull(RindexVec)) {
        if (Rf_isNull(RNumSamp))
            Rcpp::stop("n and sampleVec cannot both be NULL");
        
        if (Rf_length(RNumSamp) > 1)
            Rcpp::stop("length of n must be 1. For specific combinations, use sampleVec.");
        
        int nPass;
        CleanConvert::convertPrimitive(RNumSamp, nPass, "n");
        sampSize = static_cast<std::size_t>(nPass);
        
        if (!IsGmp) {
            if (nPass > computedRows)
                Rcpp::stop("n exceeds the maximum number of possible results");
            
            Rcpp::NumericVector tempSamp = baseSample(computedRows, nPass);
            mySample = Rcpp::as<std::vector<double>>(tempSamp);
        }
    } else {
        if (IsGmp) {
            switch (TYPEOF(RindexVec)) {
            case RAWSXP: {
                const char* raw = (char*) RAW(RindexVec);
                sampSize = ((int*) raw)[0];
                break;
            }
            default: {
                sampSize = LENGTH(RindexVec);
            }
            }
        } else {                                             // numOnly = false
            CleanConvert::convertVector(RindexVec, mySample, "sampleVec", false);
            sampSize = mySample.size();
            
            double myMax = *std::max_element(mySample.cbegin(), mySample.cend());
            
            if (myMax > computedRows) {
                Rcpp::stop("One or more of the requested values in sampleVec "
                               "exceeds the maximum number of possible results");
            }
        }
        
        if (sampSize > std::numeric_limits<int>::max())
            Rcpp::stop("The number of rows cannot exceed 2^31 - 1");
    }
    
    // Get zero base index
    for (auto &s: mySample)
        --s;
}
