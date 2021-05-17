#include "Eratosthenes.h"
#include "PrimeFactorizeSieve.h"
#include "EulerPhiSieve.h"
#include <RcppThread.h>
#include "CleanConvert.h"
#include <RcppThread.h>

template <typename typeInt, typename typeReturn, typename typeRcpp>
void MotleyMaster(typeInt myMin, typeReturn myMax, bool isEuler,
                  typeRcpp &EulerPhis, std::vector<typeInt> &numSeq,
                  std::vector<std::vector<typeInt>> &primeList,
                  int nThreads, int maxThreads) {
    
    bool Parallel = false;
    std::int_fast64_t myRange = (myMax - myMin) + 1;
    typeInt offsetStrt = 0;
    
    if (nThreads > 1 && maxThreads > 1 && myRange >= 20000) {
        Parallel = true;
        
        if (nThreads > maxThreads)
            nThreads = maxThreads;
        
        // Ensure that each thread has at least 10000
        if ((myRange / nThreads) < 10000)
            nThreads = myRange / 10000;
    }
    
    int sqrtBound = std::sqrt(static_cast<double>(myMax));
    std::vector<typeInt> primes;
    PrimeSieve::sqrtBigPrimes(sqrtBound, false, true, true, primes);
    
    if (Parallel) {
        RcppThread::ThreadPool pool(nThreads);
        typeInt lower = myMin;
        typeInt chunkSize = myRange / nThreads;
        typeReturn upper = lower + chunkSize - 1;
        
        for (int ind = 0; ind < (nThreads - 1); offsetStrt += chunkSize, 
             lower = (upper + 1), upper += chunkSize, ++ind) {
            if (isEuler) {
                pool.push(std::cref(MotleyPrimes::EulerPhiSieve<typeInt, typeReturn, typeRcpp>), lower,
                          upper, offsetStrt, std::ref(primes), std::ref(numSeq), std::ref(EulerPhis));
            } else {
                pool.push(std::cref(MotleyPrimes::PrimeFactorizationSieve<typeInt>), 
                          lower, static_cast<typeInt>(upper), offsetStrt,
                          std::ref(primes), std::ref(primeList));
            }
        }
        
        if (isEuler) {
            pool.push(std::cref(MotleyPrimes::EulerPhiSieve<typeInt, typeReturn, typeRcpp>), lower,
                      myMax, offsetStrt, std::ref(primes), std::ref(numSeq), std::ref(EulerPhis));
        } else {
            pool.push(std::cref(MotleyPrimes::PrimeFactorizationSieve<typeInt>), 
                      lower, static_cast<typeInt>(myMax), offsetStrt,
                      std::ref(primes), std::ref(primeList));
        }
        
        pool.join();
        
    } else {
        if (isEuler) {
            MotleyPrimes::EulerPhiSieve(myMin, myMax, offsetStrt, primes, numSeq, EulerPhis);
        } else {
            MotleyPrimes::PrimeFactorizationSieve(myMin, static_cast<typeInt>(myMax),
                                                  offsetStrt, primes, primeList);
        }
    }
}

template <typename typeInt, typename typeReturn, typename typeRcpp>
SEXP GlueMotley(typeInt myMin, typeReturn myMax, bool isEuler,
                typeRcpp temp, bool keepNames, int nThreads, int maxThreads) {
    
    std::size_t myRange = (myMax - myMin) + 1;
    std::vector<typeReturn> myNames;
    
    if (keepNames) {
        myNames.resize(myRange);
        typeReturn retM = myMin;
        for (std::size_t k = 0; retM <= myMax; ++retM, ++k)
            myNames[k] = retM;
    }
    
    if (isEuler) {
        std::vector<std::vector<typeInt>> tempList;
        typeRcpp EulerPhis(myRange);
        std::vector<typeInt> numSeq(myRange);
        MotleyMaster(myMin, myMax, isEuler, EulerPhis,
                     numSeq, tempList, nThreads, maxThreads);
        
        if (keepNames)
            EulerPhis.attr("names") = myNames;
        
        return EulerPhis;
    } else {
        std::vector<std::vector<typeInt>> 
            primeList(myRange, std::vector<typeInt>());
        typeRcpp tempRcpp;
        std::vector<typeInt> tempVec;
        MotleyMaster(myMin, myMax, isEuler, tempRcpp,
                     tempVec, primeList, nThreads, maxThreads);
        
        Rcpp::List myList = Rcpp::wrap(primeList);
        
        if (keepNames)
            myList.attr("names") = myNames;
        
        return myList;
    }
}

// [[Rcpp::export]]
SEXP MotleyContainer(SEXP Rb1, SEXP Rb2, bool isEuler, 
                     SEXP RNamed, SEXP RNumThreads, int maxThreads) {
    
    double bound1, bound2, myMin, myMax;
    const std::string namedObject = (isEuler) ? "namedVector" : "namedList";
    bool isNamed = CleanConvert::convertLogical(RNamed, namedObject);
    CleanConvert::convertPrimitive(Rb1, bound1, "bound1");
    
    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, "bound2");
    }
    
    if (bound1 > bound2) {
        myMax = std::floor(bound1);
        myMin = std::ceil(bound2);
    } else {
        myMax = std::floor(bound2);
        myMin = std::ceil(bound1);
    }
    
    if (myMax < 2) {
        if (isEuler) {
            Rcpp::IntegerVector z(1, 1);
            if (isNamed)
                z.attr("names") = 1;
            return z;
        } else {
            std::vector<std::vector<int>> trivialRet(1, std::vector<int>());
            Rcpp::List z = Rcpp::wrap(trivialRet);
            if (isNamed)
                z.attr("names") = 1;
            return z;
        }
    }
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads");
    
    if (myMax > std::numeric_limits<int>::max()) {
        std::int64_t intMin = static_cast<std::int64_t>(myMin);
        Rcpp::NumericVector temp;
        return GlueMotley(intMin, myMax, isEuler, temp, isNamed, nThreads, maxThreads);
    } else {
        int intMin = static_cast<int>(myMin);
        int intMax = static_cast<int>(myMax);
        Rcpp::IntegerVector temp;
        return GlueMotley(intMin, intMax, isEuler, temp, isNamed, nThreads, maxThreads);
    }
}
