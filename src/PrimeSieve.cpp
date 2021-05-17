#include "Eratosthenes.h"
#include "CleanConvert.h"

// [[Rcpp::export]]
SEXP EratosthenesRcpp(SEXP Rb1, SEXP Rb2, SEXP RNumThreads, int maxCores, int maxThreads) {
    
    double bound1, bound2;
    std::int_fast64_t myMax, myMin;
    CleanConvert::convertPrimitive(Rb1, bound1, "bound1", true, false);
    
    if (Rf_isNull(Rb2)) {
        bound2 = 1;
    } else {
        CleanConvert::convertPrimitive(Rb2, bound2, "bound2", true, false);
    }
    
    if (bound1 > bound2) {
        myMax = static_cast<std::int_fast64_t>(std::floor(bound1));
        myMin = static_cast<std::int_fast64_t>(std::ceil(bound2));
    } else {
        myMax = static_cast<std::int_fast64_t>(std::floor(bound2));
        myMin = static_cast<std::int_fast64_t>(std::ceil(bound1));
    }
    
    if (myMax <= 1)
        return Rcpp::IntegerVector();
    
    if (myMin <= 2) myMin = 1;
    
    if (myMin == myMax) {
        if (myMax % 2) {
            ++myMax;
        } else {
            if (myMax > std::numeric_limits<int>::max())
                return Rcpp::NumericVector();
            else
                return Rcpp::IntegerVector();
        }
    }
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads");
    
    std::size_t numPrimes = 0u;
    std::vector<std::size_t> runningCount;
    runningCount.push_back(0u);
    std::size_t numSects = nThreads;
    bool Parallel = false;
    
    if (myMax > std::numeric_limits<int>::max()) {
        std::vector<std::vector<double>> primeList(numSects, std::vector<double>());
        std::vector<double> tempPrimes;
        
        PrimeSieve::PrimeSieveMaster(myMin, myMax, tempPrimes, primeList,
                                     Parallel, nThreads, maxThreads, maxCores);
        
        if (Parallel) {
            for (std::size_t i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }
            
            Rcpp::NumericVector primes(numPrimes);
            Rcpp::NumericVector::iterator priBeg = primes.begin();
            
            for (std::size_t i = 0; i < numSects; ++i)
                std::move(primeList[i].cbegin(), primeList[i].cend(), priBeg + runningCount[i]);
            
            return primes;
        } else {
            Rcpp::NumericVector primes(tempPrimes.cbegin(), tempPrimes.cend());
            return primes;
        }
    } else {
        std::vector<std::vector<int>> primeList(numSects, std::vector<int>());
        std::vector<int> tempPrimes;
        
        PrimeSieve::PrimeSieveMaster(myMin, myMax, tempPrimes, primeList,
                                     Parallel, nThreads, maxThreads, maxCores);
        
        if (Parallel) {
            for (std::size_t i = 0; i < numSects; ++i) {
                numPrimes += primeList[i].size();
                runningCount.push_back(numPrimes);
            }
            
            Rcpp::IntegerVector primes(numPrimes);
            Rcpp::IntegerVector::iterator priBeg = primes.begin();
            
            for (std::size_t i = 0; i < numSects; ++i)
                std::move(primeList[i].cbegin(), primeList[i].cend(), priBeg + runningCount[i]);
            
            return primes;
        } else {
            Rcpp::IntegerVector primes(tempPrimes.cbegin(), tempPrimes.cend());
            return primes;
        }
    }
}
