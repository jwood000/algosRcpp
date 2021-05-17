#ifndef ERATOSTHENES_H
#define ERATOSTHENES_H

#include <cmath>
#include <vector>
#include <cstdint>

namespace PrimeSieve{
    template <typename typePrime>
    void PrimeSieveMaster(std::int_fast64_t minNum, std::int_fast64_t maxNum, std::vector<typePrime> &primes, 
                          std::vector<std::vector<typePrime>> &primeList, bool &Parallel,
                          int nThreads = 1, int maxThreads = 1, int maxCores = 1);
    
    template <typename typePrime>
    void sqrtBigPrimes(int sqrtBound, bool bAddZero, bool bAddExtraPrime,
                       bool bAddTwo, std::vector<typePrime> &sievePrimes);
}

#endif
