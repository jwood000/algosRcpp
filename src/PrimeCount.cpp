#include "CleanConvert.h"
#include "PrimesSegSieve.h"
#include "PhiTinyLookup.h"
#include "Eratosthenes.h"
#include <RcppThread.h>
#include <mutex>
#include <cmath>

// "MasterPrimeCount" based off of the highly optimized
// "pi_legendre.cpp" algorithm by Kim Walisch, which calculates
// the numbers of primes less than n using Legendre's formula.
// Kim Walisch's official github repo for "pi_legendre" is:
//                      https://github.com/kimwalisch/primecount

namespace PrimeCounting {
    
    // PiPrime is very similar to the PrimeSieveSmall only we are not
    // considering a range. That is, we are only concerned with finding
    // the number of primes less than maxNum. We are also only counting
    // primes instead of generating them.
    std::int64_t PiPrime (std::int64_t maxNum) {
        
        constexpr std::int_fast64_t segSize = Almost210L1Cache;
        constexpr std::size_t nWheels = N_WHEELS210_PER_SEG;
        constexpr std::size_t szWheel210 = SZ_WHEEL210;
        const int sqrtBound = static_cast<int>(std::sqrt(static_cast<double>(maxNum)));
        
        // the wheel already has the first 4 primes marked as
        // false, so we need to account for them here. N.B.
        // the calling functions has checks for cases where
        // maxNum is less than 11, so no need to check here.
        std::int64_t count = 4;
        
        std::vector<std::int64_t> smallPrimes, nextStrt;
        const std::int64_t flrMaxNum = segSize * std::floor(maxNum / segSize);
        
        std::size_t ind = 1;
        for (; smallPrimeBase[ind] <= sqrtBound; ++ind)
            smallPrimes.push_back(smallPrimeBase[ind]);
        
        smallPrimes.push_back(smallPrimeBase[ind]);
        std::vector<char> sieve(segSize, 1);
        sieve[1] = 0;
        
        std::int64_t sqrPrime = 9;
        std::int64_t lowerBnd = 0, upperBnd = segSize, myNum = 1;
        std::size_t p = 1;
        
        for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
            upperBnd = lowerBnd + segSize;
            
            for (; sqrPrime <= upperBnd; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = smallPrimes[p] * smallPrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                std::int64_t j = nextStrt[i];
                for (std::int64_t k = smallPrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            for (std::size_t q = 0; q < nWheels; ++q)
                for (std::size_t w = 0; w < szWheel210; myNum += ARR_WHEEL210[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        ++count;
                    
            std::fill(sieve.begin(), sieve.end(), 1);
        }
        
        if (lowerBnd < maxNum) {
            for (; sqrPrime <= maxNum; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = smallPrimes[p] * smallPrimes[p];
            }
            
            for (std::size_t i = 3; i < nextStrt.size(); ++i) {
                std::int64_t j = nextStrt[i];
                for (std::int64_t k = smallPrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = 0;
                
                nextStrt[i] = j - segSize;
            }
            
            for (std::size_t q = 0; q < nWheels && myNum <= maxNum; ++q)
                for (std::size_t w = 0; w < szWheel210 && myNum <= maxNum; myNum += ARR_WHEEL210[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        ++count;
        }
        
        return count;
    }

    constexpr int MAX_A = 100;
    std::array<std::vector<uint16_t>, MAX_A> phiCache;
    std::mutex theBlocker;
    
    // Increment MAX_A, so we can have easier access to indexing. E.g.
    // when a = 3 (i.e. the number of primes is 3), instead of accessing
    // the third entry in phiTiny like phiTiny[a - 1], we simply use "a"
    // as our index (i.e. phiTiny[a]). The first entry in phiTiny will
    // be empty. The same goes for phiPrimes and phiPi.
    std::vector<std::int64_t> phiPrimes;
    std::vector<std::int64_t> phiPi;
    
    // The arrays below are utilized by phiTiny
    // primeProds[n] = \prod_{i=1}^{n} primes[i]
    // myTotients[n] = \prod_{i=1}^{n} (primes[i] - 1)
    std::int64_t phiTinyCalc(std::int64_t x, std::int64_t a) {
        static constexpr std::array<int, 7> myTotients = {{1, 1, 2, 8, 48, 480, 5760}};
        static constexpr std::array<int, 7> primeProds = {{1, 2, 6, 30, 210, 2310, 30030}};
        
        std::int64_t pp = primeProds[a];
        return (x / pp) * myTotients[a] + phiTiny[a][x % pp];
    }
    
    void updateCache(std::uint64_t x, std::uint64_t a, std::int64_t mySum) {
        if (a < phiCache.size() &&
            x <= std::numeric_limits<uint16_t>::max()) {
            // Protect phiCache while its being updated
            std::lock_guard<std::mutex> guard(theBlocker);
            if (x >= phiCache[a].size()) phiCache[a].resize(x + 1, 0);
            phiCache[a][x] = static_cast<uint16_t>(std::abs(mySum));
        }
    }
    
    bool isCached(std::uint64_t x, std::uint64_t a) {
        return a < phiCache.size() &&
               x < phiCache[a].size() &&
               phiCache[a][x];
    }
    
    bool isPix(std::int64_t x, std::int64_t a) {
        return x < static_cast<std::int64_t>(phiPi.size()) &&
               x < (phiPrimes[a + 1] * phiPrimes[a + 1]);
    }
    
    std::int64_t getStrt(std::int64_t y) {
        static constexpr std::array<int, 7> myTinyPrimes = {{0, 2, 3, 5, 7, 11, 13}};
        static constexpr std::array<int, 13> myTinyPi = {{0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5}};
        
        if (y >= myTinyPrimes.back())
            return phiTinySize;
        else
            return myTinyPi[y];
    }
    
    template <int SIGN>
    std::int64_t phiWorker(std::int64_t x, std::int64_t a) {
        if (x <= phiPrimes[a]) {
            return SIGN;
        } else if (a <= phiTinySize) {
            return phiTinyCalc(x, a) * SIGN;
        } else if (isPix(x, a)) {
            return (phiPi[x] - a + 1) * SIGN;
        } else if (isCached(x, a)) {
            return phiCache[a][x] * SIGN;
        } else {
            std::int64_t sqrtx = static_cast<std::int64_t>(std::sqrt(static_cast<double>(x)));
            std::int64_t piSqrtx = a;
            std::int64_t strt = getStrt(sqrtx);
            
            if (sqrtx < static_cast<std::int64_t>(phiPi.size()))
                piSqrtx = std::min(static_cast<std::int64_t>(phiPi[sqrtx]), a);
            
            std::int64_t mySum = 0;
            mySum += (piSqrtx - a) * SIGN;
            mySum += phiTinyCalc(x, strt) * SIGN;
            
            for (std::int64_t i = strt; i < piSqrtx; ++i) {
                std::int64_t x2 = x / phiPrimes[i + 1];
                
                if (isPix(x2, i)) {
                    mySum += (phiPi[x2] - i + 1) * -SIGN;
                } else {
                    mySum += phiWorker<-SIGN>(x2, i);
                }
            }
            
            updateCache(x, a, mySum);
            return mySum;
        }
    }
    
    std::int64_t phiForeman(std::int64_t lowerBound, std::int64_t upperBound, std::int64_t x) {
        
        std::int64_t mySum = 0;
        
        for (std::int64_t i = lowerBound; i < upperBound; ++i)
            mySum += phiWorker<-1>(x / phiPrimes[i + 1], i);
        
        return mySum;
    }
    
    const double getChunkFactor(std::int64_t x) {
        static constexpr std::array<double, 9> nums = {{1e10, 1e12, 2e13, 5e13, 8e13, 1e14, 5e14, 1e15, 1e16}};
        static constexpr std::array<double, 9> factor = {{1.3, 1.2, 1.1, 1.07, 1.05, 1.01, 1.007, 1.006, 1.005}};
        
        const auto it = std::upper_bound(nums.cbegin(), nums.cend(), static_cast<double>(x));
        return std::log(factor[it - nums.cbegin()]);
    }
    
    std::int64_t phiMaster(std::int64_t x, std::int64_t a, int nThreads, bool Parallel) {
        
        const std::int64_t sqrtx = static_cast<std::int64_t>(std::sqrt(static_cast<double>(x)));
        const std::int64_t piSqrtx = std::min(static_cast<std::int64_t>(phiPi[sqrtx]), a);
        const std::int64_t strt = getStrt(sqrtx);
        std::int64_t mySum = phiTinyCalc(x, strt) + piSqrtx - a;
        
        if (Parallel) {
            constexpr double divLim = 10000000000.0;    //  1e10
            
            // We know x, t, and we need to solve for n.
            // x / (1.5^t)^n < 1e10 -->>> 1e10 * (1.5^t)^n = x --->>> x / 1e10 = (1.5^t)^n
            // log(x / 1e10) = n * log(1.5^t) --->> log(x / 1e10) / log(1.5^t) :
            //                 log(x / 1e10) / (t * log(1.5))
            std::int64_t myRange = (piSqrtx - strt) + 1;
            std::int64_t lower = strt;
            std::int64_t upper;
            RcppThread::ThreadPool pool(nThreads);
            std::vector<std::future<std::int64_t>> myFutures;
            
            if (x > divLim) {
                const double chunkFactor = getChunkFactor(x);
                const int nLoops = 1 + std::ceil(std::log(x / divLim) / (nThreads * chunkFactor));
                const double multOne = std::exp(std::log(myRange) / (nThreads * nLoops));
                const int firstStep = std::ceil(std::log10(Significand53 / x)) + 1;
                std::int64_t base = lower;
                
                //*********************** Begin of firstThreads *****************************
                // Typically, the size of multOne will be very small resulting in not very
                // meaningful chunk sizes in the first few loops. Because of this, we 
                // gurantee that each of the first threads will have at least firstStep 
                // iterations and that the threads following will have at least firstStep + 1
                
                // Here we have m = multOne, p = power, and s = firstStep:
                // m^(p +  1) - m^p = m^p * (m - 1) > (s + 1)  -->>  m^p > (s + 1) / (m - 1)  -->>
                // p * log(m) > log((s + 1) / (m - 1))  -->>  p > log((s + 1) / (m - 1)) / log(m)
                double power = std::ceil(std::log((firstStep + 1) / (multOne - 1)) / std::log(multOne)) + 1;
                
                // Same as above: lower + (upper - lower) / s   and given:   m^p + lower = upper
                //   -->>  (lower + m^p / s)  -->>  firstLoops = ceil(((m^p + lower) / s) / nThreads)
                const int firstLoops = std::ceil(((std::pow(multOne, power) + 
                                                 base) / firstStep) / nThreads);
                upper = lower + firstStep;
                
                for (int i = 0; i < firstLoops; ++i) {
                    for (int j = 0; j < nThreads; lower = upper, upper += firstStep, ++j)
                        myFutures.push_back(pool.pushReturn(phiForeman, lower, upper, x));
                    
                    pool.wait();
                }
                
                //*********************** End of firstThreads *****************************
                //*************************************************************************
                //*********************** Begin of midThreads *****************************
                
                while (std::pow(multOne, power) < (upper - base))
                    ++power;
                
                double dblChunk = std::pow(multOne, power);
                upper = static_cast<std::int64_t>(dblChunk) + base;

                while ((dblChunk * std::pow(multOne, nThreads - 1) + base) < piSqrtx) {
                    for (int j = 0; j < nThreads; lower = upper, ++j, 
                            dblChunk *= multOne, upper = static_cast<std::int64_t>(dblChunk) + base) {
                        myFutures.push_back(pool.pushReturn(phiForeman, lower, upper, x));
                    }

                    pool.wait();
                }

                for (int j = 0; j < (nThreads - 1) && upper < piSqrtx; lower = upper, ++j, 
                        dblChunk *= multOne, upper = static_cast<std::int64_t>(dblChunk) + base) {
                    myFutures.push_back(pool.pushReturn(phiForeman, lower, upper, x));
                }

                myFutures.push_back(pool.pushReturn(phiForeman, lower, piSqrtx, x));

                for (std::size_t j = 0; j < myFutures.size(); ++j)
                    mySum += myFutures[j].get();

                pool.join();
                
            } else {
                std::int64_t chunk = myRange / nThreads;
                upper = lower + chunk - 1;
                
                for (int j = 0; j < (nThreads - 1); lower = upper, upper += chunk, ++j)
                    myFutures.push_back(pool.pushReturn(phiForeman, lower, upper, x));
                
                myFutures.push_back(pool.pushReturn(phiForeman, lower, piSqrtx, x));
                
                for (std::size_t j = 0; j < myFutures.size(); ++j)
                    mySum += myFutures[j].get();
                
                pool.join();
            }
        } else {
            mySum += phiForeman(strt, piSqrtx, x);
        }
        
        return mySum;
    }
    
    // All values verified by Kim Walisch's primecount library (nThreads = 8)
    //  10^9 -->>          50,847,534   -->>   811.6 microseconds
    // 10^10 -->>         455,052,511   -->>   5.362 milliseconds
    // 10^11 -->>       4,118,054,813   -->>   21.45 milliseconds
    // 10^12 -->>      37,607,912,018   -->>   119.4 milliseconds
    // 10^13 -->>     346,065,536,839   -->>   868.4 milliseconds
    // 10^14 -->>   3,204,941,750,802   -->>   6.724 seconds
    // 10^15 -->>  29,844,570,422,669   -->>  49.554 seconds
    // MAX VALUE (2^53 - 1) -->> 
    //            252,252,704,148,404   -->> 352.862 seconds
    
    std::int64_t MasterPrimeCount(std::int64_t n, int nThreads = 1, int maxThreads = 1) {
        
        const std::int64_t sqrtBound = static_cast<std::int64_t>(std::sqrt(static_cast<double>(n)));
        std::vector<std::int64_t> resetPhiPrimes;
        PrimeSieve::sqrtBigPrimes(sqrtBound, true, false, true, resetPhiPrimes);
        phiPrimes = resetPhiPrimes;
        
        phiPi.resize(sqrtBound + 1);
        std::int64_t count = 0;
        const std::int64_t maxPrime = phiPrimes.back();
        
        for (std::int64_t i = 1; i <= maxPrime; ++i) {
            if (i >= phiPrimes[count + 1])
                ++count;
            
            phiPi[i] = count;
        }
        
        for (std::int64_t i = (maxPrime + 1); i <= sqrtBound; ++i)
            phiPi[i] = count;
        
        bool Parallel = false;
        
        if (nThreads > 1 && maxThreads > 1) {
            Parallel = true;
            if (nThreads > maxThreads) {nThreads = maxThreads;}
            if ((maxThreads < 2) || (n < 1e7)) {Parallel = false;}
        }
        
        const std::int64_t piSqrt = PiPrime(sqrtBound);
        const std::int64_t phiSqrt = phiMaster(n, piSqrt, nThreads, Parallel);
        const std::int64_t int64result = piSqrt + phiSqrt - 1;
        
        return int64result;
    }
}

//[[Rcpp::export]]
SEXP PrimeCountRcpp(SEXP Rn, SEXP RNumThreads, int maxThreads) {
    double dblNum;
    CleanConvert::convertPrimitive(Rn, dblNum, "n");
    const std::int64_t n = static_cast<std::int64_t>(dblNum);
    
    if (n < 100000) {
        if (n < 10) {
            if (n == 1)
                return Rcpp::wrap(0);
            else if (n == 2)
                return Rcpp::wrap(1);
            else if (n < 5)
                return Rcpp::wrap(2);
            else if (n < 7)
                return Rcpp::wrap(3);
            else
                return Rcpp::wrap(4);
        }
        
        return Rcpp::wrap(static_cast<int>(PrimeCounting::PiPrime(n)));
    }
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads");
    
    std::int64_t result = PrimeCounting::MasterPrimeCount(n, nThreads, maxThreads);
    
    if (result > std::numeric_limits<int>::max()) {
        return Rcpp::wrap(static_cast<double>(result));
    } else {
        return Rcpp::wrap(static_cast<int>(result));
    }
}
