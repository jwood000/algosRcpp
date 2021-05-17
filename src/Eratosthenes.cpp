#include "PrimesSegSieve.h"
#include "Wheel.h"
#include <cmath>
#include <deque>
#include <RcppThread.h>

// "PrimeSieve" implements a simple segmented version of the Sieve of 
// Eratosthenes (original implementation authored by Kim Walisch). An
// overview of this method can be found here: 
//                      http://primesieve.org/segmented_sieve.html
// Kim Walisch's official github repo for the primesieve is:
//                      https://github.com/kimwalisch/primesieve

namespace PrimeSieve {

    constexpr int L1_CACHE_SIZE = 32768;
    
    //************************ Prime Counting Estimates ************************
    // These numbers were obtained empirically using the prime number theorem
    // along with a prime counting function. Values were computed in each
    // range below in cutPoints. For example in the range [40000, 120000),
    // if we increment by 10 and calculate the following: 
    //                       pi(x) - x / log(x) = e
    // We then calculate the percentage of x / log(x) or : e / (x / log(x))
    // for every value and find the max:
    //             a <- sapply(seq(40000, 120000, 10), primeCount)
    //             b <- sapply(seq(40000, 120000, 10), function(x) x/log(x))
    //             myDiff <- a - b;  max(myDiff / b) ## [1] 0.1153694
    
    constexpr std::array<double, 15> PERCINC = {{0.2500, 0.1160, 0.1030, 0.0850, 0.0712,
                                                 0.0614, 0.0538, 0.0495, 0.0480, 0.0431, 
                                                 0.0392, 0.0360, 0.0332, 0.0309, 0.0292}};
    
    constexpr std::array<double, 15> CUTPOINTS = {{40000.0,            120000.0,
                                                   1000000.0,          10000000.0,
                                                   100000000.0,        1000000000.0,
                                                   5000000000.0,       10000000000.0,
                                                   100000000000.0,     1000000000000.0,
                                                   10000000000000.0,   100000000000000.0,
                                                   1000000000000000.0, 5000000000000000.0,
                                                   10000000000000000.0}};
    
    // The following function is based off of the prime number theorem
    inline std::size_t EstimatePiPrime(double minNum, double maxNum) {
        const auto it = std::upper_bound(CUTPOINTS.cbegin(), CUTPOINTS.cend(), maxNum);
        const std::size_t myIndex = it - CUTPOINTS.cbegin();
        double dblRes = std::ceil((maxNum / std::log(maxNum)) * (1 + PERCINC[myIndex]));
        
        if (minNum > 1000)
            dblRes -= std::floor((minNum / std::log(minNum)) * (1 + PERCINC[myIndex]));
        
        std::size_t result = dblRes;
        return result;
    }
    
    template <typename typePrime>
    void PrimeSieveSmall(std::int_fast64_t minNum, std::int_fast64_t maxNum,
                         const std::vector<std::int_fast64_t> &sievePrimes,
                         std::vector<typePrime> &myPrimes) {
        
        const std::size_t szWheel2310 = SZ_WHEEL2310;
        const std::size_t nWheels = N_WHEELS2310_PER_SEG;
        const std::size_t myReserve = EstimatePiPrime(static_cast<double>(minNum),
                                                      static_cast<double>(maxNum));
        myPrimes.reserve(myReserve);
        
        if (maxNum <= smallPrimeBase.back()) {
            std::size_t ind = 0u;
            for (; smallPrimeBase[ind] < minNum; ++ind) {}
            
            for (; smallPrimeBase[ind] <= maxNum; ++ind) 
                myPrimes.push_back(static_cast<typePrime>(smallPrimeBase[ind]));
        } else {
            
            const std::int_fast64_t segSize = Almost2310L1Cache;
            const std::int_fast64_t flrMaxNum = static_cast<std::int_fast64_t>(segSize * 
                                                                     std::floor(static_cast<double>(maxNum) / segSize));
            
            std::int_fast64_t lowerBnd = static_cast<std::int_fast64_t>(segSize * 
                std::floor(static_cast<double>(minNum) / segSize));
            
            std::vector<char> sieve(segSize, 1);
            
            if (minNum < wheel2310PrimeLim) {
                sieve[1] = 0;
                std::size_t ind = 0u;
                for (; smallPrimeBase[ind] < minNum; ++ind) {}
                
                for (; smallPrimeBase[ind] < wheel2310PrimeLim; ++ind)
                    myPrimes.push_back(static_cast<typePrime>(smallPrimeBase[ind]));
            }
            
            std::int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
            std::int_fast64_t myNum = 1 + lowerBnd;
            
            std::size_t p = 1u;
            std::int_fast64_t sqrPrime = 9;
            std::vector<std::int_fast64_t> nextStrt;
            
            if (minNum > 2) {
                for (std::int_fast64_t myIndex = 0; sqrPrime <= upperBnd; ++p) {
                    if (lowerBnd > sqrPrime) {
                        const std::int_fast64_t remTest = lowerBnd % sievePrimes[p - 1u];
                        if (remTest == 0) {
                            myIndex = sievePrimes[p - 1u];
                        } else {
                            myIndex = sievePrimes[p - 1u] - remTest;
                            if ((myIndex % 2) == 0) {myIndex += sievePrimes[p - 1u];}
                        }
                    } else {
                        myIndex = sqrPrime - lowerBnd;
                    }
                    
                    nextStrt.push_back(myIndex);
                    sqrPrime = sievePrimes[p] * sievePrimes[p];
                }
                
                for (std::size_t i = 4u; i < nextStrt.size(); ++i) {
                    std::int_fast64_t j = nextStrt[i];
                    for (std::int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                        sieve[j] = 0;
                    
                    nextStrt[i] = j - segSize;
                }
                
                std::size_t strtChunk = 0u;
                
                if (minNum > myNum) {
                    strtChunk = (minNum - myNum) / NUM2310;
                    myNum += strtChunk * NUM2310;
                    
                    for (std::size_t w = 0u; w < szWheel2310 && myNum <= maxNum; myNum += ARR_WHEEL2310[w], ++w) {
                        if (myNum >= minNum) {
                            if (sieve[myNum - lowerBnd]) {
                                myPrimes.push_back(static_cast<typePrime>(myNum));
                            }
                        }
                    }
                    
                    ++strtChunk;
                }
                
                if (upperBnd < flrMaxNum) {
                    std::int_fast64_t sieveInd = myNum - lowerBnd;
                    
                    for (std::size_t q = strtChunk; q < nWheels; ++q) {
                        for (std::size_t w = 0u; w < szWheel2310; sieveInd += ARR_WHEEL2310[w + 7u], w += 8u) {
                            if (sieve[sieveInd]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                            if (sieve[sieveInd += ARR_WHEEL2310[w]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 1u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 2u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 3u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 4u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 5u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                            if (sieve[sieveInd += ARR_WHEEL2310[w + 6u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        }
                    }
                    
                    myNum = sieveInd + lowerBnd;
                } else {
                    for (std::size_t q = strtChunk; q < nWheels && myNum <= maxNum; ++q)
                        for (std::size_t w = 0u; w < szWheel2310 && myNum <= maxNum; myNum += ARR_WHEEL2310[w], ++w)
                            if (sieve[myNum - lowerBnd])
                                myPrimes.push_back(static_cast<typePrime>(myNum));
                }
                
                std::fill(sieve.begin(), sieve.end(), 1);
                lowerBnd += segSize;
            }
            
            for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
                upperBnd = lowerBnd + segSize;
                
                for (; sqrPrime <= upperBnd; ++p) {
                    nextStrt.push_back(sqrPrime - lowerBnd);
                    sqrPrime = sievePrimes[p] * sievePrimes[p];
                }
                
                for (std::size_t i = 4u; i < nextStrt.size(); ++i) {
                    std::int_fast64_t j = nextStrt[i];
                    for (std::int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                        sieve[j] = 0;
                    
                    nextStrt[i] = j - segSize;
                }
                
                std::int_fast64_t sieveInd = myNum - lowerBnd;
                
                for (std::size_t q = 0u; q < nWheels; ++q) {
                    for (std::size_t w = 0u; w < szWheel2310; sieveInd += ARR_WHEEL2310[w + 7u], w += 8u) {
                        if (sieve[sieveInd]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL2310[w]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL2310[w + 1u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL2310[w + 2u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL2310[w + 3u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL2310[w + 4u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL2310[w + 5u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL2310[w + 6u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    }
                }
                
                myNum = sieveInd + lowerBnd;
                std::fill(sieve.begin(), sieve.end(), 1);
            }
            
            // Get remaining primes that are greater than flrMaxNum and less than maxNum
            if (lowerBnd < maxNum) {
                for (; sqrPrime <= maxNum; ++p) {
                    nextStrt.push_back(sqrPrime - lowerBnd);
                    sqrPrime = sievePrimes[p] * sievePrimes[p];
                }
                
                for (std::size_t i = 4u; i < nextStrt.size(); ++i) {
                    std::int_fast64_t j = nextStrt[i];
                    for (std::int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                        sieve[j] = 0;
                }
                
                for (std::size_t q = 0u; q < nWheels && myNum <= maxNum; ++q)
                    for (std::size_t w = 0u; w < szWheel2310 && myNum <= maxNum; myNum += ARR_WHEEL2310[w], ++w)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back(static_cast<typePrime>(myNum));
            }
        }
    }
    
    template <typename typePrime>
    void PrimeSieveMedium(std::int_fast64_t minNum, std::int_fast64_t maxNum, 
                          const std::vector<std::int_fast64_t> &sievePrimes,
                          std::vector<typePrime> &myPrimes) {
        
        const std::size_t szWheel30030 = SZ_WHEEL30030;
        const std::int_fast64_t sz30030 = NUM30030;
        const std::size_t nWheels = static_cast<std::size_t>(std::max(1.0, std::ceil(std::sqrt(static_cast<double>(maxNum)) / sz30030)));
        const std::int_fast64_t segSize = static_cast<std::int_fast64_t>(nWheels * sz30030);
        
        // We have to add myPrimes.size() because we could have started in PrimeSieveSmall
        const std::size_t myReserve = EstimatePiPrime(static_cast<double>(minNum),
                                                      static_cast<double>(maxNum)) + myPrimes.size();
        myPrimes.reserve(myReserve);
        const std::int_fast64_t flrMaxNum = segSize * std::floor(maxNum / segSize);
        std::vector<bool> sieve(segSize, true);
        
        std::int_fast64_t lowerBnd = segSize * std::floor(static_cast<double>(minNum) / segSize);
        std::int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
        std::int_fast64_t myNum = 1 + lowerBnd;
        
        std::size_t p = 1u;
        std::int_fast64_t sqrPrime = 9;
        std::vector<std::int_fast64_t> nextStrt;
        
        for (std::int_fast64_t myIndex = 0; sqrPrime <= upperBnd; ++p) {
            if (lowerBnd > sqrPrime) {
                std::int_fast64_t remTest = lowerBnd % sievePrimes[p - 1u];
                if (remTest == 0) {
                    myIndex = sievePrimes[p - 1u];
                } else {
                    myIndex = sievePrimes[p - 1u] - remTest;
                    if ((myIndex % 2) == 0) {myIndex += sievePrimes[p - 1u];}
                }
            } else {
                myIndex = sqrPrime - lowerBnd;
            }
            
            nextStrt.push_back(myIndex);
            sqrPrime = sievePrimes[p] * sievePrimes[p];
        }
        
        for (std::size_t i = 5u; i < nextStrt.size(); ++i) {
            std::int_fast64_t j = nextStrt[i];
            for (std::int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                sieve[j] = false;
            
            nextStrt[i] = j - segSize;
        }
        
        std::size_t strtChunk = 0u;
        
        if (minNum > myNum) {
            strtChunk = (minNum - myNum) / sz30030;
            myNum += strtChunk * sz30030;
            
            for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum; myNum += ARR_WHEEL30030[w], ++w) {
                if (myNum >= minNum) {
                    if (sieve[myNum - lowerBnd]) {
                        myPrimes.push_back(static_cast<typePrime>(myNum));
                    }
                }
            }
            
            ++strtChunk;
        }
        
        if (upperBnd < flrMaxNum) {
            std::int_fast64_t sieveInd = myNum - lowerBnd;
            
            for (std::size_t q = strtChunk; q < nWheels; ++q) {
                for (std::size_t w = 0u; w < szWheel30030; sieveInd += ARR_WHEEL30030[w + 7u], w += 8u) {
                    if (sieve[sieveInd]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 1u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 2u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 3u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 4u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 5u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 6u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                }
            }
            
            myNum = sieveInd + lowerBnd;
        } else {
            for (std::size_t q = strtChunk; q < nWheels && myNum <= maxNum; ++q)
                for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum; myNum += ARR_WHEEL30030[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back(static_cast<typePrime>(myNum));
        }
        
        std::fill(sieve.begin(), sieve.end(), true);
        lowerBnd += segSize;
        
        for (; lowerBnd < flrMaxNum; lowerBnd += segSize) {
            upperBnd = lowerBnd + segSize;
            
            for (; sqrPrime <= upperBnd; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = sievePrimes[p] * sievePrimes[p];
            }
            
            for (std::size_t i = 5u; i < nextStrt.size(); ++i) {
                std::int_fast64_t j = nextStrt[i];
                for (std::int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = false;
                
                nextStrt[i] = j - segSize;
            }
            
            std::int_fast64_t sieveInd = myNum - lowerBnd;
            
            for (std::size_t q = 0u; q < nWheels; ++q) {
                for (std::size_t w = 0u; w < szWheel30030; sieveInd += ARR_WHEEL30030[w + 7u], w += 8u) {
                    if (sieve[sieveInd]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 1u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 2u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 3u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 4u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 5u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 6u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                }
            }
            
            myNum = sieveInd + lowerBnd;
            std::fill(sieve.begin(), sieve.end(), true);
        }
        
        // Get remaining primes that are greater than flrMaxNum and less than maxNum
        if (lowerBnd < maxNum) {
            for (; sqrPrime <= maxNum; ++p) {
                nextStrt.push_back(sqrPrime - lowerBnd);
                sqrPrime = sievePrimes[p] * sievePrimes[p];
            }
            
            for (std::size_t i = 5u; i < nextStrt.size(); ++i) {
                std::int_fast64_t j = nextStrt[i];
                for (std::int_fast64_t k = sievePrimes[i] * 2; j < segSize; j += k)
                    sieve[j] = false;
            }
            
            for (std::size_t q = 0u; q < nWheels && myNum <= maxNum; ++q)
                for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum; myNum += ARR_WHEEL30030[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back(static_cast<typePrime>(myNum));
        }
    }
    
    template <typename typePrime>
    void PrimeSieveBig(std::int_fast64_t minNum, std::int_fast64_t maxNum, 
                       const std::vector<std::int_fast64_t> &svPriOne,
                       const std::vector<std::int_fast64_t> &svPriTwo, 
                       std::vector<typePrime> &myPrimes, std::size_t nBigSegs,
                       const std::vector<char> &check30030) {
        
        const std::int_fast64_t sz30030 = NUM30030;
        const std::int_fast64_t segSize = static_cast<std::int_fast64_t>(nBigSegs * sz30030);
        const std::size_t numWheelSegs = nBigSegs, szWheel30030 = SZ_WHEEL30030;
        
        // We have to add myPrimes.size() because we could have started in PrimeSieveMedium
        const std::size_t myReserve = EstimatePiPrime(static_cast<double>(minNum),
                                                      static_cast<double>(maxNum)) + myPrimes.size();;
        myPrimes.reserve(myReserve);
        std::int_fast64_t remTest, divTest, myIndex, lowerBnd;
        lowerBnd = static_cast<std::int_fast64_t>(segSize * (minNum / segSize));
        
        const std::size_t svPriOneSize = svPriOne.size();
        std::vector<std::int_fast64_t> nextStrtOne(svPriOneSize);
        
        for (std::size_t i = 0u; i < svPriOneSize; ++i) {
            remTest = lowerBnd % svPriOne[i];
            
            if (remTest == 0) {
                myIndex = svPriOne[i];
            } else {
                myIndex = svPriOne[i] - remTest;
                if (myIndex % 2 == 0) {myIndex += svPriOne[i];}
            }
            
            nextStrtOne[i] = myIndex;
        }
        
        const std::int_fast64_t myRange = (maxNum - lowerBnd) + 1;
        std::int_fast64_t numCacheSegs = myRange / segSize;
        
        double wholeTest = static_cast<double>(myRange) / segSize;
        
        if (wholeTest != static_cast<int>(wholeTest))
            ++numCacheSegs;
        
        std::int_fast64_t remPrime, timesTwo;
        std::int_fast64_t tempInd, maxIndex = myRange + 1;
        bool bKeepGoing;
        
        // Keeps track of which primes will be used in each interval
        std::deque<std::vector<std::size_t>> myBuckets(numCacheSegs,
                                                       std::vector<std::size_t>());
        
        for (std::size_t i = 0u; i < svPriTwo.size(); ++i) {
            remTest = (lowerBnd % svPriTwo[i]);
            
            if (remTest == 0) {
                myIndex = svPriTwo[i];
            } else {
                myIndex = (svPriTwo[i] - remTest);
                if ((myIndex % 2) == 0) {myIndex += svPriTwo[i];}
            }
            
            remTest = (myIndex % sz30030) - 1;
            timesTwo = (2 * svPriTwo[i]);
            remPrime = (timesTwo % sz30030);
            bKeepGoing = (myIndex <= maxIndex);
            
            // Populate rest of the buckets
            while (bKeepGoing) {
                while (check30030[remTest] && bKeepGoing) {
                    myIndex += timesTwo;
                    bKeepGoing = (myIndex <= maxIndex);
                    remTest += remPrime;
                    if (remTest >= sz30030) {remTest -= sz30030;}
                }
                
                tempInd = myIndex;
                
                if (bKeepGoing) {
                    divTest = (myIndex / segSize);
                    myIndex -= (divTest * segSize);
                    myBuckets[divTest].push_back(static_cast<std::size_t>(myIndex));
                }
                
                myIndex = tempInd + timesTwo;
                bKeepGoing = (myIndex <= maxIndex);
                remTest += remPrime;
                if (remTest >= sz30030) remTest -= sz30030;
            }
        }
        
        const std::int_fast64_t flrMaxNum = segSize * std::floor(static_cast<double>(maxNum) / segSize);
        std::int_fast64_t upperBnd = std::min(lowerBnd + segSize, maxNum);
        std::int_fast64_t myNum = 1 + lowerBnd;
        
        // vector used for sieving
        std::vector<bool> sieve(segSize, true);
        std::size_t strt = 0u;
        
        if (minNum != lowerBnd) {
            for (std::size_t i = 5u; i < svPriOneSize; ++i) {
                std::int_fast64_t j = nextStrtOne[i];
                for (std::int_fast64_t k = (svPriOne[i] * 2); j < segSize; j += k)
                    sieve[j] = false;
                
                nextStrtOne[i] = (j - segSize);
            }
            
            for (const auto &elem: myBuckets[0])
                sieve[elem] = false;
            
            myBuckets.pop_front();
            std::size_t strtChunk = 0u;
            
            if (minNum > myNum) {
                strtChunk = (minNum - myNum) / NUM30030;
                myNum += strtChunk * NUM30030;
                
                for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum; myNum += ARR_WHEEL30030[w], ++w)
                    if (myNum >= minNum)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back(static_cast<typePrime>(myNum));
                        
                ++strtChunk;
            }
            
            if (upperBnd < flrMaxNum) {
                std::int_fast64_t sieveInd = myNum - lowerBnd;
                
                for (std::size_t q = strtChunk; q < numWheelSegs; ++q) {
                    for (std::size_t w = 0u; w < szWheel30030; sieveInd += ARR_WHEEL30030[w + 7u], w += 8u) {
                        if (sieve[sieveInd]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL30030[w]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL30030[w + 1u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL30030[w + 2u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL30030[w + 3u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL30030[w + 4u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL30030[w + 5u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                        if (sieve[sieveInd += ARR_WHEEL30030[w + 6u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    }
                }
                
                myNum = sieveInd + lowerBnd;
            } else {
                for (std::size_t q = strtChunk; q < numWheelSegs && myNum <= maxNum; ++q)
                    for (std::size_t w = 0u; w < szWheel30030 && myNum <= maxNum; myNum += ARR_WHEEL30030[w], ++w)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back(static_cast<typePrime>(myNum));
            }
            
            std::fill(sieve.begin(), sieve.end(), true);
            lowerBnd += segSize;
            ++strt;
        }
        
        for (std::size_t v = strt; lowerBnd < flrMaxNum; ++v, lowerBnd += segSize) {
            for (std::size_t i = 5u; i < svPriOneSize; ++i) {
                std::int_fast64_t j = nextStrtOne[i];
                for (std::int_fast64_t k = (svPriOne[i] * 2); j < segSize; j += k)
                    sieve[j] = false;
                
                nextStrtOne[i] = (j - segSize);
            }
            
            for (const auto &elem: myBuckets[0])
                sieve[elem] = false;
            
            myBuckets.pop_front();
            std::int_fast64_t sieveInd = myNum - lowerBnd;
            
            for (std::size_t q = 0u; q < numWheelSegs; ++q) {
                for (std::size_t w = 0u; w < szWheel30030; sieveInd += ARR_WHEEL30030[w + 7u], w += 8u) {
                    if (sieve[sieveInd]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 1u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 2u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 3u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 4u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 5u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                    if (sieve[sieveInd += ARR_WHEEL30030[w + 6u]]) myPrimes.push_back(static_cast<typePrime>(sieveInd + lowerBnd));
                }
            }
            
            myNum = sieveInd + lowerBnd;
            std::fill(sieve.begin(), sieve.end(), true);
        }
        
        // Get remaining primes that are greater than flrMaxNum and less than maxNum
        if (lowerBnd < maxNum) {
            for (std::size_t i = 5u; i < svPriOneSize; ++i) {
                std::int_fast64_t j = nextStrtOne[i];
                for (std::int_fast64_t k = (svPriOne[i] * 2); j < segSize; j += k)
                    sieve[j] = false;
            }
            
            for (const auto &elem: myBuckets[0])
                sieve[elem] = false;
            
            for (std::size_t q = 0u; (myNum <= maxNum) && (q < numWheelSegs); ++q)
                for (std::size_t w = 0u; (myNum <= maxNum) && (w < szWheel30030); myNum += ARR_WHEEL30030[w], ++w)
                    if (sieve[myNum - lowerBnd])
                        myPrimes.push_back(static_cast<typePrime>(myNum));
        }
    }
    
    inline void sqrtSmallPrimes(int sqrtBound, std::vector<std::int_fast64_t> &sievePrimes) {
        std::size_t ind = 1u;
        for (; smallPrimeBase[ind] <= sqrtBound; ++ind)
            sievePrimes.push_back(smallPrimeBase[ind]);
        
        sievePrimes.push_back(smallPrimeBase[ind]);
    }
    
    template <typename typePrime>
    void sqrtBigPrimes(int sqrtBound, bool bAddZero, bool bAddExtraPrime,
                       bool bAddTwo, std::vector<typePrime> &sievePrimes) {
        
        if (sqrtBound < smallPrimeBase.back()) {
            if (bAddZero) sievePrimes.push_back(0);
            std::size_t ind = (bAddTwo) ? 0u : 1u;
            
            for (; smallPrimeBase[ind] <= sqrtBound; ++ind)
                sievePrimes.push_back(smallPrimeBase[ind]);
            
            if (bAddExtraPrime)
                sievePrimes.push_back(smallPrimeBase[ind]);
        } else {
            int sqrtSqrtBound = static_cast<int>(std::sqrt(static_cast<double>(sqrtBound)));
            std::vector<std::int_fast64_t> sqrtSievePrimes;
            sqrtSmallPrimes(sqrtSqrtBound, sqrtSievePrimes);
            
            // The number, 225, comes from the observation that the largest prime
            // gap less than 100 million is 219 @ 47,326,693. This is important
            // because we need to guarantee that we obtain the smallest prime
            // greater than sqrt(n). We know that the largest number that this
            // algo can except is 2^53 - 1, which gives a sqrt of 94,906,266
            
            std::int_fast64_t myLower = 3, myUpper = sqrtBound;
            if (bAddExtraPrime) {myUpper += 225;}
            const std::size_t sqrtReserve = EstimatePiPrime(1.0, static_cast<double>(myUpper));
            sievePrimes.reserve(sqrtReserve);
            
            if (bAddZero) {sievePrimes.push_back(0);}
            if (bAddTwo) {myLower = 1;}
            PrimeSieveSmall(myLower, myUpper, sqrtSievePrimes, sievePrimes);
        }
    }
    
    template <typename typePrime>
    void PrimeWorker(std::int_fast64_t minNum, std::int_fast64_t maxNum, const std::vector<std::int_fast64_t> &svPriMain,
                     const std::vector<std::int_fast64_t> &svPriOne, const std::vector<std::int_fast64_t> &svPriTwo, 
                     std::vector<typePrime> &primes, std::int_fast64_t smallCut, std::int_fast64_t medCut, 
                     std::size_t nBigSegs, const std::vector<char> &check30030) {
        
        if (maxNum > medCut) {
            // PrimeSieveBig is not meant for small values, so we
            // first get them with PrimeSieveMed/Small
            if (minNum < smallCut) {
                PrimeSieveSmall(minNum, smallCut, svPriMain, primes);
                PrimeSieveMedium(smallCut, medCut, svPriMain, primes);
                PrimeSieveBig(medCut, maxNum, svPriOne, svPriTwo, primes, nBigSegs, check30030);
            } else if (minNum < medCut) {
                PrimeSieveMedium(minNum, medCut, svPriMain, primes);
                PrimeSieveBig(medCut, maxNum, svPriOne, svPriTwo, primes, nBigSegs, check30030);
            } else {
                PrimeSieveBig(minNum, maxNum, svPriOne, svPriTwo, primes, nBigSegs, check30030);
            }
        } else if (maxNum > smallCut) {
            if (minNum < smallCut) {
                PrimeSieveSmall(minNum, smallCut, svPriMain, primes);
                PrimeSieveMedium(smallCut, maxNum, svPriMain, primes);
            } else {
                PrimeSieveMedium(minNum, maxNum, svPriMain, primes);
            }
        } else {
            PrimeSieveSmall(minNum, maxNum, svPriMain, primes);
        }
    }
    
    template <typename typePrime>
    void PrimeSieveMaster(std::int_fast64_t minNum, std::int_fast64_t maxNum, std::vector<typePrime> &primes, 
                          std::vector<std::vector<typePrime>> &primeList, bool &Parallel,
                          int nThreads = 1, int maxThreads = 1, int maxCores = 1) {
        
        const std::int_fast64_t myRange = maxNum - minNum;
        std::int_fast64_t smallCut = Almost2310L1Cache - 100;
        smallCut *= smallCut;
        
        // Here, we estimate the maximum sieving interval. Based off of empirical
        // evidence we start seeing improved performance in PrimeSieveBig after
        // (sqrtMedCut)^2 vs PrimeSieveMedium. We suspect this has to do with the
        // amount of L3 Cache available. On our main testing device, we have 6MiB.
        // If you were to calculate the sieving size using the formula below, we
        // would obtain (4 * 1024^2 * 6) / 8 (std::vector<bool> possible memory
        // optimization of packing elements into individual bits) ~= 3MiB which
        // gives us room to store other important objects for sieving.
        std::int_fast64_t sqrtMedCut = maxCores * (1024 * 1024) * 6;
        
        const int sqrtBound = std::sqrt(static_cast<double>(maxNum));
        std::size_t nBigSegs = 1u;
        
        if (nThreads > 1 && maxThreads > 1) {
            Parallel = true;
            if (nThreads > maxThreads) {nThreads = maxThreads;}
            sqrtMedCut /= nThreads;
        }
        
        const std::int_fast64_t medCut = sqrtMedCut * sqrtMedCut;
        std::vector<std::int_fast64_t> svPriMain, svPriOne, svPriTwo;
        std::vector<char> check30030(NUM30030, 1);
        std::int_fast64_t segUnitSize = NUM30030;
        
        //       bool bAddZero, bool bAddExtraPrime, bool bAddTwo
        sqrtBigPrimes(sqrtBound, false, true, false, svPriMain);
        
        if (maxNum > medCut) {
            
            // Similar to sqrtMedCut above (see comment). This formulation
            // was determined empirically. The idea is that as the number
            // of threads increases, we want to ensure that we don't pollute
            // the cache memory. Each thread will need enough space for all
            // of the data structures necessary to produce primes, thus
            // we reduce the segment size as nThreads increases.
            nBigSegs = ((maxThreads / nThreads) * 2 * 1024 * 1024) /  L1_CACHE_SIZE;
            
            std::size_t ind = 0u;
            const int limitOne = static_cast<int>(nBigSegs * segUnitSize);
            const std::size_t svMainSize = svPriMain.size();
            
            // Get the primes that are guaranteed to mark an index in every segment interval
            for (; (2 * svPriMain[ind]) < limitOne && ind < svMainSize; ++ind)
                svPriOne.push_back(svPriMain[ind]);
            
            // Get the rest
            for (; ind < svMainSize; ++ind)
                svPriTwo.push_back(svPriMain[ind]);
            
            for (std::size_t i = 0u, ind = 0u; i < SZ_WHEEL30030; ind += ARR_WHEEL30030[i], ++i)
                check30030[ind] = 0;
            
        } else if (maxNum > smallCut) {
            nBigSegs = static_cast<std::size_t>(std::ceil(
                static_cast<double>(sqrtBound) / segUnitSize));
        } else {
            segUnitSize = NUM2310;
            nBigSegs = N_WHEELS2310_PER_SEG;
        }
        
        nBigSegs = (nBigSegs < 1) ? 1u : nBigSegs;
        const std::int_fast64_t segSize = nBigSegs * segUnitSize;
        
        // There is some overhead for setting up the data structures
        // for each thread, so we ensure each thread has enough work
        // to outweigh adding another thread. Determined empirically.
        // The tests that were performed showed that one thread was
        // just as efficient as two threads when myRange was equal
        // to 4 * segSize. This is especially true for PrimeSieveBig.
        const std::int_fast64_t testNumThreads = myRange / (segSize * 4);
        
        if (maxThreads < 2) {Parallel = false;}
        
        if (Parallel) {
            if (testNumThreads > 1) {
                if (testNumThreads < nThreads) {nThreads = testNumThreads;}
            } else {
                Parallel = false;
            }
        }
        
        if (Parallel) {
            int ind = 0;
            RcppThread::ThreadPool pool(nThreads);
            std::int_fast64_t upperBnd, lowerBnd = minNum;
            const std::int_fast64_t chunkSize = segSize * std::floor(static_cast<double>(myRange / nThreads) / segSize);
            upperBnd = segSize * std::floor(minNum / segSize) + chunkSize;
            
            for (; ind < (nThreads - 1); lowerBnd = upperBnd, upperBnd += chunkSize, ++ind) {
                pool.push(std::cref(PrimeWorker<typePrime>), lowerBnd, upperBnd, std::ref(svPriMain), std::ref(svPriOne),
                          std::ref(svPriTwo), std::ref(primeList[ind]), smallCut, medCut, nBigSegs, std::ref(check30030));
            }
            
            pool.push(std::cref(PrimeWorker<typePrime>), lowerBnd, maxNum, std::ref(svPriMain), std::ref(svPriOne),
                      std::ref(svPriTwo), std::ref(primeList[ind]), smallCut, medCut, nBigSegs, std::ref(check30030));
            
            pool.join();
        } else {
            PrimeWorker(minNum, maxNum, svPriMain, svPriOne, svPriTwo, primes, smallCut, medCut, nBigSegs, check30030);
        }
    }
}

template void PrimeSieve::PrimeSieveMaster<double>(std::int_fast64_t, std::int_fast64_t, std::vector<double>&, 
                                                   std::vector<std::vector<double>>&, bool&, int, int, int);

template void PrimeSieve::PrimeSieveMaster<std::int64_t>(std::int_fast64_t, std::int_fast64_t, std::vector<std::int64_t>&, 
                                                    std::vector<std::vector<std::int64_t>>&, bool&, int, int, int);

template void PrimeSieve::PrimeSieveMaster<int>(std::int_fast64_t, std::int_fast64_t, std::vector<int>&, 
                                                std::vector<std::vector<int>>&, bool&, int, int, int);

template void PrimeSieve::sqrtBigPrimes<double>(int, bool, bool, bool, std::vector<double>&);
template void PrimeSieve::sqrtBigPrimes<std::int64_t>(int, bool, bool, bool, std::vector<std::int64_t>&);
template void PrimeSieve::sqrtBigPrimes<int>(int, bool, bool, bool, std::vector<int>&);
