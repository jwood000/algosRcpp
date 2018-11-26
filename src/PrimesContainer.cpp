#include <Rcpp.h>
#include <cmath>
#include <array>
#include <libdivide.h>
#include <PrimesSegSieve.h>
#include <PhiTinyLookup.h>
#include <GetFacsUtils.h>

// "AllPrimesCpp" implements a simple segmented version of the
// Sieve of Eratosthenes (original implementation authored by 
// Kim Walisch). An overview of this method can be found here: 
//                      http://primesieve.org/segmented_sieve.html
// Kim Walisch's official github repo for the primesieve is:
//                      https://github.com/kimwalisch/primesieve
// "MasterPrimeCount" is a simple variation of the highly optimized
// "pi_legendre.cpp" algorithm by Kim Walisch, which calculates
// the numbers of primes less than n using Legendre's formula.
// Kim Walisch's official github repo for "pi_legendre" is:
//                      https://github.com/kimwalisch/primecount

// This is the largest multiple of 2*3*5*7 = 210
// that is less than 2^15 = 32768 = 32KB. This
// is the typical size of most CPU's L1 cache
const int L1CacheSize = 32760;
const unsigned int wheelSize = 48;
const int maxVal210 = 210;

static const int_fast64_t wheel210[wheelSize] = {10, 2, 4, 2, 4, 6, 2, 6, 4, 2,
                                                  4, 6, 6, 2, 6, 4, 2, 6, 4, 6,
                                                  8, 4, 2, 4, 2, 4, 8, 6, 4, 6,
                                                  2, 4, 6, 2, 6, 6, 4, 2, 4, 6,
                                                  2, 6, 4, 2, 4, 2, 10, 2};

constexpr unsigned long int firstPriSize = (sizeof(firstPrimes) / sizeof(firstPrimes[0]));

template <typename typeReturn, typename typePrime>
std::vector<typeReturn> AllPrimesCpp (typePrime minNum,
                                     typePrime maxNum,
                                     bool bAddZero = false) {
    // bAddZero is used for counting primes. More
    // specifically, it makes it easier to access a
    // particular index and avoids countless
    // subtractions (i.e. to access the 2nd, simply
    // use 2 as your index instead of 2 - 1).

    int sqrtBound = (int) std::sqrt((double) maxNum);
    int_fast64_t segmentSize = (int_fast64_t) L1CacheSize;
    unsigned int numSegments = L1CacheSize / maxVal210;
    std::vector<typeReturn> myPrimes;
    typePrime myRange = maxNum - minNum + 1;
    
    // Percentages obtained here: https://en.wikipedia.org/wiki/Prime-counting_function
    typePrime myReserve;
    if (maxNum < 100000) {
        myReserve = (typePrime) floor((double) 2*myRange/log((double)myRange));
    } else if (maxNum < 10000000) {
        myReserve = (typePrime) floor((double) myRange/(log((double)myRange) * 0.905));
    } else if (maxNum < 100000000) {
        myReserve = (typePrime) floor((double) myRange/(log((double)myRange) * 0.933));
    } else if (maxNum < 1000000000) {
        myReserve = (typePrime) floor((double) myRange/(log((double)myRange) * 0.942));
    } else {
        myReserve = (typePrime) floor((double) myRange/(log((double)myRange) * 0.948));
    }
    
    myPrimes.reserve(myReserve);
    if (bAddZero) myPrimes.push_back(0);
    std::size_t i = 0;

    if (maxNum < firstPrimes[firstPriSize - 1] || minNum <= 13)
        for (; firstPrimes[i] < minNum; ++i) {}

    if (maxNum < firstPrimes[firstPriSize - 1]) {
        for (; firstPrimes[i] <= maxNum; ++i) 
            myPrimes.push_back((typeReturn) firstPrimes[i]);
        
        return myPrimes;
    } else if (minNum <= 13) {
        for (; firstPrimes[i] < 10; ++i)
            myPrimes.push_back((typeReturn) firstPrimes[i]);
    }

    // ensure segmentSize is greater than sqrt(n)
    // as well as multiple of L1CacheSize
    if (segmentSize < sqrtBound) {
        numSegments = (unsigned int) ceil((double) sqrtBound / maxVal210);
        segmentSize = numSegments * maxVal210;
    }
    
    std::vector<int_fast64_t> smallPrimes, nextStrt;
    int_fast64_t flrMaxNum = segmentSize * floor((double) maxNum / segmentSize);

    // Create vector of primes up to the sqrt(n)
    // discarding 2, and tacking on the next
    // prime after the sqrt(n)
    if (maxNum <= 1000000000) {
        int i = 1;
        for (; firstPrimes[i] <= sqrtBound; ++i)
            smallPrimes.push_back(firstPrimes[i]);
        
        smallPrimes.push_back(firstPrimes[i]);
    } else {
        for (std::size_t i = 1; i < (firstPriSize - 1); ++i)
            smallPrimes.push_back(firstPrimes[i]);

        // The number, 225, comes from the observation that
        // the largest prime gap less than 100 million is
        // 219 @ 47,326,693. This is important because we
        // need to guarantee that we obtain the smallest
        // prime greater than sqrt(n). We know that the
        // largest number that this algorithm can except
        // is 2^53 - 1, which gives a sqrt of 94,906,266
        int_fast64_t myLimit = sqrtBound + 225;
        std::vector<char> sqrtIsPrime(myLimit + 1, 1);
        int_fast64_t step, lastP = smallPrimes[0];

        for (std::size_t p = 0; lastP * lastP <= myLimit; ++p) {
            lastP = smallPrimes[p];
            step = lastP * 2;
            for (int_fast64_t j = lastP * lastP; j <= myLimit; j += step)
                sqrtIsPrime[j] = 0;
        }

        unsigned int currentSize = smallPrimes.size();
        lastP = smallPrimes[currentSize - 1];

        for (int_fast64_t j = lastP + 2; j <= myLimit; j += 2)
            if (sqrtIsPrime[j])
                smallPrimes.push_back(j);
    }
    
    // vector used for sieving
    std::vector<char> sieve(segmentSize, 1);
    sieve[1] = 0;

    int_fast64_t sqrPrime = (int_fast64_t) (3 * 3);
    int_fast64_t lowerBnd = 0, upperBnd = segmentSize, myNum = 1;
    int p = 1;

    if (minNum > 2) {
        lowerBnd = (typePrime) (segmentSize * floor((double) minNum / segmentSize));
        upperBnd = std::min(lowerBnd + segmentSize, (int_fast64_t) maxNum);
        myNum += lowerBnd;
        int_fast64_t myStart;

        for (; sqrPrime <= upperBnd; ++p) {
            if (lowerBnd > sqrPrime) {
                int_fast64_t remTest = lowerBnd % smallPrimes[p - 1];
                if (remTest == 0) {
                    myStart = 0;
                } else {
                    myStart = smallPrimes[p - 1] - remTest;
                }
                if ((myStart % 2) == 0)
                    myStart += smallPrimes[p - 1];
            } else {
                myStart = sqrPrime - lowerBnd;
            }
            
            nextStrt.push_back(myStart);
            sqrPrime = smallPrimes[p] * smallPrimes[p];
        }

        for (std::size_t i = 3; i < nextStrt.size(); ++i) {
            int_fast64_t j = nextStrt[i];
            for (int_fast64_t k = smallPrimes[i] * 2; j < segmentSize; j += k)
                sieve[j] = 0;
            nextStrt[i] = j - segmentSize;
        }

        if (upperBnd < flrMaxNum) {
            for (std::size_t q = 0; q < numSegments; ++q) {
                for (std::size_t w = 0; w < wheelSize; ++w) {
                    if (myNum >= minNum)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back((typeReturn) myNum);
                    myNum += wheel210[w];
                }
            }
        } else {
            for (std::size_t q = 0; q < numSegments && myNum <= maxNum; ++q) {
                for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                    if (myNum >= minNum)
                        if (sieve[myNum - lowerBnd])
                            myPrimes.push_back((typeReturn) myNum);
                    myNum += wheel210[w];
                }
            }
        }

        std::fill(sieve.begin(), sieve.end(), 1);
        lowerBnd += segmentSize;
    }

    for (; lowerBnd < flrMaxNum; lowerBnd += segmentSize) {
        // current segment = interval [lowerBnd, upperBnd]
        upperBnd = lowerBnd + segmentSize;

        // sieve the current segment && sieving primes <= sqrt(upperBnd)
        for (; sqrPrime <= upperBnd; ++p) {
            nextStrt.push_back(sqrPrime - lowerBnd);
            sqrPrime = smallPrimes[p] * smallPrimes[p];
        }

        for (std::size_t i = 3; i < nextStrt.size(); ++i) {
            int_fast64_t j = nextStrt[i];
            for (int_fast64_t k = smallPrimes[i] * 2; j < segmentSize; j += k)
                sieve[j] = 0;
            nextStrt[i] = j - segmentSize;
        }

        for (std::size_t q = 0; q < numSegments; ++q) {
            for (std::size_t w = 0; w < wheelSize; ++w) {
                if (sieve[myNum - lowerBnd])
                    myPrimes.push_back((typeReturn) myNum);
                myNum += wheel210[w];
            }
        }

        std::fill(sieve.begin(), sieve.end(), 1);
    }

    // Get remaining primes that are greater than flrMaxNum and less than maxNum
    if (lowerBnd < maxNum) {
        for (; sqrPrime <= maxNum; ++p) {
            nextStrt.push_back(sqrPrime - lowerBnd);
            sqrPrime = smallPrimes[p] * smallPrimes[p];
        }

        for (std::size_t i = 3; i < nextStrt.size(); ++i) {
            int_fast64_t j = nextStrt[i];
            for (int_fast64_t k = smallPrimes[i] * 2; j < segmentSize; j += k)
                sieve[j] = 0;
            nextStrt[i] = j - segmentSize;
        }

        for (std::size_t q = 0; q < numSegments && myNum <= maxNum; ++q) {
            for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                if (sieve[myNum - lowerBnd])
                    myPrimes.push_back((typeReturn) myNum);
                myNum += wheel210[w];
            }
        }
    }

    return myPrimes;
}

// PiPrime is very similar to the algorithm above
// only we are not considering a range. That is,
// we are only concerned with finding the number
// of primes less than maxNum. We are also only
// counting primes instead of generating them.
int64_t PiPrime (int64_t maxNum) {
    
    int sqrtBound = (int) std::sqrt((double) maxNum);
    int64_t segmentSize = (int64_t) L1CacheSize;
    unsigned int numSegments = L1CacheSize / maxVal210;
    
    // the wheel already has the first 4 primes marked as
    // false, so we need to account for them here. N.B.
    // the calling functions has checks for cases where
    // maxNum is less than 11, so need to check here.
    int64_t count = 4;
    
    if (segmentSize < sqrtBound) {
        numSegments = (unsigned int) ceil((double) sqrtBound / maxVal210);
        segmentSize = numSegments * maxVal210;
    }
    
    std::vector<int64_t> smallPrimes, nextStrt;
    int64_t flrMaxNum = segmentSize * floor((double) maxNum / segmentSize);
    
    int i = 1;
    for (; firstPrimes[i] <= sqrtBound; ++i)
        smallPrimes.push_back(firstPrimes[i]);
    
    smallPrimes.push_back(firstPrimes[i]);
    std::vector<char> sieve(segmentSize, 1);
    sieve[1] = 0;
    
    int64_t sqrPrime = (int64_t) (3 * 3);
    int64_t lowerBnd = 0, upperBnd = segmentSize, myNum = 1;
    int p = 1;
    
    for (; lowerBnd < flrMaxNum; lowerBnd += segmentSize) {
        upperBnd = lowerBnd + segmentSize;
        
        for (; sqrPrime <= upperBnd; ++p) {
            nextStrt.push_back(sqrPrime - lowerBnd);
            sqrPrime = smallPrimes[p] * smallPrimes[p];
        }
        
        for (std::size_t i = 3; i < nextStrt.size(); ++i) {
            int64_t j = nextStrt[i];
            for (int64_t k = smallPrimes[i] * 2; j < segmentSize; j += k)
                sieve[j] = 0;
            nextStrt[i] = j - segmentSize;
        }
        
        for (std::size_t q = 0; q < numSegments; ++q) {
            for (std::size_t w = 0; w < wheelSize; ++w) {
                if (sieve[myNum - lowerBnd])
                    ++count;
                myNum += wheel210[w];
            }
        }
        
        std::fill(sieve.begin(), sieve.end(), 1);
    }
    
    if (lowerBnd < maxNum) {
        for (; sqrPrime <= maxNum; ++p) {
            nextStrt.push_back(sqrPrime - lowerBnd);
            sqrPrime = smallPrimes[p] * smallPrimes[p];
        }
        
        for (std::size_t i = 3; i < nextStrt.size(); ++i) {
            int64_t j = nextStrt[i];
            for (int64_t k = smallPrimes[i] * 2; j < segmentSize; j += k)
                sieve[j] = 0;
            nextStrt[i] = j - segmentSize;
        }
        
        for (std::size_t q = 0; q < numSegments && myNum <= maxNum; ++q) {
            for (std::size_t w = 0; w < wheelSize && myNum <= maxNum; ++w) {
                if (sieve[myNum - lowerBnd])
                    ++count;
                myNum += wheel210[w];
            }
        }
    }
    
    return count;
}

const int MAX_A = 100;
const int phiTinySize = 6;
std::array<std::vector<uint16_t>, MAX_A> phiCache;

// Increment MAX_A, so we can have easier access
// to indexing.  E.g when a = 3 (i.e. the number)
// of primes is 3), instead of accessing the third
// entry in phiTiny like phiTiny[a - 1], we simply
// use "a" as our index (i.e. phiTiny[a]). The
// first entry in phiTiny will be empty. The same
// goes for phiPrimes and phiPi.
std::array<std::vector<int16_t>, phiTinySize + 1> phiTiny;
std::vector<int64_t> phiPrimes;
std::vector<int64_t> phiPi;

// The arrays below are utilized by phiTiny
// primeProds[n] = \prod_{i=1}^{n} primes[i]
// myTotients[n] = \prod_{i=1}^{n} (primes[i] - 1)
const std::array<int, 7> myTinyPrimes = {{0, 2, 3, 5, 7, 11, 13}};
const std::array<int, 7> primeProds = {{1, 2, 6, 30, 210, 2310, 30030}};
const std::array<int, 7> myTotients = {{1, 1, 2, 8, 48, 480, 5760}};
const std::array<int, 13> myTinyPi = {{0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5}};

inline int64_t phiTinyCalc (int64_t x, int64_t a) {
    int64_t pp = primeProds[a];
    return (x / pp) * myTotients[a] + phiTiny[a][x % pp];
}

void updateCache (uint64_t x, uint64_t a, int64_t mySum) {
    if (a < phiCache.size() &&
        x <= std::numeric_limits<uint16_t>::max()) {
        if (x >= phiCache[a].size())
            phiCache[a].resize(x + 1, 0);
        
        phiCache[a][x] = (uint16_t) std::abs(mySum);
    }
}

inline bool isCached (uint64_t x, uint64_t a) {
    return a < phiCache.size() &&
           x < phiCache[a].size() &&
           phiCache[a][x];
}

inline bool isPix(int64_t x, int64_t a) {
    return x < (int) phiPi.size() &&
           x < (phiPrimes[a + 1] * phiPrimes[a + 1]);
}

inline int64_t getStrt(int64_t y) {
    if (y >= myTinyPrimes.back())
        return phiTinySize;
    else
        return myTinyPi[y];
}

template <int SIGN>
int64_t phiWorker (int64_t x, int64_t a) {
    if (x <= phiPrimes[a])
        return SIGN;
    else if (a <= phiTinySize)
        return phiTinyCalc(x, a) * SIGN;
    else if (isPix(x, a))
        return (phiPi[x] - a + 1) * SIGN;
    else if (isCached(x, a))
        return phiCache[a][x] * SIGN;
        
    int64_t sqrtx = (int64_t) std::sqrt((double) x);
    int64_t piSqrtx = a;
    int64_t strt = getStrt(sqrtx);
    int64_t mySum = 0;
    
    if (sqrtx < (int) phiPi.size())
        piSqrtx = std::min((int64_t) phiPi[sqrtx], a);
    
    mySum += (piSqrtx - a) * SIGN;
    mySum += phiTinyCalc(x, strt) * SIGN;
    
    for (int64_t i = strt; i < piSqrtx; ++i) {
        int64_t x2 = x / phiPrimes[i + 1];
        
        if (isPix(x2, i))
            mySum += (phiPi[x2] - i + 1) * -SIGN;
        else
            mySum += phiWorker<-SIGN>(x2, i);
    }
    
    updateCache(x, a, mySum);
    return mySum;
}

int64_t phiMaster (int64_t x, int64_t a) {
    
    int64_t sqrtx = (int64_t) std::sqrt((double) x);
    int64_t piSqrtx = std::min((int64_t) phiPi[sqrtx], a);
    int64_t strt = getStrt(sqrtx);
    int64_t mySum = phiTinyCalc(x, strt) + piSqrtx - a;
    
    for (int64_t i = strt; i < piSqrtx; ++i)
        mySum += phiWorker<-1>(x / phiPrimes[i + 1], i);
    
    return mySum;
}

// 10^9  -->> 50,847,534
// 10^10 -->> 455,052,511
// 10^11 -->> 4,118,054,813
// 10^12 -->> 37,607,912,018
// 10^13 -->> 346,065,536,839
// 10^14 -->> 3,204,941,750,802
// 10^15 -->> 29,844,570,422,669

//[[Rcpp::export]]
SEXP MasterPrimeCount (SEXP Rn) {
    double bound1;
    
    switch(TYPEOF(Rn)) {
        case REALSXP: {
            bound1 = Rcpp::as<double>(Rn);
            break;
        }
        case INTSXP: {
            bound1 = Rcpp::as<double>(Rn);
            break;
        }
        default: {
            Rcpp::stop("bound1 must be of type numeric or integer");
        }
    }

    if (bound1 < 1 || bound1 > Significand53)
        Rcpp::stop("n must be a positive number less than 2^53");

    int64_t n = (int64_t) bound1;
    if (n < 100000) {
        if (n < 10) {
            if (n == 1) return Rcpp::wrap(0);
            std::vector<int> smallPs = AllPrimesCpp<int>((int) 1, (int) n, false);
            return (Rcpp::wrap(smallPs.size()));
        }
        return Rcpp::wrap((int) PiPrime(n));
    }
    
    // populate phiTiny... phiTiny[0] should not be accessed
    phiTiny[1].resize(1 + 1);
    phiTiny[1][0] = 0;
    phiTiny[1][1] = 1;
    int16_t arr6[] = {0,1,1,1,1,2};
    std::vector<int16_t> phi6 (arr6, arr6 + sizeof(arr6) / sizeof(arr6[0]) );
    phiTiny[2] = phi6;
    std::vector<int16_t> phi30 (arr30, arr30 + sizeof(arr30) / sizeof(arr30[0]) );
    phiTiny[3] = phi30;
    std::vector<int16_t> phi210 (arr210, arr210 + sizeof(arr210) / sizeof(arr210[0]) );
    phiTiny[4] = phi210;
    std::vector<int16_t> phi2310 (arr2310, arr2310 + sizeof(arr2310) / sizeof(arr2310[0]) );
    phiTiny[5] = phi2310;
    
    std::vector<int16_t> phi30030;
    phi30030.reserve(30031);
    phi30030.push_back(0);
    
    for (std::size_t i = 1; i <= 5760; ++i)
        for (int16_t j = 0; j < arr30030freq[i]; ++j)
            phi30030.push_back(i);
    phiTiny[6] = phi30030;
    
    int64_t sqrtBound = (int64_t) std::sqrt((double) n);
    phiPrimes = AllPrimesCpp<int64_t>((int64_t) 1, sqrtBound, true);
    
    phiPi.resize(sqrtBound + 1);
    int64_t count = 0;
    int64_t maxPrime = phiPrimes[phiPrimes.size() - 1];
    
    for (int64_t i = 1; i <= maxPrime; ++i) {
        if (i >= phiPrimes[count + 1])
            ++count;
        phiPi[i] = count;
    }
    
    for (int64_t i = (maxPrime + 1); i <= sqrtBound; ++i)
        phiPi[i] = count;
    
    int64_t piSqrt = PiPrime(sqrtBound);
    int64_t phiSqrt = phiMaster(n, piSqrt);
    int64_t int64result = piSqrt + phiSqrt - 1;
    
    if (int64result > INT_MAX) {
        double dblResult = (double) int64result;
        return Rcpp::wrap(dblResult);
    }
    
    int intResult = (int) int64result;
    return Rcpp::wrap(intResult);
}

// This function is slightly different than the getStartingIndex
// in the DivisorsContainer.cpp file. The step passed in this
// function is a power of prime and requires an additional check
// (i.e. else if (myPrime < lowerB)).
template <typename typeInt>
inline typeInt getStartIndexPowP (typeInt lowerB,
                                 typeInt step, typeInt myPrime) {
    
    typeInt retStrt, remTest = lowerB % step;
    
    if (remTest == 0) {
        retStrt = 0;
    } else if (myPrime < lowerB) {
        retStrt = step - remTest;
    } else {
        retStrt = step - lowerB;
    }
    
    return retStrt;
}

template <typename typeInt, typename typeReturn>
Rcpp::List PrimeFactorizationSieve (typeInt m, typeReturn retN, bool keepNames) {
    
    typeInt n = (typeInt) retN;
    typeInt myRange = n;
    myRange += (1 - m);
    
    std::vector<std::vector<typeReturn> > 
        MyPrimeList(myRange, std::vector<typeReturn>());
    
    typename std::vector<std::vector<typeReturn> >::iterator it2d, itEnd;
    itEnd = MyPrimeList.end();
    typeInt sqrtBound = (typeInt) floor(sqrt((double)n));
    
    typeInt j, myStep, myStart, myNum = m;
    double myLogN = log((double)n);
    unsigned int mySize, limit;
    
    std::vector<typeReturn> myNames;
    if (keepNames) {
        myNames.resize(myRange);
        typeReturn retM = (typeReturn) m;
        for (std::size_t k = 0; retM <= retN; ++retM, ++k)
            myNames[k] = retM;
    }
    
    if (n > 3) {
        std::vector<typeInt> primes = AllPrimesCpp<typeInt>((typeInt) 1, sqrtBound);
        typename std::vector<typeInt>::iterator p, primesEnd;
        std::vector<int> myMemory(myRange, 1);
        std::vector<int>::iterator myMalloc;
        primesEnd = primes.end();
        
        for (p = primes.begin(); p < primesEnd; ++p) {
            limit = (unsigned long int) trunc(myLogN/log((double)*p));
            if (m < 2) {
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = (typeInt) pow((double) *p, (double) i);
                    for (j = (myStep - 1); j < n; j += myStep)
                        ++myMemory[j];
                }
            } else {
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = (typeInt) pow((double) *p, (double) i);
                    myStart = getStartIndexPowP(m, myStep, *p);
                    for (j = myStart; j < myRange; j += myStep)
                        ++myMemory[j];
                }
            }
        }
            
        if (myNum < 2){
            ++myNum;
            it2d = MyPrimeList.begin() + 1;
            myMalloc = myMemory.begin() + 1;
        } else {
            it2d = MyPrimeList.begin();
            myMalloc = myMemory.begin();
        }
        
        for (; it2d < itEnd; ++it2d, ++myNum, ++myMalloc) {
            it2d->reserve(*myMalloc);
            it2d->push_back((typeReturn) myNum);
        }
    
        if (m < 2) {
            for (p = primes.begin(); p < primesEnd; ++p) {
                limit = (unsigned long int) trunc(myLogN/log((double)*p));
                libdivide::divider<typeInt> fastDiv(*p);
                
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = (typeInt) pow((double) *p, (double) i);
                    
                    for (j = (myStep - 1); j < n; j += myStep) {
                        mySize = MyPrimeList[j].size() - 1;
                        if (MyPrimeList[j][mySize] > *p) {
                            myNum = (typeInt) MyPrimeList[j][mySize];
                            myNum /= fastDiv;
                            MyPrimeList[j][mySize] = (typeReturn) myNum;
                            MyPrimeList[j].insert(MyPrimeList[j].end() - 1, (typeReturn) *p);
                        }
                    }
                }
            }
        } else {
            for (p = primes.begin(); p < primesEnd; ++p) {
                limit = (unsigned long int) trunc(myLogN/log((double)*p));
                libdivide::divider<typeInt> fastDiv(*p);
                
                for (std::size_t i = 1; i <= limit; ++i) {
                    myStep = (typeInt) pow((double)*p, (double) i);
                    myStart = getStartIndexPowP(m, myStep, *p);
    
                    for (j = myStart; j < myRange; j += myStep) {
                        mySize = MyPrimeList[j].size() - 1;
                        if (MyPrimeList[j][mySize] > *p) {
                            myNum = (typeInt) MyPrimeList[j][mySize];
                            myNum /= fastDiv;
                            MyPrimeList[j][mySize] = (typeReturn) myNum;
                            MyPrimeList[j].insert(MyPrimeList[j].end() - 1, (typeReturn) *p);
                        }
                    }
                }
            }
        }
    } else { // edge case where m,n = 2 or 3
        int strt = 0;
        myNum = m;
        if (m == 1) {
            ++strt;
            ++myNum;
        }
        for (int i = strt; i < myRange; ++i, ++myNum)
            MyPrimeList[i].push_back(myNum);
    }
    
    Rcpp::List myList = Rcpp::wrap(MyPrimeList);
    if (keepNames)
        myList.attr("names") = myNames;
    
    return myList;
}

template <typename typeRcpp, typename typeInt, typename typeReturn>
typeRcpp EulerPhiSieveCpp (typeInt m, typeReturn retN, bool keepNames) {
    
    typeInt n = (typeInt) retN;
    typeInt myRange = n;
    myRange += (1 - m);
    
    typeInt j, priTypeInt, myNum = m;
    double myLogN = log((double) n);
    unsigned long int limit;
    std::vector<typeReturn> EulerPhis(myRange);
    std::vector<typeInt> numSeq(myRange);
    
    std::vector<typeReturn> myNames;
    if (keepNames) {
        myNames.resize(myRange);
        typeReturn retM = (typeReturn) m;
        for (std::size_t k = 0; retM <= retN; ++retM, ++k)
            myNames[k] = retM;
    }
    
    typeReturn retM = (typeReturn) m;
    for (std::size_t i = 0; retM <= retN; ++retM, ++i) {
        EulerPhis[i] = retM;
        numSeq[i] = (typeInt) retM;
    }
    
    std::vector<typeInt> primes;
    typename std::vector<typeInt>::iterator p;
    
    if (m < 2) {
        primes  = AllPrimesCpp<typeInt>((typeInt) 1, n);
        
        for (p = primes.begin(); p < primes.end(); ++p) {
            libdivide::divider<typeInt> fastDiv(*p);
            for (j = (*p - 1); j < n; j += *p) {
                myNum = (typeInt) EulerPhis[j];
                myNum /= fastDiv;
                EulerPhis[j] -= (typeReturn) myNum;
            }
        }
    } else if (n > 3) {
        typeInt sqrtBound = floor(sqrt((double) n));
        primes = AllPrimesCpp<typeInt>((typeInt) 1, sqrtBound);
        typeInt myStart, myStep;
    
        for (p = primes.begin(); p < primes.end(); ++p) {
            limit = (unsigned long int) trunc(myLogN / log((double) *p));
            priTypeInt = *p;
            myStart = getStartIndexPowP(m, *p, priTypeInt);
            libdivide::divider<typeInt> fastDiv(priTypeInt);
    
            for (j = myStart; j < myRange; j += *p) {
                numSeq[j] /= fastDiv;
                myNum = (typeInt) EulerPhis[j];
                myNum /= fastDiv;
                EulerPhis[j] -= (typeReturn) myNum;
            }
    
            for (std::size_t i = 2; i <= limit; ++i) {
                myStep = (typeInt) pow((double)*p, i);
                myStart = getStartIndexPowP(m, myStep, priTypeInt);
                
                for (j = myStart; j < myRange; j += myStep)
                    numSeq[j] /= fastDiv;
            }
        }

        for (typeInt i = 0; i < myRange; ++i) {
            if (numSeq[i] > 1) {
                myNum = (typeInt) EulerPhis[i];
                myNum /= numSeq[i];
                EulerPhis[i] -= (typeReturn) myNum;
            }
        }
    } else { // edge case where m,n = 2 or 3
        for (int i = 0; i < myRange; ++i)
            --EulerPhis[i];
    }
    
    typeRcpp myVector = Rcpp::wrap(EulerPhis);
    if (keepNames)
        myVector.attr("names") = myNames;
        
    return myVector;
}

// [[Rcpp::export]]
SEXP EratosthenesRcpp (SEXP Rb1, SEXP Rb2, 
                       SEXP RIsList, SEXP RIsEuler, SEXP RNamed) {
    double bound1, bound2, myMax, myMin;
    bool isList = false, isEuler = false, isNamed = false;

    switch(TYPEOF(Rb1)) {
        case REALSXP: {
            bound1 = Rcpp::as<double>(Rb1);
            break;
        }
        case INTSXP: {
            bound1 = Rcpp::as<double>(Rb1);
            break;
        }
        default: {
            Rcpp::stop("bound1 must be of type numeric or integer");
        }
    }
    
    isList = Rcpp::as<bool>(RIsList);
    isEuler = Rcpp::as<bool>(RIsEuler);
    isNamed = Rcpp::as<bool>(RNamed);
    
    if (bound1 <= 0 || bound1 > Significand53)
        Rcpp::stop("bound1 must be a positive number less than 2^53");

    if (Rf_isNull(Rb2)) {
        myMax = floor(bound1);
        
        if (isList) {
            if (myMax < 2){
                std::vector<std::vector<int> > trivialRet(1, std::vector<int>());
                Rcpp::List z = Rcpp::wrap(trivialRet);
                if (isNamed)
                    z.attr("names") = 1;
                return z;
            }
            
            if (myMax > (INT_MAX - 1))
                return PrimeFactorizationSieve((int64_t) 1, (double) myMax, isNamed);
            
            return PrimeFactorizationSieve((int32_t) 1, (int32_t) myMax, isNamed);
        } else if (isEuler) {
            if (myMax <= 1) {
                Rcpp::IntegerVector z(1, 1);
                if (isNamed)
                    z.attr("names") = 1;
                return z;
            }
            
            if (myMax > (INT_MAX - 1))
                return EulerPhiSieveCpp<Rcpp::NumericVector>((int64_t) 1, (double) myMax, isNamed);
            
            return EulerPhiSieveCpp<Rcpp::IntegerVector>((int32_t) 1, (int32_t) myMax, isNamed);
        } else {
            if (myMax < 2)
                return Rcpp::IntegerVector();

            if (myMax > (INT_MAX - 1)) {
                return Rcpp::wrap(AllPrimesCpp<double>((int_fast64_t) 1,
                                                   (int_fast64_t) myMax));
            }
            return Rcpp::wrap(AllPrimesCpp<int_fast32_t>((int_fast32_t) 1,
                                               (int_fast32_t) myMax));
        }
    } else {
        switch(TYPEOF(Rb2)) {
            case REALSXP: {
                bound2 = Rcpp::as<double>(Rb2);
                break;
            }
            case INTSXP: {
                bound2 = Rcpp::as<double>(Rb2);
                break;
            }
            default: {
                Rcpp::stop("bound2 must be of type numeric or integer");
            }
        }
        
        if (bound2 <= 0 || bound2 > Significand53)
            Rcpp::stop("bound2 must be a positive number less than 2^53");

        if (bound1 > bound2) {
            myMax = bound1;
            myMin = bound2;
        } else {
            myMax = bound2;
            myMin = bound1;
        }
        
        myMin = ceil(myMin);
        myMax = floor(myMax);
        
        if (isList) {
            if (myMax < 2) {
                std::vector<std::vector<int> > trivialRet(1, std::vector<int>());
                Rcpp::List z = Rcpp::wrap(trivialRet);
                if (isNamed)
                    z.attr("names") = 1;
                return z;
            }
            
            if (myMax > (INT_MAX - 1))
                return PrimeFactorizationSieve((int64_t) myMin, (double) myMax, isNamed);
            
            return PrimeFactorizationSieve((int32_t) myMin, (int32_t) myMax, isNamed);
        } else if (isEuler) {
            if (myMax <= 1) {
                Rcpp::IntegerVector z(1, 1);
                if (isNamed)
                    z.attr("names") = 1;
                return z;
            }
            
            if (myMax > (INT_MAX - 1))
                return EulerPhiSieveCpp<Rcpp::NumericVector>((int64_t) myMin, (double) myMax, isNamed);
            
            return EulerPhiSieveCpp<Rcpp::IntegerVector>((int32_t) myMin, (int32_t) myMax, isNamed);
        } else {
            if (myMax <= 1)
                return Rcpp::IntegerVector();
            
            if (myMin <= 2) myMin = 1;
            if (myMin == myMax) {++myMax;}
            
            if (myMax > (INT_MAX - 1)) {
                    return Rcpp::wrap(AllPrimesCpp<double>((int_fast64_t) myMin,
                                                       (int_fast64_t) myMax));
            }
            return Rcpp::wrap(AllPrimesCpp<int_fast32_t>((int_fast32_t) myMin, 
                                               (int_fast32_t) myMax));
        }
    }
}
