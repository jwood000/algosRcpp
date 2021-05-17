#include "PollardRhoUtils.h"

template <typename typeReturn>
void PrimeFacList(std::size_t m, std::size_t n, std::vector<double> myNums,
                  std::vector<std::vector<typeReturn>> &MyPrimeList) {
    
    for (std::size_t i = m; i < n; ++i) {
        std::vector<typeReturn> factors;
        
        std::int64_t mPass = static_cast<std::int64_t>(myNums[i]);
        
        if (mPass < 0) {
            mPass = std::abs(mPass);
            factors.push_back(-1);
        }
        
        if (mPass > 0) {
            getPrimeFactors(mPass, factors);
            MyPrimeList[i] = factors;
        }
    }
}

template <typename typeReturn>
std::vector<typeReturn> Factorize(std::vector<typeReturn> &factors) {
    
    std::size_t n = factors.size();
    
    if (n == 1) {
        std::vector<typeReturn> primeReturn(2, 1);
        primeReturn[1] = factors[0];
        return primeReturn;
    } else {
        std::vector<std::size_t> lengths;
        typeReturn prev = factors[0];
        
        std::size_t numUni = 0;
        std::vector<typeReturn> uniFacs(n);
        uniFacs[0] = factors[0];
        lengths.push_back(1);
        
        for(auto it = factors.cbegin() + 1; it < factors.cend(); ++it) {
            if (prev == *it) {
                ++lengths[numUni];
            } else {
                ++numUni;
                prev = *it;
                lengths.push_back(1);
                uniFacs[numUni] = *it;
            }
        }
        
        std::size_t numFacs = 1;
        
        for (std::size_t i = 0; i <= numUni; ++i)
            numFacs *= (lengths[i] + 1);
        
        std::vector<typeReturn> myFacs(numFacs);
        
        for (std::size_t i = 0; i <= lengths[0]; ++i)
            myFacs[i] = static_cast<typeReturn>(std::pow(uniFacs[0], i));
        
        if (numUni > 0) {
            std::size_t fSz = 1;
            typeReturn temp;
            
            for (std::size_t j = 1; j <= numUni; ++j) {
                fSz *= (lengths[j - 1] + 1);
                for (std::size_t i = 1; i <= lengths[j]; ++i) {
                    for (std::size_t k = 0, ind = (i * fSz); k < fSz; ++k, ++ind) {
                        temp = static_cast<typeReturn>(std::pow(uniFacs[j], i));
                        temp *= myFacs[k];
                        myFacs[ind] = temp;
                    }
                }
            }
        }
        
        std::sort(myFacs.begin(), myFacs.end());
        return myFacs;
    }
}

template <typename typeReturn>
void FactorList(std::size_t m, std::size_t n, std::vector<double> &myNums,
                std::vector<std::vector<typeReturn>> &MyDivList) {
    
    bool isNegative = false;
    
    for (std::size_t j = m; j < n; ++j) {
        std::vector<typeReturn> myDivisors;
        std::int64_t mPass = static_cast<std::int64_t>(myNums[j]);
        
        if (mPass < 0) {
            mPass = std::abs(mPass);
            isNegative = true;
        } else {
            isNegative = false;
        }
        
        if (mPass > 1) {
            std::vector<typeReturn> factors;
            getPrimeFactors(mPass, factors);
            myDivisors = Factorize<typeReturn>(factors);
            
            if (isNegative) {
                const std::size_t facSize = myDivisors.size();
                std::vector<typeReturn> tempInt(2 * facSize);
                std::size_t posInd = facSize, negInd = facSize - 1;
                
                for (std::size_t i = 0; i < facSize; ++i, ++posInd, --negInd) {
                    tempInt[negInd] = -1 * myDivisors[i];
                    tempInt[posInd] = myDivisors[i];
                }
                
                myDivisors = tempInt;
            }
        } else {
            if (isNegative)
                myDivisors.push_back(-1);
            if (mPass > 0)
                myDivisors.push_back(1);
        }
        
        MyDivList[j] = myDivisors;
    }
}

void IsPrimeVec(std::size_t m, std::size_t n, std::vector<double> &myNums,
                Rcpp::LogicalVector &primeTest) {
    
    mpz_t testMpzt;
    mpz_init(testMpzt);
    
    for (std::size_t j = m; j < n; ++j) {
        std::int64_t testVal = static_cast<std::int64_t>(myNums[j]);
        
        if (testVal == 1) {
            primeTest[j] = false;
        } else if ((testVal & 1) == 0) {
            if (testVal > 2)
                primeTest[j] = false;
        } else {
            int p = 3;
            for (std::size_t i = 1; i < primesDiffPR.size();) {
                if ((testVal % p) != 0) {
                    p += primesDiffPR[i++];
                    if (testVal < (p * p))
                        break;
                } else {
                    if (testVal > p)
                        primeTest[j] = false;
                    break;
                }
            }
        }
        
        if (primeTest[j]) {
            // 1e9 was determined empirically. Creating an mpz_t and calling
            // mpz_probab_prime_p isn't as efficient as calling a similar
            // function that only deals with primitive types.
            if (myNums[j] < 1000000000) {
                primeTest[j] = IsPrime(testVal);
            } else {
                mpz_set_d(testMpzt, myNums[j]);
                
                if (mpz_probab_prime_p(testMpzt, MR_REPS) == 0)
                    primeTest[j] = false;
            }
        }
    }
    
    mpz_clear(testMpzt);
}


template void PrimeFacList(std::size_t, std::size_t, std::vector<double>,
                           std::vector<std::vector<int>>&);

template void PrimeFacList(std::size_t, std::size_t, std::vector<double>,
                           std::vector<std::vector<double>>&);

template std::vector<int> Factorize(std::vector<int> &factors);
template std::vector<double> Factorize(std::vector<double> &factors);

template void FactorList(std::size_t, std::size_t, std::vector<double>&,
                         std::vector<std::vector<int>>&);

template void FactorList(std::size_t, std::size_t, std::vector<double>&,
                         std::vector<std::vector<double>>&);

