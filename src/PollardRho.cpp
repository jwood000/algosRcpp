#include "PollardRhoDepends.h"
#include "PollardRhoUtils.h"
#include "CleanConvert.h"
#include <RcppThread.h>
#include <cmath>

template <typename typeReturn>
void PollardRhoMaster(std::vector<double> &myNums, typeReturn myMax, bool bPrimeFacs,
                      bool bAllFacs, std::vector<std::vector<typeReturn>> &MyList,
                      Rcpp::LogicalVector &primeTest, std::size_t myRange, 
                      int nThreads = 1, int maxThreads = 1) {
    
    bool Parallel = false;
    
    if (nThreads > 1 && myRange > 1 && maxThreads > 1) {
        Parallel = true;
        if (nThreads > maxThreads) {nThreads = maxThreads;}
        if ((myRange / nThreads) < 1) {nThreads = myRange;}
    }
    
    if (Parallel) {
        RcppThread::ThreadPool pool(nThreads);
        const std::size_t chunkSize = myRange / nThreads;
        std::size_t n = chunkSize - 1;
        std::size_t m = 0;
        
        for (int j = 0; j < (nThreads - 1); m = n, n += chunkSize, ++j) {
            if (bPrimeFacs)
                pool.push(std::cref(PrimeFacList<typeReturn>), m, n, std::ref(myNums), std::ref(MyList));
            else if (bAllFacs)
                pool.push(std::cref(FactorList<typeReturn>), m, n, std::ref(myNums), std::ref(MyList));
            else
                pool.push(std::cref(IsPrimeVec), m, n, std::ref(myNums), std::ref(primeTest));
        }
        
        if (bPrimeFacs)
            pool.push(std::cref(PrimeFacList<typeReturn>), m, myRange, std::ref(myNums), std::ref(MyList));
        else if (bAllFacs)
            pool.push(std::cref(FactorList<typeReturn>), m, myRange, std::ref(myNums), std::ref(MyList));
        else
            pool.push(std::cref(IsPrimeVec), m, myRange, std::ref(myNums), std::ref(primeTest));
        
        pool.join();
        
    } else {
        if (bPrimeFacs)
            PrimeFacList(0u, myRange, myNums, MyList);
        else if (bAllFacs)
            FactorList(0u, myRange, myNums, MyList);
        else
            IsPrimeVec(0u, myRange, myNums, primeTest);
    }
}

template <typename typeReturn>
SEXP TheGlue(std::vector<double> &myNums, typeReturn myMax, bool bPrimeFacs,
             bool bAllFacs, bool keepNames, int nThreads, int maxThreads) {
    
    std::size_t myRange = myNums.size();
    Rcpp::LogicalVector tempVec;
    
    if (bPrimeFacs) {
        if (myRange == 1) {
            std::int64_t mPass = static_cast<std::int64_t>(myNums[0]);
            if (mPass == 0) {return Rcpp::IntegerVector();}
            std::vector<typeReturn> factors;
            
            if (mPass < 0) {
                mPass = std::abs(mPass);
                factors.push_back(-1);
            }
            
            getPrimeFactors(mPass, factors);
            return Rcpp::wrap(factors);
        } else {
            std::vector<std::vector<typeReturn>> 
                MyPrimeList(myRange, std::vector<typeReturn>());
            
            PollardRhoMaster(myNums, myMax, bPrimeFacs, bAllFacs,
                             MyPrimeList, tempVec, myRange, nThreads, maxThreads);
            
            Rcpp::List myList = Rcpp::wrap(MyPrimeList);
            if (keepNames)
                myList.attr("names") = myNums;
            
            return myList;
        }
    } else if (bAllFacs) {
        if (myRange == 1) {
            std::int64_t mPass = static_cast<std::int64_t>(myNums[0]);
            std::vector<typeReturn> factors;
            std::vector<typeReturn> myDivisors;
            bool isNegative = false;
            
            if (mPass < 0) {
                mPass = std::abs(mPass);
                isNegative = true;
            }
            
            if (mPass > 1) {
                getPrimeFactors(mPass, factors);
                myDivisors = Factorize<typeReturn>(factors);
                
                if (isNegative) {
                    const std::size_t facSize = myDivisors.size();
                    std::vector<typeReturn> negPosFacs(2 * facSize);
                    std::size_t posInd = facSize, negInd = facSize - 1;
                    
                    for (std::size_t i = 0; i < facSize; ++i, ++posInd, --negInd) {
                        negPosFacs[negInd] = -1 * myDivisors[i];
                        negPosFacs[posInd] = myDivisors[i];
                    }
                    
                    return Rcpp::wrap(negPosFacs);
                }
                
                return Rcpp::wrap(myDivisors);
            } else  {
                if (isNegative)
                    myDivisors.push_back(-1);
                if (mPass > 0)
                    myDivisors.push_back(1);
                
                return Rcpp::wrap(myDivisors);
            }
        }
        
        std::vector<std::vector<typeReturn>> 
            MyPrimeList(myRange, std::vector<typeReturn>());
        
        PollardRhoMaster(myNums, myMax, bPrimeFacs, bAllFacs,
                         MyPrimeList, tempVec, myRange, nThreads, maxThreads);
        
        Rcpp::List myList = Rcpp::wrap(MyPrimeList);
        if (keepNames)
            myList.attr("names") = myNums;
        
        return myList;
    } else {
        Rcpp::LogicalVector isPrimeVec(myRange, true);
        std::vector<std::vector<typeReturn>> tempList;
        
        PollardRhoMaster(myNums, myMax, bPrimeFacs, bAllFacs,
                         tempList, isPrimeVec, myRange, nThreads, maxThreads);
        if (keepNames)
            isPrimeVec.attr("names") = myNums;
        
        return isPrimeVec;
    }
}

// [[Rcpp::export]]
SEXP PollardRhoContainer(SEXP Rv, SEXP RNamed, bool bPrimeFacs,
                         bool bAllFacs, SEXP RNumThreads, int maxThreads) {
    
    std::vector<double> myNums;
    bool isNamed = CleanConvert::convertLogical(RNamed, "namedList");
    
    if (bPrimeFacs || bAllFacs)  // numOnly = true, checkWhole = true, negPoss = true
        CleanConvert::convertVector(Rv, myNums, "v", true, true, true);
    else
        CleanConvert::convertVector(Rv, myNums, "v");
    
    double myMax = *std::max_element(myNums.cbegin(), myNums.cend());
    double myMin = *std::min_element(myNums.cbegin(), myNums.cend());
    
    if (std::abs(myMin) > myMax)
        myMax = std::abs(myMin);
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads");
    
    if (myMax > std::numeric_limits<int>::max()) {
        return TheGlue(myNums, myMax, bPrimeFacs, bAllFacs, isNamed, nThreads, maxThreads);
    } else {
        int intMax = static_cast<int>(myMax);
        return TheGlue(myNums, intMax, bPrimeFacs, bAllFacs, isNamed, nThreads, maxThreads);
    }
}

