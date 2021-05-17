#include <RcppThread.h>
#include <libdivide.h>
#include "CleanConvert.h"
#include <cmath>

template <typename typeInt>
inline typeInt getStartingIndex (typeInt lowerB, typeInt step) {
    
    if (step >= lowerB)
        return (2 * step - lowerB);
        
    const typeInt remTest = lowerB % step;
    const typeInt retStrt = (remTest == 0) ? 0 : (step - remTest);

    return retStrt;
}

template <typename typeInt, typename typeReturn>
void NumDivisorsSieve(typeInt m, typeInt n, typeInt offsetStrt, typeReturn &numFacs) {
    
    const typeInt myRange = offsetStrt + (n - m) + 1;
    const typeInt sqrtBound = static_cast<typeInt>(std::sqrt(n));
    
    for (typeInt i = 2; i <= sqrtBound; ++i) {
        const typeInt myNum = offsetStrt + (i * sqrtBound) - m;
        typeInt j = offsetStrt + getStartingIndex(m, i);
        
        for (; j <= myNum; j += i)
            ++numFacs[j];
        
        for (; j < myRange; j += i)
            numFacs[j] += 2;
    }
    
    // Subtract 1 from the first entry as 1 has only itself
    // as a divisor. N.B. myRange was initialized with 2
    if (m < 2)
        --numFacs[0];
}

template <typename typeInt, typename typeReturn>
void DivisorsSieve(typeInt m, typeReturn retN, typeInt offsetStrt,
                   std::vector<std::vector<typeReturn>> &MyDivList) {
    
    const typeInt n = static_cast<typeInt>(retN);
    constexpr typeInt zeroOffset = 0;
    const typeInt myRange = (n - m) + 1;
    
    typename std::vector<std::vector<typeReturn>>::iterator it2d, itEnd;
    itEnd = MyDivList.begin() + offsetStrt + myRange;
    
    std::vector<int> myMemory(myRange, 2);
    NumDivisorsSieve(m, n, zeroOffset, myMemory);
    
    if (m < 2)
        it2d = MyDivList.begin() + 1;
    else
        it2d = MyDivList.begin() + offsetStrt;
    
    if (m < 2) {
        for (std::size_t i = 1; it2d < itEnd; ++it2d, ++i) {
            it2d->reserve(myMemory[i]);
            it2d->push_back(1);
        }
        
        MyDivList[0].push_back(1);
        
        for (typeInt i = 2; i <= n; ++i)
            for (typeInt j = i; j <= n; j += i)
                MyDivList[j - 1].push_back(static_cast<typeReturn>(i));
        
    } else {
        typeReturn numRet = static_cast<typeReturn>(m);
        std::vector<int> begIndex(myRange, 0);
        
        for (std::size_t i = 0; it2d < itEnd; ++it2d, ++i, ++numRet) {
            it2d->resize(myMemory[i]);
            it2d->back() = numRet;
            it2d->front() = 1;
            --myMemory[i];
        }
        
        typeInt sqrtBound = static_cast<typeInt>(std::sqrt(n));
        typeInt offsetRange = myRange + offsetStrt;

        for (typeInt i = 2; i <= sqrtBound; ++i) {

            const typeInt myStart = getStartingIndex(m, i);
            const libdivide::divider<typeInt> fastDiv(i);

            for (typeInt j = myStart + offsetStrt, myNum = m + myStart,
                     memInd = myStart; j < offsetRange; j += i, myNum += i, memInd += i) {

                MyDivList[j][++begIndex[memInd]] = static_cast<typeReturn>(i);
                const typeInt testNum = myNum / fastDiv;

                // Ensure we won't duplicate adding an element. If
                // testNum <= sqrtBound, it will be added in later
                // iterations. Also, we insert this element in the
                // pentultimate position as it will be the second
                // to the largest element at the time of inclusion.
                // E.g. let i = 5, myNum = 100, so the current
                // vectors looks like so: v = 1, 5, 10, 100 (5 was
                // added to the second position above). With i = 5,
                // testNum = 100 / 5 = 20, thus we add it to the
                // pentultimate position to give v = 1 5 10 20 100.
                if (testNum > sqrtBound)
                    MyDivList[j][--myMemory[memInd]] = static_cast<typeReturn>(testNum);
            }
        }
    }
}

template <typename typeInt, typename typeReturn, typename typeDivCount>
void DivisorMaster(typeInt myMin, typeReturn myMax, bool bDivSieve,
                   typeDivCount &DivCountV, std::vector<std::vector<typeReturn>> &MyDivList,
                   std::size_t myRange, int nThreads = 1, int maxThreads = 1) {
    
    bool Parallel = false;
    typeInt offsetStrt = 0;
    const typeInt intMax = static_cast<typeInt>(myMax);
    
    if (nThreads > 1 && maxThreads > 1  && myRange >= 20000) {
        Parallel = true;
        
        if (nThreads > maxThreads)
            nThreads = maxThreads;
        
        // Ensure that each thread has at least 10000
        if ((myRange / nThreads) < 10000)
            nThreads = myRange / 10000;
    }
    
    if (Parallel) {
        RcppThread::ThreadPool pool(nThreads);
        typeInt lowerBnd = myMin;
        const typeInt chunkSize = myRange / nThreads;
        typeReturn upperBnd = lowerBnd + chunkSize - 1;
        
        for (int ind = 0; ind < (nThreads - 1); offsetStrt += chunkSize, 
                    lowerBnd = (upperBnd + 1), upperBnd += chunkSize, ++ind) {
            if (bDivSieve) {
                pool.push(std::cref(DivisorsSieve<typeInt, typeReturn>), 
                          lowerBnd, upperBnd, offsetStrt, std::ref(MyDivList));
            } else {
                pool.push(std::cref(NumDivisorsSieve<typeInt, typeDivCount>), lowerBnd,
                          static_cast<typeInt>(upperBnd), offsetStrt, std::ref(DivCountV));
            }
        }

        if (bDivSieve) {
            pool.push(std::cref(DivisorsSieve<typeInt, typeReturn>),
                      lowerBnd, myMax, offsetStrt, std::ref(MyDivList));
        } else {
            pool.push(std::cref(NumDivisorsSieve<typeInt, typeDivCount>),
                      lowerBnd, intMax, offsetStrt, std::ref(DivCountV));
        }

        pool.join();
        
    } else {
        if (bDivSieve) {
            DivisorsSieve(myMin, myMax, offsetStrt, MyDivList);
        } else {
            NumDivisorsSieve(myMin, intMax, offsetStrt, DivCountV);
        }     
    }
}

template <typename typeInt, typename typeReturn>
SEXP TheGlue(typeInt myMin, typeReturn myMax, bool bDivSieve,
             bool keepNames, int nThreads, int maxThreads) {
    
    const std::size_t myRange = (myMax - myMin) + 1;
    std::vector<typeReturn> myNames;
    
    if (keepNames) {
        myNames.resize(myRange);
        typeReturn retM = myMin;
        for (std::size_t k = 0; retM <= myMax; ++retM, ++k)
            myNames[k] = retM;
    }
    
    if (bDivSieve) {
        std::vector<std::vector<typeReturn>> 
            MyDivList(myRange, std::vector<typeReturn>());
        Rcpp::IntegerVector tempRcpp;
        DivisorMaster(myMin, myMax, bDivSieve, tempRcpp,
                      MyDivList, myRange, nThreads, maxThreads);
        
        Rcpp::List myList = Rcpp::wrap(MyDivList);
        if (keepNames)
            myList.attr("names") = myNames;
        
        return myList;
    } else {
        std::vector<std::vector<typeReturn>> tempList;
        Rcpp::IntegerVector facCountV(myRange, 2);
        DivisorMaster(myMin, myMax, bDivSieve, facCountV,
                      tempList, myRange, nThreads, maxThreads);
        
        if (keepNames)
            facCountV.attr("names") = myNames;
        
        return facCountV;
    }
}

// [[Rcpp::export]]
SEXP DivNumSieve(SEXP Rb1, SEXP Rb2, bool bDivSieve, 
                 SEXP RNamed, SEXP RNumThreads, int maxThreads) {
    
    double bound1, myMax, myMin;
    const std::string namedObject = (bDivSieve) ? "namedList" : "namedVector";
    bool keepNames = CleanConvert::convertLogical(RNamed, namedObject);
    CleanConvert::convertPrimitive(Rb1, bound1, "bound1");
    
    double bound2 = 1;
    if (!Rf_isNull(Rb2)) CleanConvert::convertPrimitive(Rb2, bound2, "bound2");
    
    if (bound1 > bound2) {
        myMax = std::floor(bound1);
        myMin = std::ceil(bound2);
    } else {
        myMax = std::floor(bound2);
        myMin = std::ceil(bound1);
    }
    
    if (myMax < 2) {
        if (bDivSieve) {
            std::vector<std::vector<int> > trivialRet(1, std::vector<int>(1, 1));
            Rcpp::List z = Rcpp::wrap(trivialRet);
            if (keepNames)
                z.attr("names") = 1;
            
            return z;
        } else {
            Rcpp::IntegerVector v(1, 1);
            if (keepNames)
                v.attr("names") = 1;
            
            return v;
        }
    }
    
    int nThreads = 1;
    if (!Rf_isNull(RNumThreads))
        CleanConvert::convertPrimitive(RNumThreads, nThreads, "nThreads");
    
    if (myMax > std::numeric_limits<int>::max()) {
        std::int_fast64_t intMin = static_cast<std::int_fast64_t>(myMin);
        return TheGlue(intMin, myMax, bDivSieve, keepNames, nThreads, maxThreads);
    } else {
        int intMin = static_cast<int>(myMin);
        int intMax = static_cast<int>(myMax);
        return TheGlue(intMin, intMax, bDivSieve, keepNames, nThreads, maxThreads);
    }
}
