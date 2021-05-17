#include "MotleyStartIndex.h"
#include "Eratosthenes.h"
#include <libdivide.h>

namespace MotleyPrimes {

    template <typename typeInt, typename typeReturn, typename typeRcpp>
    void EulerPhiSieve(typeInt m, typeReturn retN, typeInt offsetStrt,
                       const std::vector<typeInt> &primes,
                       std::vector<typeInt> &numSeq,
                       typeRcpp &EulerPhis) {
        
        const typeInt n = static_cast<typeInt>(retN);
        const typeInt myRange = (n - m) + 1;
        
        typeInt myNum = m;
        const double myLogN = std::log(n);
        typeReturn retNum = static_cast<typeReturn>(m);
        
        for (std::size_t i = offsetStrt; retNum <= retN; ++retNum, ++i) {
            EulerPhis[i] = retNum;
            numSeq[i] = static_cast<typeInt>(retNum);
        }
        
        if (m < 2) {
            bool tempPar = false;
            std::vector<typeInt> fullPrimes;
            const std::int_fast64_t intMin = static_cast<std::int_fast64_t>(m);
            const std::int_fast64_t intMax = static_cast<std::int_fast64_t>(retN);
            std::vector<std::vector<typeInt>> tempList;
            PrimeSieve::PrimeSieveMaster(intMin, intMax, fullPrimes, tempList, tempPar);
            typename std::vector<typeInt>::iterator p;
            
            for (p = fullPrimes.begin(); p < fullPrimes.end(); ++p) {
                const libdivide::divider<typeInt> fastDiv(*p);
                for (typeInt j = (*p - 1); j < n; j += *p) {
                    myNum = static_cast<typeInt>(EulerPhis[j]);
                    myNum /= fastDiv;
                    EulerPhis[j] -= static_cast<typeReturn>(myNum);
                }
            }
        } else if (n > 3) {
            typename std::vector<typeInt>::const_iterator p;
            const typeInt sqrtBound = static_cast<typeInt>(std::sqrt(static_cast<double>(retN)));
            const typeInt offsetRange = myRange + offsetStrt;
            
            for (p = primes.cbegin(); (*p) <= sqrtBound; ++p) {
                const std::size_t limit = static_cast<std::size_t>(myLogN / std::log(*p));
                const typeInt myStart = offsetStrt + getStartIndexPowP(m, *p, *p);
                const libdivide::divider<typeInt> fastDiv(*p);
                
                for (typeInt j = myStart; j < offsetRange; j += *p) {
                    numSeq[j] /= fastDiv;
                    myNum = static_cast<typeInt>(EulerPhis[j]);
                    myNum /= fastDiv;
                    EulerPhis[j] -= static_cast<typeReturn>(myNum);
                }
                
                for (std::size_t i = 2; i <= limit; ++i) {
                    const typeInt myStep = static_cast<typeInt>(std::pow(*p, i));
                    const typeInt myStart = offsetStrt + getStartIndexPowP(m, myStep, *p);
                    
                    for (typeInt j = myStart; j < offsetRange; j += myStep)
                        numSeq[j] /= fastDiv;
                }
            }
            
            for (typeInt i = offsetStrt; i < offsetRange; ++i) {
                if (numSeq[i] > 1) {
                    myNum = static_cast<typeInt>(EulerPhis[i]);
                    myNum /= numSeq[i];
                    EulerPhis[i] -= static_cast<typeReturn>(myNum);
                }
            }
            
        } else { // edge case where m,n = 2 or 3
            for (int i = 0; i < myRange; ++i)
                --EulerPhis[i];
        }
    }
}

template void MotleyPrimes::EulerPhiSieve(int, int, int,
                                          const std::vector<int>&,
                                          std::vector<int>&,
                                          Rcpp::IntegerVector&);

template void MotleyPrimes::EulerPhiSieve(std::int64_t, double, std::int64_t,
                                          const std::vector<std::int64_t>&,
                                          std::vector<std::int64_t>&,
                                          Rcpp::NumericVector&);
