#include <RcppAlgos.h>
using namespace Rcpp;

SEXP PermutationsRcpp(int n, int m, bool repetition, CharacterVector vStr,
                      int nRows, std::vector<int> vInt, std::vector<double> vNum,
                      bool isMult, bool isFac, bool keepRes, bool isChar, SEXP Rv,
                      bool isInt, std::vector<int> myReps, SEXP f1, SEXP f2) {
    
    if (isChar) {
        if (isMult) {
            return MultisetPermutation<CharacterMatrix>(n, m, vStr, myReps, nRows, false);
        } else {
            return PermuteGeneral<CharacterMatrix>(n, m, vStr, repetition, nRows, false);
        }
    } else {
        if (Rf_isNull(f1))
            keepRes = false;
        
        if (keepRes) {
            NumericMatrix matRes;
            
            std::string mainFun2 = as<std::string >(f1);
            if (mainFun2 != "prod" && mainFun2 != "sum" && mainFun2 != "mean"
                    && mainFun2 != "max" && mainFun2 != "min") {
                stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
            }
            
            std::vector<double> rowVec(m);
            XPtr<funcPtr> xpFun2 = putFunPtrInXPtr(mainFun2);
            funcPtr myFun2 = *xpFun2;
            
            if (isMult) {
                matRes = MultisetPermutation<NumericMatrix>(n, m, vNum, myReps, nRows, true);
            } else {
                matRes = PermuteGeneral<NumericMatrix>(n, m, vNum, repetition, nRows, true);
            }
            
            for (std::size_t i = 0; i < nRows; i++) {
                for (std::size_t j = 0; j < m; j++)
                    rowVec[j] = matRes(i, j);
                
                matRes(i, m) = myFun2(rowVec);
            }
            
            return matRes;
        } else {
            if (isFac) {
                IntegerMatrix factorMat;
                IntegerVector testFactor = as<IntegerVector>(Rv);
                CharacterVector myClass = testFactor.attr("class");
                CharacterVector myLevels = testFactor.attr("levels");
                
                if (isMult) {
                    factorMat = MultisetPermutation<IntegerMatrix>(n, m, vInt, myReps, nRows, false);
                } else {
                    factorMat = PermuteGeneral<IntegerMatrix>(n, m, vInt, repetition, nRows, false);
                }
                
                factorMat.attr("class") = myClass;
                factorMat.attr("levels") = myLevels;
                
                return factorMat;
            } else {
                if (isInt) {
                    if (isMult) {
                        return MultisetPermutation<IntegerMatrix>(n, m, vInt, myReps, nRows, false);
                    } else {
                        return PermuteGeneral<IntegerMatrix>(n, m, vInt, repetition, nRows, false);
                    }
                } else {
                    if (isMult) {
                        return MultisetPermutation<NumericMatrix>(n, m, vNum, myReps, nRows, false);
                    } else {
                        return PermuteGeneral<NumericMatrix>(n, m, vNum, repetition, nRows, false);
                    }
                }
            }
        }
    }
}