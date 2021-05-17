#include "StandardUtils.h"

enum returnType : int {
    constrained = 0,
    standard = 1,
    applyFun = 2
};

bool CheckConstrnd(SEXP f1, SEXP f2, SEXP Rtarget) {
    // No need to check myType, as we have already done 
    // so in CheckStdRet. Same goes for IsFactor. 
    bool result = !Rf_isNull(f1) && !Rf_isNull(f2) && !Rf_isNull(Rtarget);
    
    if (result) {
        if (!Rf_isString(f1))
            Rcpp::stop("constraintFun must be passed as a character");
        
        if (!Rf_isString(f2))
            Rcpp::stop("comparisonFun must be passed as a character");
    }
    
    return result;
}

// [[Rcpp::export]]
int CheckReturn(SEXP Rv, SEXP f1, SEXP f2, SEXP Rtarget,
                bool IsFactor, SEXP RKeepRes, SEXP stdFun) {
    
    int res = returnType::constrained;
    
    if (Rf_isNull(f1)) {
        res = returnType::standard;
    } else {
        if (IsFactor) {
            res = returnType::standard;
        } else {
            VecType myType = VecType::Integer;
            SetType(myType, Rv);
            
            if (myType > VecType::Numeric) {
                res = returnType::standard;
            } else {
                if (!Rf_isNull(f2) && !Rf_isNull(Rtarget)) {
                    res = returnType::constrained;  // This is a constrained output
                } else if (Rf_isNull(f2) && Rf_isNull(Rtarget)) {
                    // This is applying a constrained func only
                    if (Rf_isNull(RKeepRes)) {
                        res = returnType::constrained;
                    } else {
                        bool keepRes = CleanConvert::convertLogical(RKeepRes, "keepResults");
                        
                        if (keepRes) {
                            res = returnType::constrained;
                        } else {
                            res = returnType::standard;
                        }
                    }
                } else {
                    // This means the input is non-sensible... std res = only
                    res = returnType::standard;
                }
            }
        }
    }
    
    if (res) {
        const bool applyFun = !Rf_isNull(stdFun) && !IsFactor;
        
        if (applyFun) {
            if (!Rf_isFunction(stdFun))
                Rcpp::stop("FUN must be a function!");
            
            res = returnType::applyFun;
        }
    }
    
    return res;
}
