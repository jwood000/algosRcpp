#include "FunAssign.h"

// See function do_vapply found here:
//     https://github.com/wch/r-source/blob/trunk/src/main/apply.c
void VapplyAssign(SEXP ans, SEXP vectorPass, 
                  SEXP sexpFun, SEXP rho, int commonType,
                  int commonLen, int count, int nRows) {
    
    SEXPTYPE valType;
    PROTECT_INDEX indx;
    SETCADR(sexpFun, vectorPass);
    
    SEXP val = Rf_eval(sexpFun, rho);
    PROTECT_WITH_INDEX(val, &indx);
    
    if (Rf_length(val) != commonLen) {
        Rf_error("values must be length %d,\n but FUN(X[[%d]]) result is length %d",
                 commonLen, count + 1, Rf_length(val));
    }
    
    valType = TYPEOF(val);
    if (valType != commonType) {
        bool okay = false;
        switch (commonType) {
            case CPLXSXP: okay = (valType == REALSXP) || (valType == INTSXP)
                || (valType == LGLSXP); break;
            case REALSXP: okay = (valType == INTSXP) || (valType == LGLSXP); break;
            case INTSXP:  okay = (valType == LGLSXP); break;
        }
        if (!okay)
            Rf_error("values must be type '%s',\n but FUN(X[[%d]]) result is type '%s'",
                  Rf_type2char(commonType), count + 1, Rf_type2char(valType));
        REPROTECT(val = Rf_coerceVector(val, commonType), indx);
    }
    
    if(commonLen == 1) { // common case
        switch (commonType) {
            case CPLXSXP: COMPLEX(ans)[count] = COMPLEX(val)[0]; break;
            case REALSXP: REAL(ans)   [count] = REAL   (val)[0]; break;
            case INTSXP:  INTEGER(ans)[count] = INTEGER(val)[0]; break;
            case LGLSXP:  LOGICAL(ans)[count] = LOGICAL(val)[0]; break;
            case RAWSXP:  RAW(ans)    [count] = RAW    (val)[0]; break;
            case STRSXP:  SET_STRING_ELT(ans, count, STRING_ELT(val, 0)); break;
            case VECSXP:  SET_VECTOR_ELT(ans, count, VECTOR_ELT(val, 0)); break;
        }
    } else { // commonLen > 1 (typically, or == 0) :
        switch (commonType) {
            case REALSXP: {
                double* ans_dbl = REAL(ans);
                double* val_dbl = REAL(val);
                
                for (int j = 0; j < commonLen; j++)
                    ans_dbl[count + j * nRows] = val_dbl[j];
                
                break;
            } case INTSXP: {
                int* ans_int = INTEGER(ans);
                int* val_int = INTEGER(val);
                
                for (int j = 0; j < commonLen; j++)
                    ans_int[count + j * nRows] = val_int[j];
                
                break;
            } case LGLSXP: {
                int* ans_bool = LOGICAL(ans);
                int* val_bool = LOGICAL(val);
                
                for (int j = 0; j < commonLen; j++)
                    ans_bool[count + j * nRows] = val_bool[j];
                
                break;
            } case RAWSXP: {
                Rbyte* ans_raw = RAW(ans);
                Rbyte* val_raw = RAW(val);
                
                for (int j = 0; j < commonLen; j++)
                    ans_raw[count + j * nRows] = val_raw[j];
                
                break;
            } case CPLXSXP: {
                Rcomplex* ans_cmplx = COMPLEX(ans);
                Rcomplex* val_cmplx = COMPLEX(val);
                
                for (int j = 0; j < commonLen; j++)
                    ans_cmplx[count + j * nRows] = val_cmplx[j];
                
                break;
            } case STRSXP: {
                for (int j = 0; j < commonLen; j++)
                    SET_STRING_ELT(ans, count + j * nRows, STRING_ELT(val, j));
                
                break;
            } case VECSXP: {
                for (int j = 0; j < commonLen; j++)
                    SET_VECTOR_ELT(ans, count + j * nRows, VECTOR_ELT(val, j));
                
                break;
            }
        }
    }
    
    UNPROTECT(1);
}

void FunAssign(SEXP res, SEXP vectorPass, SEXP sexpFun,
               SEXP rho, int commonType, int commonLen,
               int count, int nRows, int retType) {
    
    if (retType == VECSXP) {
        SETCADR(sexpFun, vectorPass);
        SET_VECTOR_ELT(res, count, Rf_eval(sexpFun, rho));
    } else {
        VapplyAssign(res, vectorPass, sexpFun, rho,
                     commonType, commonLen, count, nRows);
    }
    
}