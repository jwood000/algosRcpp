#include "Constraints/ConstraintsUtils.h"
#include "Partitions/PartitionsDesign.h"
#include "Partitions/PartitionsCount.h"
#include "CombinatoricsCount.h"
#include "ImportExportMPZ.h"
#include "ComputedCount.h"
#include "SetUpUtils.h"

SEXP GetCount(bool IsGmp, mpz_t computedRowsMpz, double computedRows) {

    if (IsGmp) {
        constexpr std::size_t numb = 8 * intSize;
        const std::size_t sizeNum = intSize *
            (2 + (mpz_sizeinbase(computedRowsMpz, 2) + numb - 1) / numb);
        const std::size_t size = intSize + sizeNum;

        SEXP ans = PROTECT(Rf_allocVector(RAWSXP, size));
        char* rPos = (char*) RAW(ans);
        ((int*) rPos)[0] = 1; // first int is vector-size-header

        // current position in rPos[] (starting after vector-size-header)
        myRaw(&rPos[intSize], computedRowsMpz, sizeNum);
        Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));
        UNPROTECT(1);
        return(ans);
    } else {
        if (computedRows > std::numeric_limits<int>::max()) {
            return Rf_ScalarReal(computedRows);
        } else {
            return Rf_ScalarInteger(static_cast<int>(computedRows));
        }
    }
}

SEXP CombinatoricsCount(SEXP Rv, SEXP Rm, SEXP RisRep,
                        SEXP RFreqs, SEXP RIsComb) {

    int n = 0;
    int m = 0;

    bool IsMult = false;
    VecType myType = VecType::Integer;

    std::vector<double> vNum;
    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;

    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    bool IsComb = CleanConvert::convertLogical(RIsComb, "IsComb");

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum,
              Rv, RFreqs, Rm, n, m, IsMult, IsRep);

    const double computedRows = GetComputedRows(IsMult, IsComb, IsRep,
                                                n, m, Rm, freqs, myReps);
    const bool IsGmp = (computedRows > Significand53);

    mpz_t computedRowsMpz;
    mpz_init(computedRowsMpz);

    if (IsGmp) {
        GetComputedRowMpz(computedRowsMpz, IsMult,
                          IsComb, IsRep, n, m, Rm, freqs, myReps);
    }

    return GetCount(IsGmp, computedRowsMpz, computedRows);
}

SEXP PartitionsCount(SEXP Rtarget, SEXP Rv, SEXP Rm,
                     SEXP RisRep, SEXP RFreqs, SEXP RcompFun,
                     SEXP Rlow, SEXP Rtolerance,
                     SEXP RPartDesign, SEXP Rshow) {
    int n = 0;
    int m = 0;

    bool IsMult = false;
    VecType myType = VecType::Integer;

    std::vector<double> vNum;
    std::vector<int> vInt;
    std::vector<int> myReps;
    std::vector<int> freqs;

    const bool IsConstrained = true;
    const std::string mainFun = "sum";
    bool IsRep = CleanConvert::convertLogical(RisRep, "repetition");
    const bool bDesign = CleanConvert::convertLogical(RPartDesign,
                                                      "PartitionsDesign");

    SetType(myType, Rv);
    SetValues(myType, myReps, freqs, vInt, vNum, Rv,
              RFreqs, Rm, n, m, IsMult, IsRep, IsConstrained);

    // Must be defined inside IsInteger check as targetVals could be
    // outside integer data type range which causes undefined behavior
    std::vector<int> targetIntVals;
    const funcPtr<double> funDbl = GetFuncPtr<double>(mainFun);

    std::vector<std::string> compFunVec;
    std::vector<double> targetVals;

    ConstraintType ctype;
    PartDesign part;

    part.isRep = IsRep;
    part.isMult = IsMult;
    part.mIsNull = Rf_isNull(Rm);

    ConstraintSetup(vNum, myReps, targetVals, targetIntVals,
                    funDbl, part, ctype, n, m, compFunVec,
                    mainFun, myType, Rtarget, RcompFun, Rtolerance,
                    Rlow, IsConstrained, true);

    if (part.ptype != PartitionType::NotPartition) {
        if (bDesign) {
            bool Verbose = CleanConvert::convertLogical(Rshow, "showDetail");
            return GetDesign(part, n, Verbose);
        } else {
            mpz_t mpzFiller;
            mpz_init(mpzFiller);
            return GetCount(false, mpzFiller, part.count);
        }
    } else {
        if (m == 1) {
            const auto it = std::find(vNum.cbegin(), vNum.cend(),
                                      targetVals.front());
            if (it != vNum.cend()) {
                return Rf_ScalarInteger(1);
            } else {
                return Rf_ScalarInteger(0);
            }
        } else {
            return Rf_ScalarInteger(0);
        }
    }
}
