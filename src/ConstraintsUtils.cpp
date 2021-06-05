#include "Constraints/ConstraintsUtils.h"
#include "CleanConvert.h"

bool CheckSpecialCase(bool bLower, const std::string &mainFun,
                      const std::vector<double> &vNum) {
    
    bool result = false;
    
    // If bLower, the user is looking to test a particular range. Otherwise, the constraint algo
    // will simply return (upper - lower) # of combinations/permutations that meet the criteria
    if (bLower) {
        result = true;
    } else if (mainFun == "prod") {
        for (auto v_i : vNum) {
            if (v_i < 0) {
                result = true;
                break;
            }
        }
    }
    
    return result;
}

void SetConstraintType(const std::vector<double> &vNum,
                       const PartDesign &part, const std::string &mainFun,
                       ConstraintType &ctype, bool IsConstrained, bool bLower) {

    if (IsConstrained) {
        if (CheckSpecialCase(bLower, mainFun, vNum)) {
            ctype = ConstraintType::SpecialCnstrnt;
        } else if (part.ptype == PartitionType::CoarseGrained) {
            ctype = ConstraintType::PartitionEsque;
        } else if (ctype < ConstraintType::PartMapping) {
            ctype = ConstraintType::General;
        } else {
            // If we get here, ctype has already been set to
            // PartMapping or PartStandard in SetPartitionDesign
        }
    } else {
        ctype = ConstraintType::NoConstraint;
    }
}

void ConstraintStructure(std::vector<std::string> &compFunVec,
                         std::vector<double> &targetVals,
                         bool &IsBetweenComp) {

    if (targetVals.size() > 2) {
        Rf_error("there cannot be more than 2 limitConstraints");
    } else if (targetVals.size() == 2 && targetVals[0] == targetVals[1]) {
        Rf_error("The limitConstraints must be different");
    }

    if (compFunVec.size() > 2)
        Rf_error("there cannot be more than 2 comparison operators");

    // The first 5 are "standard" whereas the 6th and 7th
    // are written with the equality first. Converting
    // them here makes it easier to deal with later.
    for (std::size_t i = 0; i < compFunVec.size(); ++i) {
        auto itComp = compForms.find(compFunVec[i]);

        if (itComp == compForms.end()) {
            Rf_error("comparison operators must be one of the"
                     " following: '>", ">=", "<", "<=', or '=='");
        } else {
            compFunVec[i] = itComp->second;
        }
    }

    if (compFunVec.size() == 2) {
        if (targetVals.size() == 1) {
            compFunVec.pop_back();
        } else {
            if (compFunVec[0] == "==" || compFunVec[1] == "==") {
                Rf_error("If comparing against two limitConstraints, the "
                         "equality comparisonFun (i.e. '==') cannot be used."
                         " Instead, use '>=' or '<='.");
            }

            if (compFunVec[0].substr(0, 1) == compFunVec[1].substr(0, 1)) {
                Rf_error("Cannot have two 'less than' comparisonFuns or two"
                        " 'greater than' comparisonFuns (E.g. c('<", "<=')"
                        " is not allowed).");
            }

            // The two cases below are for when we are looking for all combs/perms such that when
            // myFun is applied, the result is between 2 values. These comparisons are defined in
            // compSpecial in ConstraintsMain.h. If we look at the definitions of these comparison
            // functions in ConstraintsUtils.h, we see the following trend for all 4 cases:
            //
            // bool greaterLess(stdType x, const std::vector<stdType> &y) {return (x < y[0]) && (x > y[1]);}
            //
            // That is, we need the maixmum value in targetVals to be the first value and the second
            // value needs to be the minimum. At this point, the constraint algorithm will be
            // identical to when comparisonFun = "==" (i.e. we allow the algorithm to continue
            // while the results are less than (or equal to in cases where strict inequalities are
            // enforced) the target value and stop once the result exceeds that value).

            if (compFunVec[0].substr(0, 1) == ">" &&
                std::min(targetVals[0], targetVals[1]) == targetVals[0]) {

                compFunVec[0] = compFunVec[0] + "," + compFunVec[1];
                IsBetweenComp = true;
            } else if (compFunVec[0].substr(0, 1) == "<" &&
                std::max(targetVals[0], targetVals[1]) == targetVals[0]) {

                compFunVec[0] = compFunVec[1] + "," + compFunVec[0];
                IsBetweenComp = true;
            }

            if (IsBetweenComp) {
                compFunVec.pop_back();

                if (std::max(targetVals[0], targetVals[1]) == targetVals[1]) {
                    std::swap(targetVals[0], targetVals[1]);
                }
            }
        }
    } else {
        if (targetVals.size() == 2) {
            targetVals.pop_back();
        }
    }
}

void SetTolerance(const std::vector<double> &vNum,
                  const std::vector<double> &targetVals,
                  const std::string &mainFun,
                  SEXP Rtolerance, double &tolerance) {
    
    if (Rf_isNull(Rtolerance)) {
        bool IsWhole = true;
        
        for (std::size_t i = 0; i < vNum.size() && IsWhole; ++i) {
            if (static_cast<std::int64_t>(vNum[i]) != vNum[i]) {
                IsWhole = false;
            }
        }
        
        for (std::size_t i = 0; i < targetVals.size() && IsWhole; ++i) {
            if (static_cast<std::int64_t>(targetVals[i]) != targetVals[i]) {
                IsWhole = false;
            }
        }
        
        tolerance = (IsWhole && mainFun != "mean") ? 0 : defaultTolerance;
    } else {
        // numOnly = true, checkWhole = false, negPoss = false, decimalFraction = true
        CleanConvert::convertPrimitive(Rtolerance, tolerance, VecType::Numeric,
                                       "tolerance", true, false, false, true);
    }
}

void AdjustTargetVals(VecType myType, std::vector<double> &targetVals,
                      std::vector<int> &targetIntVals, SEXP Rtolerance,
                      std::vector<std::string> &compFunVec, double &tolerance,
                      const std::string &mainFun, const std::vector<double> &vNum) {

    if (myType == VecType::Integer) {
        targetIntVals.assign(targetVals.cbegin(), targetVals.cend());
    } else {
        // We first check if we are getting double precision. If so, for the
        // non-strict inequalities, we must alter the limit by some epsilon:
        //
        //                 x <= y   --->>>   x <= y + e
        //                 x >= y   --->>>   x >= y - e
        //
        // Equality is a bit tricky as we need to check whether the absolute
        // value of the difference is less than epsilon. As a result, we
        // can't alter the limit with one alteration. Observe:
        //
        //          x == y  --->>>  |x - y| <= e , which gives:
        //
        //                      -e <= x - y <= e
        //
        //                    1.     x >= y - e
        //                    2.     x <= y + e
        //
        // As a result, we must define a specialized equality check for double
        // precision. It is 'equalDbl' and can be found in ConstraintsUtils.h
        
        SetTolerance(vNum, targetVals, mainFun, Rtolerance, tolerance);
        const auto itComp = std::find(compSpecial.cbegin(),
                                      compSpecial.cend(), compFunVec[0]);

        if (compFunVec[0] == "==") {
            targetVals.push_back(targetVals[0] - tolerance);
            targetVals[0] += tolerance;
        } else if (itComp != compSpecial.end()) {
            targetVals[0] += tolerance;
            targetVals[1] -= tolerance;
        } else if (compFunVec[0] == "<=") {
            targetVals[0] += tolerance;
        } else if (compFunVec[0] == ">=") {
            targetVals[0] -= tolerance;
        }

        if (compFunVec.size() > 1) {
            if (compFunVec[1] == "<=") {
                targetVals[1] += tolerance;
            } else if (compFunVec[1] == ">=") {
                targetVals[1] -= tolerance;
            }
        }
    }
}

// Check if our function operating on the rows of our matrix can possibly produce elements
// greater than std::numeric_limits<int>::max(). We need a NumericMatrix in this case. We also need to check
// if our function is the mean as this can produce non integral values.
bool CheckIsInteger(const std::string &funPass, int n,
                    int m, const std::vector<double> &vNum,
                    const std::vector<double> &targetVals,
                    const funcPtr<double> myFunDbl, bool checkLim) {

    if (funPass == "mean") {
        return false;
    }

    std::vector<double> vAbs;

    for (int i = 0; i < n; ++i) {
        vAbs.push_back(std::abs(vNum[i]));
    }

    const double vecMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
    const std::vector<double> rowVec(m, vecMax);
    const double testIfInt = myFunDbl(rowVec, static_cast<std::size_t>(m));

    if (testIfInt > std::numeric_limits<int>::max()) {
        return false;
    }

    if (checkLim) {
        vAbs.clear();

        for (std::size_t i = 0; i < targetVals.size(); ++i) {
            if (static_cast<std::int64_t>(targetVals[i]) != targetVals[i]) {
                return false;
            } else {
                vAbs.push_back(std::abs(targetVals[i]));
            }
        }

        const double limMax = *std::max_element(vAbs.cbegin(), vAbs.cend());

        if (limMax > std::numeric_limits<int>::max()) {
            return false;
        }
    }

    return true;
}

void ConstraintSetup(const std::vector<double> &vNum,
                     const std::vector<int> &Reps,
                     std::vector<double> &targetVals,
                     std::vector<int> &targetIntVals,
                     const funcPtr<double> funDbl, PartDesign &part,
                     ConstraintType &ctype, int lenV, int m,
                     std::vector<std::string> &compFunVec,
                     const std::string &mainFun, VecType &myType,
                     SEXP Rtarget, SEXP RcompFun, SEXP Rtolerance,
                     SEXP Rlow, bool IsConstrained, bool bCalcMulti) {

    if (IsConstrained) {
        // numOnly = true, checkWhole = false, negPoss = true
        CleanConvert::convertVector(Rtarget, targetVals,
                                    VecType::Numeric,
                                    "limitConstraints",
                                    true, false, true);

        int len_comp = Rf_length(RcompFun);

        for (int i = 0; i < len_comp; ++i) {
            const std::string temp(CHAR(STRING_ELT(RcompFun, i)));
            compFunVec.push_back(temp);
        }

        bool IsBetweenComp = false;
        ConstraintStructure(compFunVec, targetVals, IsBetweenComp);

        if (myType == VecType::Integer &&
            !CheckIsInteger(mainFun, lenV, m, vNum, targetVals, funDbl, true)) {

            myType = VecType::Numeric;
        }

        double tolerance = 0;
        AdjustTargetVals(myType, targetVals, targetIntVals,
                         Rtolerance, compFunVec, tolerance, mainFun, vNum);
        
        const bool IsPartition = CheckPartition(compFunVec, vNum, mainFun,
                                                targetVals, part, Rlow, lenV,
                                                m, tolerance, IsBetweenComp);

        if (IsPartition) {
            SetPartitionDesign(Reps, vNum, part, ctype, lenV, m, bCalcMulti);
        }
    }
}
