#include "ConstraintsUtils.h"
#include <cmath>

template <typename typeVector>
inline void PopulateVec(int m, const std::vector<typeVector> &v,
                        std::vector<int> &z, int &count, int maxRows,
                        bool IsComb, std::vector<typeVector> &combinatoricsVec) {
    
    if (IsComb) {
        for (int k = 0; k < m; ++k)
            combinatoricsVec.push_back(v[z[k]]);
        
        ++count;
    } else {
        do {
            for (int k = 0; k < m; ++k)
                combinatoricsVec.push_back(v[z[k]]);
            
            ++count;
        } while (std::next_permutation(z.begin(), z.end()) && count < maxRows);
    }
}

template <typename typeVector>
void SectionOne(const std::vector<typeVector> &v, std::vector<typeVector> &testVec,
                std::vector<int> &z, const std::vector<typeVector> &targetVals,
                std::vector<typeVector> &combinatoricsVec, std::vector<typeVector> &resultsVec,
                bool &t_0, bool &t_1, int &count, partialPtr<typeVector> partialFun,
                funcPtr<typeVector> constraintFun, compPtr<typeVector> compFunOne,
                compPtr<typeVector> compFunTwo, int m, int m1, int maxRows,
                int maxZ, bool IsComb, bool xtraCol) {
    
    for (int i = 0; i < m; ++i)
        testVec[i] = v[z[i]];
    
    const typeVector partialVal = constraintFun(testVec, m1);
    typeVector testVal = partialFun(partialVal, testVec.back(), m);
    t_0 = compFunTwo(testVal, targetVals);
    
    while (t_0 && t_1) {
        if (compFunOne(testVal, targetVals)) {
            int myStart = count;
            PopulateVec(m, v, z, count, maxRows, IsComb, combinatoricsVec);
            
            if (xtraCol)
                for (int i = myStart; i < count; ++i)
                    resultsVec.push_back(testVal);
            
            t_1 = count < maxRows;
        }
        
        t_0 = z[m1] != maxZ;
        
        if (t_0) {
            ++z[m1];
            testVec[m1] = v[z[m1]];
            testVal = partialFun(partialVal, testVec.back(), m);
            t_0 = compFunTwo(testVal, targetVals);
        }
    }
}

distinctType DistinctAttr(int lenV, int m, bool IsRep, bool IsMult, int64_t target,
                          const std::vector<int> &Reps, bool IncludeZero, bool mIsNull) {
    int limit = 0;
    bool getAll = false;
    
    if (IsMult || !IsRep) {
        // The eqn below can be derived by taking note that the
        // smallest number of elements whose sum is at least 
        // the target will be comprised of the first x numbers.
        // That is, we need to solve for x such that:
        //
        //        sum(1:(x - 1)) <= target <= sum(1:x)
        //
        // These are triangle numbers which have the form:
        //
        //              sum(1:x) = x * (x + 1) / 2
        //
        // Given n = target, we have:
        //
        //    x * (x + 1) / 2 >= n  -->>  x^2 + x - 2n >= 0
        //
        // Finally, using the quadratic formula, we obtain:
        //
        //      x = (-1 + sqrt(1 + 4 * 1 * 2n)) / 2 * 1
        //
        // After solving for x, if sum(1:x) > target, we know
        // that the solution with the fewest number of elements
        // will contain x - 1 elements, hence std::floor.
        
        limit = static_cast<int>(std::floor((-1 + std::sqrt(1.0 + 8.0 * target)) / 2));
        
        if (IsMult) {
            // Ensure all elements except the first element are 1. The first
            // element should be zero and thus have a higher frequency in
            // order to test for partitions of different length.
            
            const bool allOne = std::all_of(Reps.cbegin() + 1, Reps.cend(), 
                                            [](int v_i) {return v_i == 1;});
            
            if (IncludeZero && lenV == target + 1 && allOne) {
                if (m >= limit) {
                    getAll = (Reps.front() >= (limit - 1)) ? true : false;
                } else {
                    limit = m;
                }
            } else {
                // N.B. In the calling function we have ensured that if the 
                // freqs arg is invoked with all ones, we set IsMult to false.
                limit = 0;
            }
        } else if (!IsRep) {
            if (m < limit) {
                limit = m;
            } else if (!mIsNull) {
                limit = 0;
            }
        }
    }
    
    // N.B. if limit = 0, this means we either have IsRep = true,
    // or we are not going to use the optimized algorithm. In this
    // case, we revert to the general algorithm.
    
    distinctType res;
    res.limit = limit;
    res.getAll = getAll;
    
    return res;
}

// Check if our function operating on the rows of our matrix can possibly produce elements
// greater than std::numeric_limits<int>::max(). We need a NumericMatrix in this case. We also need to check
// if our function is the mean as this can produce non integral values.
bool CheckIsInteger(const std::string &funPass, int n, int m,
                    const std::vector<double> &vNum, const std::vector<double> &targetVals,
                    funcPtr<double> myFunDbl, bool checkLim) {
    
    if (funPass == "mean")
        return false;
    
    std::vector<double> vAbs;
    
    for (int i = 0; i < n; ++i)
        vAbs.push_back(std::abs(vNum[i]));
    
    const double vecMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
    const std::vector<double> rowVec(m, vecMax);
    const double testIfInt = myFunDbl(rowVec, static_cast<std::size_t>(m));
    
    if (testIfInt > std::numeric_limits<int>::max())
        return false;
    
    if (checkLim) {
        vAbs.clear();
        
        for (std::size_t i = 0; i < targetVals.size(); ++i) {
            if (static_cast<int64_t>(targetVals[i]) != targetVals[i])
                return false;
            else
                vAbs.push_back(std::abs(targetVals[i]));
        }
        
        const double limMax = *std::max_element(vAbs.cbegin(), vAbs.cend());
        
        if (limMax > std::numeric_limits<int>::max())
            return false;
    }
    
    return true;
}

void SetStartPartitionZ(PartitionType PartType, distinctType distinctTest,
                        std::vector<int> &z, const std::vector<int> &Reps,
                        int target, int lenV, int m, bool IncludeZero) {
    
    const bool IsDistinct = (PartType > PartitionType::PartTradNoZero);
    const int lastCol =  (IsDistinct) ? distinctTest.limit - 1 : m - 1;
    const std::size_t partWidth = lastCol + 1;
    z.assign(lastCol + 1 ,0);
    
    switch (PartType) {
        case PartitionType::PartTraditional: {
            z[lastCol] = target;
            break;
        }
        case PartitionType::PartTradNoZero: {
            std::fill(z.begin(), z.end(), 1);
            z[lastCol] = target - partWidth + 1;
            break;
        }
        case PartitionType::PartDstctStdAll: {
            z[lastCol] = target;
            break;
        }
        case PartitionType::PartDstctShort: {
            z[lastCol] = target;
            break;
        }
        case PartitionType::PartDstctSpecial: {
            std::iota(z.begin() + Reps[0] - 1, z.end(), 0);
            z[lastCol] = target - (partWidth - Reps[0]) * (partWidth - Reps[0] - 1) / 2;
            break;
        }
        case PartitionType::PartDstctOneZero: {
            std::iota(z.begin(), z.end(), 0);
            z[lastCol] = target - (partWidth - 1) * (partWidth - 2) / 2;
            break;
        }
        case PartitionType::PartDstctNoZero: {
            std::iota(z.begin(), z.end(), 1);
            z[lastCol] = target - partWidth * (partWidth - 1) / 2;
            break;
        }
        default: {
            z[lastCol] = target;
        }
    }
}

void ConstraintSetup(std::vector<std::string> &compFunVec,
                     std::vector<double> &targetVals, bool &IsBetweenComp) {
    
    if (targetVals.size() > 2) {
        Rcpp::stop("there cannot be more than 2 limitConstraints");
    } else if (targetVals.size() == 2 && targetVals[0] == targetVals[1]) {
        Rcpp::stop("The limitConstraints must be different");
    }
    
    std::vector<std::string>::const_iterator itComp;
    
    if (compFunVec.size() > 2)
        Rcpp::stop("there cannot be more than 2 comparison operators");
    
    for (std::size_t i = 0; i < compFunVec.size(); ++i) {
        itComp = std::find(compForms.cbegin(), compForms.cend(), compFunVec[i]);
        
        if (itComp == compForms.end()) {
            Rcpp::stop("comparison operators must be one of the following: "
                           "'>', '>=', '<', '<=', or '=='");
        }
        
        int myIndex = std::distance(compForms.cbegin(), itComp);
        
        // The first 5 are "standard" whereas the 6th and 7th
        // are written with the equality first. Converting
        // them here makes it easier to deal with later.
        if (myIndex > 4) {myIndex -= 3;}
        compFunVec[i] = compForms[myIndex];
    }
    
    if (compFunVec.size() == 2) {
        if (targetVals.size() == 1) {
            compFunVec.pop_back();
        } else {
            if (compFunVec[0] == "==" || compFunVec[1] == "==") {
                Rcpp::stop("If comparing against two limitConstraints, the "
                               "equality comparisonFun (i.e. '==') cannot be used. "
                               "Instead, use '>=' or '<='.");
            }
            
            if (compFunVec[0].substr(0, 1) == compFunVec[1].substr(0, 1)) {
                Rcpp::stop("Cannot have two 'less than' comparisonFuns or two 'greater than' "
                               "comparisonFuns (E.g. c('<', '<=') is not allowed).");
            }
            
            // The two cases below are for when we are looking for all combs/perms such that when
            // myFun is applied, the result is between 2 values. These comparisons are defined in
            // compSpecial in ConstraintsMaster.h. If we look at the definitions of these comparison
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
                
                if (std::max(targetVals[0], targetVals[1]) == targetVals[1])
                    std::swap(targetVals[0], targetVals[1]);
            }
        }
    } else {
        if (targetVals.size() == 2)
            targetVals.pop_back();
    }
}

void AdjustTargetVals(int n, VecType myType, std::vector<double> &targetVals,
                      std::vector<int> &targetIntVals, const SEXP &Rtolerance,
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
        
        if (Rf_isNull(Rtolerance)) {
            bool IsWhole = true;
            
            for (int i = 0; i < n && IsWhole; ++i)
                if (static_cast<int64_t>(vNum[i]) != vNum[i])
                    IsWhole = false;
                
            for (std::size_t i = 0; i < targetVals.size() && IsWhole; ++i)
                if (static_cast<int64_t>(targetVals[i]) != targetVals[i])
                    IsWhole = false;
                    
            tolerance = (IsWhole && mainFun != "mean") ? 0 : defaultTolerance;
        } else {
            CleanConvert::convertPrimitive(Rtolerance, tolerance, 
                                           "tolerance", true, false, false, true);
        }
        
        
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

template <typename typeVector>
void GetPartitionCase(const std::vector<std::string> &compFunVec, std::vector<typeVector> &v,
                      const std::string &mainFun, const std::vector<typeVector> &target,
                      PartitionType &PartType, distinctType &distinctTest, const SEXP &Rlow,
                      std::vector<int> &Reps, int lenV, int &m, double tolerance, bool IsMult,
                      bool IsRep, bool IsBet, bool mIsNull) {
    
    bool PartitionCase = false;
    PartType = PartitionType::NotPartition;
    bool bLower = false;
    
    if (!Rf_isNull(Rlow)) {
        auto tempLower = FromCpp14::make_unique<mpz_t[]>(1);
        mpz_init(tempLower[0]);
        
        createMPZArray(Rlow, tempLower.get(), 1, "lower");
        bLower = mpz_cmp_si(tempLower[0], 1) > 0;
    }
    
    if (!compFunVec.empty() && !bLower) {
        if (IsMult) {
            for (int i = 0; i < (lenV - 1); ++i) {
                for (int j = (i + 1); j < lenV; ++j) {
                    if (v[i] > v[j]) {
                        std::swap(v[i], v[j]);
                        std::swap(Reps[i], Reps[j]);
                    }
                }
            }
        } else {
            std::sort(v.begin(), v.end());
        }
        
        int64_t tarTest = 0;
        
        // We don't need to check the edge case when lenV == target and there is a
        // zero. Remember, lenV is the length of the vector v, so we could have a
        // situation where v = c(0, 2, 3, 4, 5) -->> length(v) = 5. This would
        // cause a problem if we were to allow this through, however, in the
        // calling code (i.e. Combinatorics.cpp), we ensure that the distance
        // between each element is the same. This means for the example we gave,
        // we would have length(unique(diff(v))) > 1, which means PartitionCase
        // (see Combinatorics.cpp) would be false, and thus the general
        // algorithm would be executed.
        //
        // We do have to ensure that the smallest element is non-negative, othe-
        // rwise, cases like v = seq(-8, 10, 2), m = 7, rep = TRUE, & limit = 10
        // would pass as v = 0:9, m = 7, rep = TRUE, & limit = 9, --or--
        // v = 1:10, m = 7, rep = TRUE, & limit = 10
        
        if (compFunVec[0] == "==" && mainFun == "sum" && lenV > 1 && m > 1) {
            std::vector<typeVector> pTest(v.cbegin(), v.cend());
            std::sort(pTest.begin(), pTest.end());
            const typeVector tarDiff = pTest[1] - pTest[0];
            
            if (static_cast<int64_t>(pTest[0]) == pTest[0]) {
                PartitionCase = true;
                
                for (int i = 1; i < lenV; ++i) {
                    const typeVector testDiff = pTest[i] - pTest[i - 1];
                    
                    if (std::abs(testDiff - tarDiff) > tolerance
                            || static_cast<int64_t>(pTest[i]) != pTest[i]) {
                        PartitionCase = false;
                        break;
                    }
                }
            }
            
            if (target.size() == 1 || target.front() == target.back()) {
                tarTest = static_cast<int64_t>(target.front());
                
                if (PartitionCase)
                    PartitionCase = (tarTest == target.front());
            } else {
                PartitionCase = false;
            }
        }
        
        if (PartitionCase) {
            const typeVector myMax = v.back();
            const bool IncludeZero = (v.front() == 0);
            PartType = PartitionType::PartGeneral;
            
            if (myMax == tarTest && (lenV == tarTest || lenV == (tarTest + 1)) && v.front() >= 0) {
                distinctTest = DistinctAttr(lenV, m, IsRep, IsMult,
                                            tarTest, Reps, IncludeZero, mIsNull);
                
                if (distinctTest.limit > 0) {
                    m = distinctTest.limit;
                    
                    if (IncludeZero) {
                        if (IsMult) {
                            if (distinctTest.getAll) {
                                PartType = PartitionType::PartDstctStdAll;
                            } else if (Reps[0] >= (m - 1)) {
                                PartType = PartitionType::PartDstctShort;
                            } else {
                                PartType = PartitionType::PartDstctSpecial;
                            }
                        } else {
                            PartType = PartitionType::PartDstctOneZero;
                        }
                    } else {
                        PartType = PartitionType::PartDstctNoZero;
                    }
                } else if (IsRep) {
                    if (IncludeZero) {
                        PartType = PartitionType::PartTraditional;
                    } else {
                        PartType = PartitionType::PartTradNoZero;
                    }
                    
                    if (m >= lenV) { 
                        if (IncludeZero) {
                            m = lenV - 1;
                        } else {
                            PartType = PartitionType::NotPartition;
                        }
                    }
                }
            }
        } else if ((compFunVec[0] == "==" || IsBet) && lenV > 1
                && m > 1 && mainFun != "max" && mainFun != "min") {
            PartType = PartitionType::PartitonEsque;
        }
    }
}

template <typename typeVector>
bool CheckSpecialCase(int n, bool bLower, const std::string &mainFun,
                      const std::vector<typeVector> &vNum) {
    
    bool result = false;
    
    // If bLower, the user is looking to test a particular range. Otherwise, the constraint algo
    // will simply return (upper - lower) # of combinations/permutations that meet the criteria
    if (bLower) {
        result = true;
    } else if (mainFun == "prod") {
        for (int i = 0; i < n; ++i) {
            if (vNum[i] < 0) {
                result = true;
                break;
            }
        }
    }
    
    return result;
}

template void SectionOne(const std::vector<int>&, std::vector<int> &testVec,
                         std::vector<int>&, const std::vector<int>&,
                         std::vector<int>&, std::vector<int>&,
                         bool&, bool&, int&, partialPtr<int>, funcPtr<int>,
                         compPtr<int>, compPtr<int>, int, int, int, int, bool, bool);
template void SectionOne(const std::vector<double>&, std::vector<double> &testVec,
                         std::vector<int>&, const std::vector<double>&,
                         std::vector<double>&, std::vector<double>&,
                         bool&, bool&, int&, partialPtr<double>, funcPtr<double>,
                         compPtr<double>, compPtr<double>, int, int, int, int, bool, bool);

template void GetPartitionCase(const std::vector<std::string>&, std::vector<int>&,
                               const std::string&, const std::vector<int>&, PartitionType&,
                               distinctType&, const SEXP&, std::vector<int>&, int, int&,
                               double, bool, bool, bool, bool);
template void GetPartitionCase(const std::vector<std::string>&, std::vector<double>&,
                               const std::string&, const std::vector<double>&, PartitionType&,
                               distinctType&, const SEXP&, std::vector<int>&, int, int&,
                               double, bool, bool, bool, bool);

template bool CheckSpecialCase(int, bool, const std::string&, const std::vector<int>&);
template bool CheckSpecialCase(int, bool, const std::string&, const std::vector<double>&);
