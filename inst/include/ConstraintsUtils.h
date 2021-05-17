#ifndef CONSTRAINTS_UTILS_H
#define CONSTRAINTS_UTILS_H

#include "UserConstraintFuns.h"
#include "GmpDependUtils.h"
#include <chrono>

const std::vector<std::string> compForms = {"<", ">", "<=", ">=", "==", "=<", "=>"};
const std::vector<std::string> compSpecial = {"==", ">,<", ">=,<", ">,<=", ">=,<="};
const std::vector<std::string> compHelper = {"<=", "<", "<", "<=", "<="};

// in R console: print(sqrt(.Machine$double.eps), digits = 16)
// [1] 0.00000001490116119384766
// Which is also 2^(-26)
constexpr double defaultTolerance = 0.00000001490116119384766;

// Used for checking whether user has interrupted computation
constexpr auto timeout = std::chrono::milliseconds(1000);

struct distinctType {
    int limit = 0;
    bool getAll = false;
};

distinctType DistinctAttr(int lenV, int m, bool IsRep, bool IsMult, int64_t target,
                          const std::vector<int> &Reps, bool IncludeZero);

bool CheckIsInteger(const std::string &funPass, int n, int m,
                    const std::vector<double> &vNum, const std::vector<double> &targetVals,
                    funcPtr<double> myFunDbl, bool checkLim = false);

void SetStartPartitionZ(PartitionType PartType, distinctType distinctTest,
                        std::vector<int> &z, const std::vector<int> &Reps,
                        int target, int lenV, int m, bool IncludeZero);

void ConstraintSetup(std::vector<std::string> &compFunVec,
                     std::vector<double> &targetVals, bool &IsBetweenComp);

void AdjustTargetVals(int n, VecType myType, std::vector<double> &targetVals,
                      std::vector<int> &targetIntVals, const SEXP &Rtolerance,
                      std::vector<std::string> &compFunVec, double &tolerance,
                      const std::string &mainFun, const std::vector<double> &vNum);

template <typename typeVector>
using partialReducePtr = void (*const)(int m, typeVector &partial, typeVector w);

template <typename typeVector>
void PartialReduceProd(int m, typeVector &partial, typeVector w) {
    partial /= w;
}

template <typename typeVector>
void PartialReduceSum(int m, typeVector &partial, typeVector w) {
    partial -= w;
}

template <typename typeVector>
void PartialReduceMean(int m, typeVector& partial, typeVector w) {
    partial = (partial * static_cast<double>(m) - w) / static_cast<double>(m - 1);
}

template <typename T>
Rcpp::XPtr<partialReducePtr<T>> putPartialReduceInXPtr(const std::string &myFun) {
    
    if (myFun == "prod") {
        return(Rcpp::XPtr<partialReducePtr<T>>(new partialReducePtr<T>(&PartialReduceProd)));
    } else if (myFun == "sum") {
        return(Rcpp::XPtr<partialReducePtr<T>>(new partialReducePtr<T>(&PartialReduceSum)));
    } else {
        return(Rcpp::XPtr<partialReducePtr<T>>(new partialReducePtr<T>(&PartialReduceMean)));
    }
}

template <typename typeVector>
void GetPartitionCase(const std::vector<std::string> &compFunVec, std::vector<typeVector> &v,
                      const std::string &mainFun, const std::vector<typeVector> &target,
                      PartitionType &PartType, distinctType &distinctTest, const SEXP &Rlow,
                      std::vector<int> &Reps, int lenV, int &m, double tolerance, bool IsMult,
                      bool IsRep, bool IsBet, bool mIsNull);

template <typename typeVector>
bool CheckSpecialCase(int n, bool bLower, const std::string &mainFun,
                      const std::vector<typeVector> &vNum);

template <typename typeVector>
void SectionOne(const std::vector<typeVector> &v, std::vector<typeVector> &testVec,
                std::vector<int> &z, const std::vector<typeVector> &targetVals,
                std::vector<typeVector> &combinatoricsVec, std::vector<typeVector> &resultsVec,
                bool &t_0, bool &t_1, int &count, partialPtr<typeVector> partialFun,
                funcPtr<typeVector> constraintFun, compPtr<typeVector> compFunOne,
                compPtr<typeVector> compFunTwo, int m, int m1, int maxRows,
                int maxZ, bool IsComb, bool xtraCol);

#endif
