#ifndef CONSTRAINTS_SPECIAL_H
#define CONSTRAINTS_SPECIAL_H

#include <vector>
#include <string>

template <typename typeRcpp, typename T>
typeRcpp ConstraintsSpecial(int n, int m, std::vector<T> v, bool IsRep, int nRows, 
                            bool keepRes, std::vector<int> z, double lower, std::string mainFun, 
                            bool IsMult, double computedRows, std::vector<std::string> compFunVec,
                            std::vector<T> targetVals, bool IsComb, std::vector<int> myReps,
                            std::vector<int> freqs, bool bLower, double userRows);

#endif
