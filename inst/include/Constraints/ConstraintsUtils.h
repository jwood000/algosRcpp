#ifndef CONSTRAINTS_UTILS_H
#define CONSTRAINTS_UTILS_H

#include "Constraints/UserConstraintFuns.h"
#include <string>
#include <limits>
#include <vector>
#include <cstdint>

bool CheckIsInteger(const std::string &funPass, int n,
                    int m, const std::vector<double> &vNum,
                    const std::vector<double> &targetVals,
                    const funcPtr<double> myFunDbl,
                    bool checkLim = false);

#endif