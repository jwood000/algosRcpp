#include <Combinations.h>
#include <Permutations.h>
#include <ConstraintsUtils.h>
#include <NthResult.h>
#include <CountGmp.h>
#include <CleanConvert.h>

const std::vector<std::string> compForms = {"<", ">", "<=", ">=", "==", "=<", "=>"};
const std::vector<std::string> compSpecial = {"==", ">,<", ">=,<", ">,<=", ">=,<="};
const std::vector<std::string> compHelper = {"<=", "<", "<", "<=", "<="};

const double Significand53 = 9007199254740991.0;

template <typename typeRcpp>
typeRcpp SubMat(typeRcpp m, int n) {
    int k = m.ncol();
    typeRcpp subMatrix = Rcpp::no_init_matrix(n, k);
    
    for (int i = 0; i < n; ++i)
        subMatrix(i, Rcpp::_) = m(i, Rcpp::_);
    
    return subMatrix;
}

// This function applys a constraint function to a vector v with respect
// to a constraint value "lim". The main idea is that combinations are added
// successively, until a particular combination exceeds the given constraint
// value for a given constraint function. After this point, we can safely skip
// several combinations knowing that they will exceed the given constraint value.

template <typename typeRcpp, typename typeVector>
typeRcpp CombinatoricsConstraints(int n, int r, std::vector<typeVector> &v, bool repetition,
                                  std::string myFun, std::vector<std::string> comparison, 
                                  std::vector<typeVector> lim, int numRows, bool isComb,
                                  bool xtraCol, std::vector<int> &Reps, bool isMult) {
    
    // myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max";
    // and myComparison is a comparison operator: "<", "<=", ">", ">=", or "==";
    
    typeVector testVal;
    int i, j, count = 0, vSize = v.size(), numCols;
    numCols = xtraCol ? (r + 1) : r;
    unsigned long int uR = r;
    typeRcpp combinatoricsMatrix = Rcpp::no_init_matrix(numRows, numCols);
    
    Rcpp::XPtr<funcPtr<typeVector> > xpFun = putFunPtrInXPtr<typeVector>(myFun);
    funcPtr<typeVector> constraintFun = *xpFun;
    
    // We first check if we are getting double precision.
    // If so, for the non-stric inequalities, we have
    // to alter the limit by epsilon:
    //
    //           x <= y   --->>>   x <= y + e
    //           x >= y   --->>>   x >= y - e
    //
    // Equality is a bit tricky as we need to check
    // whether the absolute value of the difference is
    // less than epsilon. As a result, we can't alter
    // the limit with one alteration. Observe:
    //
    //   x == y  --->>>  |x - y| <= e , which gives:
    //
    //             - e <= x - y <= e
    //
    //         1.     x >= y - e
    //         2.     x <= y + e
    //
    // As a result, we must define a specialized equality
    // check for double precision. It is 'equalDbl' and
    // can be found in ConstraintsUtils.h
    
    if (!std::is_integral<typeVector>::value) {
        if (comparison[0] == "<=") {
            lim[0] = lim[0] + std::numeric_limits<double>::epsilon();
        } else if (comparison[0] == ">=") {
            lim[0] = lim[0] - std::numeric_limits<double>::epsilon();
        }
        
        if (comparison.size() > 1) {
            if (comparison[1] == "<=") {
                lim[1] = lim[1] + std::numeric_limits<double>::epsilon();
            } else if (comparison[1] == ">=") {
                lim[1] = lim[1] - std::numeric_limits<double>::epsilon();
            }
        }
    }
    
    for (std::size_t nC = 0; nC < comparison.size(); ++nC) {
        
        Rcpp::XPtr<compPtr<typeVector>> xpCompOne = putCompPtrInXPtr<typeVector>(comparison[nC]);
        compPtr<typeVector> comparisonFunOne = *xpCompOne;
        
        Rcpp::XPtr<compPtr<typeVector>> xpCompTwo = xpCompOne;
        compPtr<typeVector> comparisonFunTwo;
    
        if (comparison[nC] == ">" || comparison[nC] == ">=") {
            if (isMult) {
                for (i = 0; i < (n - 1); ++i) {
                    for (j = (i + 1); j < n; ++j) {
                        if (v[i] < v[j]) {
                            std::swap(v[i], v[j]);
                            std::swap(Reps[i], Reps[j]);
                        }
                    }
                }
            } else {
                std::sort(v.begin(), v.end(), std::greater<double>());
            }
            comparisonFunTwo = *xpCompOne;
        } else {
            if (isMult) {
                for (i = 0; i < (n-1); ++i) {
                    for (j = (i+1); j < n; ++j) {
                        if (v[i] > v[j]) {
                            std::swap(v[i], v[j]);
                            std::swap(Reps[i], Reps[j]);
                        }
                    }
                }
            } else {
                std::sort(v.begin(), v.end());
            }
            
            std::vector<std::string>::const_iterator itComp = std::find(compSpecial.begin(), 
                                                                        compSpecial.end(), 
                                                                        comparison[nC]);
            if (itComp != compSpecial.end()) {
                int myIndex = std::distance(compSpecial.begin(), itComp);
                Rcpp::XPtr<compPtr<typeVector> > xpCompThree = putCompPtrInXPtr<typeVector>(compHelper[myIndex]);
                comparisonFunTwo = *xpCompThree;
            } else {
                comparisonFunTwo = *xpCompOne;
            }
        }
        
        std::vector<int> z, zCheck, zPerm(r);
        std::vector<typeVector> testVec(r);
        bool t_1, t_2, t = true, keepGoing = true;
        int r1 = r - 1, r2 = r - 2, k = 0; 
        int numIter, maxZ = n - 1;
        
        if (isMult) {
            int zExpSize = 0;
            std::vector<int> zExpand, zIndex, zGroup(r);
            
            for (i = 0; i < n; ++i)
                zExpSize += Reps[i];
            
            for (i = 0; i < n; ++i) {
                zIndex.push_back(k);
                for (j = 0; j < Reps[i]; ++j, ++k)
                    zExpand.push_back(i);
            }
            
            for (i = 0; i < r; ++i)
                z.push_back(zExpand[i]);
            
            while (keepGoing) {
                
                t_2 = true;
                for (i = 0; i < r; ++i)
                    testVec[i] = v[zExpand[zIndex[z[i]]]];
                
                testVal = constraintFun(testVec, uR);
                t = comparisonFunTwo(testVal, lim);
                
                while (t && t_2 && keepGoing) {
                    
                    testVal = constraintFun(testVec, uR);
                    t_1 = comparisonFunOne(testVal, lim);
                    
                    if (t_1) {
                        if (isComb) {
                            for (k = 0; k < r; ++k)
                                combinatoricsMatrix(count, k) = v[zExpand[zIndex[z[k]]]];
                            
                            ++count;
                        } else {
                            for (k = 0; k < r; ++k)
                                zPerm[k] = zExpand[zIndex[z[k]]];
                            
                            numIter = (int) NumPermsWithRep(zPerm);
                            if ((numIter + count) > numRows)
                                numIter = numRows - count;
                            
                            for (i = 0; i < numIter; ++i, ++count) {
                                for (k = 0; k < r; ++k)
                                    combinatoricsMatrix(count, k) = v[zPerm[k]];
                                
                                std::next_permutation(zPerm.begin(), zPerm.end());
                            }
                        }
                    }
                    
                    keepGoing = (count < numRows);
                    t_2 = (z[r1] != maxZ);
                    
                    if (t_2) {
                        ++z[r1];
                        testVec[r1] = v[zExpand[zIndex[z[r1]]]];
                        testVal = constraintFun(testVec, uR);
                        t = comparisonFunTwo(testVal, lim);
                    }
                }
                
                if (keepGoing) {
                    zCheck = z;
                    for (i = r2; i >= 0; --i) {
                        if (zExpand[zIndex[z[i]]] != zExpand[zExpSize - r + i]) {
                            ++z[i];
                            testVec[i] = v[zExpand[zIndex[z[i]]]];
                            zGroup[i] = zIndex[z[i]];
                            for (k = (i+1); k < r; ++k) {
                                zGroup[k] = zGroup[k-1] + 1;
                                z[k] = zExpand[zGroup[k]];
                                testVec[k] = v[zExpand[zIndex[z[k]]]];
                            }
                            testVal = constraintFun(testVec, uR);
                            t = comparisonFunTwo(testVal, lim);
                            if (t) {break;}
                        }
                    }
                    
                    if (!t) {
                        keepGoing = false;
                    } else if (zCheck == z) {
                        keepGoing = false;
                    }
                }
            }
            
        } else if (repetition) {
            
            v.erase(std::unique(v.begin(), v.end()), v.end());
            vSize = v.size();
            z.assign(r, 0);
            maxZ = vSize - 1;
            
            while (keepGoing) {
                
                t_2 = true;
                for (i = 0; i < r; ++i)
                    testVec[i] = v[z[i]];
                
                testVal = constraintFun(testVec, uR);
                t = comparisonFunTwo(testVal, lim);
                
                while (t && t_2 && keepGoing) {
                    
                    testVal = constraintFun(testVec, uR);
                    t_1 = comparisonFunOne(testVal, lim);
                    
                    if (t_1) {
                        if (isComb) {
                            for (k = 0; k < r; ++k)
                                combinatoricsMatrix(count, k) = v[z[k]];
                            
                            ++count;
                        } else {
                            zPerm = z;
                            numIter = (int) NumPermsWithRep(zPerm);
                            if ((numIter + count) > numRows)
                                numIter = numRows - count;
                            
                            for (i = 0; i < numIter; ++i, ++count) {
                                for (k = 0; k < r; ++k)
                                    combinatoricsMatrix(count, k) = v[zPerm[k]];
                                
                                std::next_permutation(zPerm.begin(), zPerm.end());
                            }
                        }
                        
                        keepGoing = (count < numRows);
                    }
                    
                    t_2 = (z[r1] != maxZ);
                    
                    if (t_2) {
                        ++z[r1];
                        testVec[r1] = v[z[r1]];
                        testVal = constraintFun(testVec, uR);
                        t = comparisonFunTwo(testVal, lim);
                    }
                }
                
                if (keepGoing) {
                    zCheck = z;
                    for (i = r2; i >= 0; --i) {
                        if (z[i] != maxZ) {
                            ++z[i];
                            testVec[i] = v[z[i]];
                            for (k = (i+1); k < r; ++k) {
                                z[k] = z[k-1];
                                testVec[k] = v[z[k]];
                            }
                            testVal = constraintFun(testVec, uR);
                            t = comparisonFunTwo(testVal, lim);
                            if (t) {break;}
                        }
                    }
                    
                    if (!t) {
                        keepGoing = false;
                    } else if (zCheck == z) {
                        keepGoing = false;
                    }
                }
            }
            
        } else {
            
            for (i = 0; i < r; ++i)
                z.push_back(i);
            
            int indexRows, nMinusR = (n - r), myRow;
            indexRows = isComb ? 0 : (int) NumPermsNoRep(r, r1);
            uint8_t *indexMatrix = new uint8_t[indexRows * r];
            
            if (!isComb) {
                indexRows = (int) NumPermsNoRep(r, r1);
                std::vector<uint8_t> indexVec(r);
                std::iota(indexVec.begin(), indexVec.end(), 0);
                
                for (i = 0, myRow = 0; i < indexRows; ++i, myRow += r) {
                    for (j = 0; j < r; ++j)
                        indexMatrix[myRow + j] = indexVec[j];
                    
                    std::next_permutation(indexVec.begin(), indexVec.end());
                }
            }
    
            while (keepGoing) {
    
                t_2 = true;
                for (i = 0; i < r; ++i)
                    testVec[i] = v[z[i]];
    
                testVal = constraintFun(testVec, uR);
                t = comparisonFunTwo(testVal, lim);
    
                while (t && t_2 && keepGoing) {
    
                    testVal = constraintFun(testVec, uR);
                    t_1 = comparisonFunOne(testVal, lim);
    
                    if (t_1) {
                        if (isComb) {
                            for (k=0; k < r; ++k)
                                combinatoricsMatrix(count, k) = v[z[k]];
    
                            ++count;
                        } else {
                            if (indexRows + count > numRows)
                                indexRows = numRows - count;
    
                            for (j = 0, myRow = 0; j < indexRows; ++j, ++count, myRow += r)
                                for (k = 0; k < r; ++k)
                                    combinatoricsMatrix(count, k) = v[z[indexMatrix[myRow + k]]];
                        }
    
                        keepGoing = (count < numRows);
                    }
    
                    t_2 = (z[r1] != maxZ);
    
                    if (t_2) {
                        ++z[r1];
                        testVec[r1] = v[z[r1]];
                        testVal = constraintFun(testVec, uR);
                        t = comparisonFunTwo(testVal, lim);
                    }
                }
    
                if (keepGoing) {
                    zCheck = z;
                    for (i = r2; i >= 0; --i) {
                        if (z[i] != (nMinusR + i)) {
                            ++z[i];
                            testVec[i] = v[z[i]];
                            for (k = (i+1); k < r; ++k) {
                                z[k] = z[k-1]+1;
                                testVec[k] = v[z[k]];
                            }
                            testVal = constraintFun(testVec, uR);
                            t = comparisonFunTwo(testVal, lim);
                            if (t) {break;}
                        }
                    }
                    if (!t) {
                        keepGoing = false;
                    } else if (zCheck == z) {
                        keepGoing = false;
                    }
                }
            }
            
            delete[] indexMatrix;
        }
        
        lim.erase(lim.begin());
    }
       
    return SubMat(combinatoricsMatrix, count);
}

// This is called when we can't easily produce a monotonic sequence overall,
// and we must generate and test every possible combination/permutation
template <typename typeRcpp, typename typeVector>
typeRcpp SpecCaseRet(int n, int m, std::vector<typeVector> v, bool IsRep, int nRows, bool keepRes,
                     std::vector<int> z, std::string mainFun1, bool IsMultiset, double computedRows,
                     std::vector<std::string> compFunVec, std::vector<typeVector> myLim, bool IsComb,
                     std::vector<int> myReps, std::vector<int> freqs, bool bLower, bool permNonTriv, double userRows) {
    
    std::vector<typeVector> rowVec(m);
    
    if (!bLower) {
        if (computedRows > INT_MAX)
            Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
        
        nRows = (int) computedRows;
    }
    
    bool Success;
    std::vector<int> indexMatch;
    indexMatch.reserve(nRows);
    
    Rcpp::XPtr<funcPtr<typeVector> > xpFun1 = putFunPtrInXPtr<typeVector>(mainFun1);
    funcPtr<typeVector> myFun1 = *xpFun1;
    
    typeRcpp matRes;
    typeVector testVal;
    
    std::vector<typeVector> matchVals;
    matchVals.reserve(nRows);
    
    if (IsComb) {
        if (IsMultiset) {
            matRes = Combinations::MultisetCombination<typeRcpp>(n, m, v, myReps, freqs, nRows, keepRes, z);
        } else {
            matRes = Combinations::ComboGeneral<typeRcpp>(n, m, v, IsRep, nRows, keepRes, z);
        }
    } else {
        if (IsMultiset) {
            matRes = Permutations::MultisetPermutation<typeRcpp>(n, m, v, nRows, keepRes, z);
        } else {
            matRes = Permutations::PermuteGeneral<typeRcpp>(n, m, v, IsRep, nRows, keepRes, z, permNonTriv);
        }
    }
    
    Rcpp::XPtr<compPtr<typeVector> > xpComp = putCompPtrInXPtr<typeVector>(compFunVec[0]);
    compPtr<typeVector> myComp = *xpComp;
    
    Rcpp::XPtr<compPtr<typeVector> > xpComp2 = xpComp;
    compPtr<typeVector> myComp2;
    
    unsigned long int uM = m;
    
    if (compFunVec.size() == 1) {
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < m; ++j)
                rowVec[j] = matRes(i, j);
            
            testVal = myFun1(rowVec, uM);
            Success = myComp(testVal, myLim);
            if (Success) {
                indexMatch.push_back(i);
                matchVals.push_back(testVal);
            }
        }
    } else {
        xpComp2 = putCompPtrInXPtr<typeVector>(compFunVec[1]);
        myComp2 = *xpComp2;
        std::vector<typeVector> myLim2 = myLim;
        myLim2.erase(myLim2.begin());
        
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < m; ++j)
                rowVec[j] = matRes(i, j);
            
            testVal = myFun1(rowVec, uM);
            Success = myComp(testVal, myLim) || myComp2(testVal, myLim2);
            if (Success) {
                indexMatch.push_back(i);
                matchVals.push_back(testVal);
            }
        }
    }
    
    int numCols = m;
    if (keepRes) {++numCols;}
    int numMatches = indexMatch.size();
    
    if (bLower)
        nRows = numMatches;
    else
        nRows  = (numMatches > userRows) ? userRows : numMatches;
    
    typeRcpp returnMatrix = Rcpp::no_init_matrix(nRows, numCols);
    
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < m; ++j)
            returnMatrix(i,j) = matRes(indexMatch[i],j);
        
        if (keepRes)
            returnMatrix(i,m) = matchVals[i];
    }
    
    return returnMatrix;
}

template <typename typeVector>
SEXP ApplyFunction(int n, int m, typeVector sexpVec, bool IsRep, int nRows, bool IsComb,
                       std::vector<int> myReps, std::vector<int> freqs, std::vector<int> z,
                       bool permNonTriv, bool IsMultiset, SEXP stdFun, SEXP rho) {
    if (IsComb) {
        if (IsMultiset)
            return Combinations::MultisetComboApplyFun(n, m, sexpVec, myReps, freqs, nRows, z, stdFun, rho);
        else
            return Combinations::ComboGeneralApplyFun(n , m, sexpVec, IsRep, nRows, z, stdFun, rho);
    } else {
        return Permutations::PermutationApplyFun(n, m, sexpVec, IsRep,nRows, IsMultiset, z, stdFun, rho);
    }
}


// [[Rcpp::export]]
SEXP CombinatoricsRcpp(SEXP Rv, SEXP Rm, SEXP Rrepetition, SEXP RFreqs, SEXP Rlow, SEXP Rhigh,
                       SEXP f1, SEXP f2, SEXP lim, bool IsComb, SEXP RKeepRes, bool IsFactor,
                       bool IsCount, SEXP stdFun, SEXP myEnv) {
    
    int n, m = 0, m1, m2;
    int lenFreqs = 0, nRows = 0;
    bool IsRepetition, IsLogical, keepRes;
    bool IsMultiset, IsInteger, IsCharacter;
    
    std::vector<double> vNum;
    std::vector<int> vInt, myReps, freqsExpanded;
    Rcpp::CharacterVector vStr;
    keepRes = Rcpp::as<bool>(RKeepRes);
    
    switch(TYPEOF(Rv)) {
        case LGLSXP: {
            IsLogical = true;
            keepRes = IsInteger = IsCharacter = false;
            break;
        }
        case INTSXP: {
            IsInteger = true;
            IsLogical = IsCharacter = false;
            break;
        }
        case REALSXP: {
            IsLogical = IsInteger = IsCharacter = false;
            break;
        }
        case STRSXP: {
            IsCharacter = true;
            keepRes = IsLogical = IsInteger = false;
            break;
        }
        default: {
            Rcpp::stop("Only integers, numerical, character, and factor classes are supported for v");   
        }
    }
    
    if (Rf_isNull(RFreqs)) {
        IsMultiset = false;
        myReps.push_back(1);
    } else {
        IsMultiset = true;
        CleanConvert::convertVector(RFreqs, myReps, "freqs must be of type numeric or integer");
        
        lenFreqs = (int) myReps.size();
        for (int i = 0; i < lenFreqs; ++i) {
            if (myReps[i] < 1) 
                Rcpp::stop("each element in freqs must be a positive number");
            
            for (int j = 0; j < myReps[i]; ++j)
                freqsExpanded.push_back(i);
        }
    }
    
    if (Rf_isNull(Rm)) {
        if (IsMultiset) {
            m = freqsExpanded.size();
        } else {
            Rcpp::stop("m and freqs cannot both be NULL");
        }
    } else {
        if (Rf_length(Rm) > 1)
            Rcpp::stop("length of m must be 1");
        
        CleanConvert::convertPrimitive(Rm, m, "m must be of type numeric or integer");
    }
    
    if (m < 1)
        Rcpp::stop("m must be positive");
    
    std::vector<double> rowVec(m);
    
    if (!Rf_isLogical(Rrepetition))
        Rcpp::stop("repetitions must be a logical value");
    
    IsRepetition = Rcpp::as<bool>(Rrepetition);
    
    if (!Rf_isLogical(RKeepRes))
        Rcpp::stop("keepResults must be a logical value");
    
    double seqEnd;
    
    if (IsCharacter) {
        vStr = Rcpp::as<Rcpp::CharacterVector>(Rv);
        n = vStr.size();
    } else if (IsLogical) {
        vInt = Rcpp::as<std::vector<int> >(Rv);
        n = vInt.size();
    } else {
        if (Rf_length(Rv) == 1) {
            seqEnd = Rcpp::as<double>(Rv);
            if (Rcpp::NumericVector::is_na(seqEnd)) {seqEnd = 1;}
            if (seqEnd > 1) {m1 = 1; m2 = seqEnd;} else {m1 = seqEnd; m2 = 1;}
            Rcpp::IntegerVector vTemp = Rcpp::seq(m1, m2);
            IsInteger = true;
            vNum = Rcpp::as<std::vector<double> >(vTemp);
        } else {
            vNum = Rcpp::as<std::vector<double> >(Rv);
        }
        
        n = vNum.size();
    }
    
    if (IsInteger) {
        for (int i = 0; i < n && IsInteger; ++i)
            if (Rcpp::NumericVector::is_na(vNum[i]))
                IsInteger = false;
        
        if (IsInteger)
            vInt.assign(vNum.begin(), vNum.end());
    }
        
    bool IsConstrained;
    if (IsFactor)
        keepRes = IsConstrained = IsCharacter = IsInteger = false;
    
    if (Rf_isNull(lim)) {
        IsConstrained = false;
    } else {
        if (IsCharacter || Rf_isNull(f1)) {
            IsConstrained = false;
        } else {
            if (!Rf_isString(f1))
                Rcpp::stop("constraintFun must be passed as a character");
            
            if (Rf_isNull(f2)) {
                IsConstrained = false;
            } else if (IsLogical) {
                IsConstrained = false;
            } else {
                IsConstrained = true;
                if (!Rf_isString(f2))
                    Rcpp::stop("comparisonFun must be passed as a character");
            }
        }
    }
    
    if (IsConstrained) {
        for (int i = (vNum.size() - 1); i >= 0; --i)
            if (Rcpp::NumericVector::is_na(vNum[i]))
                vNum.erase(vNum.begin() + i);
            
        n = vNum.size();
    }
    
    double computedRows;
    
    if (IsMultiset) {
        if (n != lenFreqs)
            Rcpp::stop("the length of freqs must equal the length of v");
        
        if (m > (int) freqsExpanded.size())
            m = freqsExpanded.size();
        
        if (IsComb) {
            computedRows = MultisetCombRowNum(n, m, myReps);
        } else {
            if (Rf_isNull(Rm)) {
                computedRows = NumPermsWithRep(freqsExpanded);
            } else if (m == (int) freqsExpanded.size()) {
                computedRows = NumPermsWithRep(freqsExpanded);
            } else {
                computedRows = MultisetPermRowNum(n, m, myReps);
            }
        }
    } else {
        if (IsRepetition) {
            if (IsComb)
                computedRows = NumCombsWithRep(n, m);
            else
                computedRows = std::pow((double) n, (double) m);
        } else {
            if (m > n)
                Rcpp::stop("m must be less than or equal to the length of v");
            
            if (IsComb)
                computedRows = nChooseK(n, m);
            else
                computedRows = NumPermsNoRep(n, m);
        }
    }
    
    mpz_t computedRowMpz;
    mpz_init(computedRowMpz);
    bool isGmp = false;
    
    if (computedRows > Significand53) {
        isGmp = true;
        if (IsMultiset) {
            if (IsComb) {
                MultisetCombRowNumGmp(computedRowMpz, n, m, myReps);
            } else {
                if (Rf_isNull(Rm)) {
                    NumPermsWithRepGmp(computedRowMpz, freqsExpanded);
                } else if (m == (int) freqsExpanded.size()) {
                    NumPermsWithRepGmp(computedRowMpz, freqsExpanded);
                } else {
                    MultisetPermRowNumGmp(computedRowMpz, n, m, myReps);
                }
            }
        } else {
            if (IsRepetition) {
                if (IsComb) {
                    NumCombsWithRepGmp(computedRowMpz, n, m);
                } else {
                    mpz_ui_pow_ui(computedRowMpz, n, m);
                }
            } else {
                if (IsComb)
                    nChooseKGmp(computedRowMpz, n, m);
                else
                    NumPermsNoRepGmp(computedRowMpz, n, m);
            }
        }
    }
    
    double lower = 0, upper = 0;
    bool bLower = false;
    bool bUpper = false;
    mpz_t lowerMpz[1], upperMpz[1];
    mpz_init(lowerMpz[0]); mpz_init(upperMpz[0]);
    mpz_set_ui(lowerMpz[0], 0); mpz_set_ui(upperMpz[0], 0);
    
    if (!IsCount) {
        if (isGmp) {
            if (!Rf_isNull(Rlow)) {
                bLower = true;
                createMPZArray(Rlow, lowerMpz, 1);
                mpz_sub_ui(lowerMpz[0], lowerMpz[0], 1);
                
                if (mpz_cmp_ui(lowerMpz[0], 0) < 0)
                    Rcpp::stop("bounds must be positive");
            }
            
            if (!Rf_isNull(Rhigh)) {
                bUpper = true;
                createMPZArray(Rhigh, upperMpz, 1);
                
                if (mpz_cmp_ui(upperMpz[0], 0) < 0)
                    Rcpp::stop("bounds must be positive");
            }
        } else { 
            if (!Rf_isNull(Rlow)) {
                bLower = true;
                CleanConvert::convertPrimitive(Rlow, lower, "bounds must be of type numeric or integer");
                --lower;
                
                if (lower < 0)
                    Rcpp::stop("bounds must be positive");
            }
            
            if (!Rf_isNull(Rhigh)) {
                bUpper = true;
                CleanConvert::convertPrimitive(Rhigh, upper, "bounds must be of type numeric or integer");
                
                if (upper < 0)
                    Rcpp::stop("bounds must be positive");
            }
        }
    }
    
    if (isGmp) {
        if (mpz_cmp(lowerMpz[0], computedRowMpz) >= 0 || mpz_cmp(upperMpz[0], computedRowMpz) > 0)
            Rcpp::stop("bounds cannot exceed the maximum number of possible results");
    } else {
        if (lower >= computedRows || upper > computedRows)
            Rcpp::stop("bounds cannot exceed the maximum number of possible results");
    }
    
    if (IsCount) {
        if (isGmp) {
            unsigned long int size;
            unsigned long int numb = 8 * sizeof(int);
            size = sizeof(int) * (2 + (mpz_sizeinbase(computedRowMpz, 2) + numb - 1) / numb);
            
            SEXP ansPos = PROTECT(Rf_allocVector(RAWSXP, size));
            char* rPos = (char*)(RAW(ansPos));
            ((int*)(rPos))[0] = 1; // first int is vector-size-header
            
            // current position in rPos[] (starting after vector-size-header)
            unsigned long int posPos = sizeof(int);
            posPos += myRaw(&rPos[posPos], computedRowMpz, size);
            
            Rf_setAttrib(ansPos, R_ClassSymbol, Rf_mkString("bigz"));
            UNPROTECT(1);
            return(ansPos);
        } else {
            return Rcpp::wrap(computedRows);
        }
    }
    
    std::vector<int> startZ(m);
    bool permNonTrivial = false;
    if (!isGmp)
        mpz_set_d(lowerMpz[0], lower);
    
    if (bLower && mpz_cmp_ui(lowerMpz[0], 0) > 0) {
        if (IsComb) {
            if (isGmp)
                startZ = nthCombinationGmp(n, m, lowerMpz[0], IsRepetition, IsMultiset, myReps);
            else
                startZ = nthCombination(n, m, lower, IsRepetition, IsMultiset, myReps);
        } else {
            permNonTrivial = true;
            if (isGmp)
                startZ = nthPermutationGmp(n, m, lowerMpz[0], IsRepetition, IsMultiset, myReps);
            else
                startZ = nthPermutation(n, m, lower, IsRepetition, IsMultiset, myReps);
            
            if (IsMultiset) {
                
                for (std::size_t j = 0; j < startZ.size(); ++j) {
                    for (std::size_t i = 0; i < freqsExpanded.size(); ++i) {
                        if (freqsExpanded[i] == startZ[j]) {
                            freqsExpanded.erase(freqsExpanded.begin() + i);
                            break;
                        }
                    }
                }
                
                for (std::size_t i = 0; i < freqsExpanded.size(); ++i)
                    startZ.push_back(freqsExpanded[i]);
                
            } else if (!IsRepetition) {
                
                if (m < n) {
                    for (int i = 0; i < n; ++i) {
                        bool bExist = false;
                        for (std::size_t j = 0; j < startZ.size(); ++j) {
                            if (startZ[j] == i) {
                                bExist = true;
                                break;
                            }
                        }
                        if (!bExist)
                            startZ.push_back(i);
                    }
                }
            }
        }
    } else {
        if (IsComb) {
            if (IsMultiset)
                startZ.assign(freqsExpanded.begin(), freqsExpanded.begin() + m);
            else if (IsRepetition)
                std::fill(startZ.begin(), startZ.end(), 0);
            else
                std::iota(startZ.begin(), startZ.end(), 0);
        } else {
            if (IsMultiset) {
                startZ = freqsExpanded;
            } else if (IsRepetition) {
                std::fill(startZ.begin(), startZ.end(), 0);
            } else {
                startZ.resize(n);
                std::iota(startZ.begin(), startZ.end(), 0);
            }
        }
    }
    
    double userNumRows = 0;
    
    if (isGmp) {
        if (bLower && bUpper) {
            mpz_sub(upperMpz[0], upperMpz[0], lowerMpz[0]);
            mpz_abs(lowerMpz[0], upperMpz[0]);
            if (mpz_cmp_ui(lowerMpz[0], INT_MAX) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(upperMpz[0]);
        } else if (bUpper) {
            permNonTrivial = true;
            if (mpz_cmp_d(upperMpz[0], INT_MAX) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
                
            userNumRows = mpz_get_d(upperMpz[0]);
        } else if (bLower) {
            mpz_sub(computedRowMpz, computedRowMpz, lowerMpz[0]);
            mpz_abs(lowerMpz[0], computedRowMpz);
            if (mpz_cmp_d(lowerMpz[0], INT_MAX) > 0)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = mpz_get_d(computedRowMpz);
        }
    } else {
        if (bLower && bUpper)
            userNumRows = upper - lower;
        else if (bUpper)
            userNumRows = upper;
        else if (bLower)
            userNumRows = computedRows - lower;
    }
    
    mpz_clear(lowerMpz[0]); mpz_clear(upperMpz[0]);
    mpz_clear(computedRowMpz);
    
    if (userNumRows == 0) {
        if (bLower && bUpper) {
            // Since lower is decremented and upper isn't
            // upper - lower = 0 means that lower is one
            // larger than upper as put in by the user
            
            Rcpp::stop("The number of rows must be positive. Either the lowerBound "
                           "exceeds the maximum number of possible results or the "
                           "lowerBound is greater than the upperBound.");
        } else {
            if (computedRows > INT_MAX)
                Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
            
            userNumRows = computedRows;
            nRows = (int) computedRows;
        }
    } else if (userNumRows < 0) {
        Rcpp::stop("The number of rows must be positive. Either the lowerBound "
                  "exceeds the maximum number of possible results or the "
                  "lowerBound is greater than the upperBound.");
    } else if (userNumRows > INT_MAX) {
        Rcpp::stop("The number of rows cannot exceed 2^31 - 1.");
    } else {
        nRows = (int) userNumRows;
    }
    
    unsigned long int uM = m;
    
    if (IsConstrained) {
        std::vector<double> myLim;
        CleanConvert::convertVector(lim, myLim, "limitConstraints must be of type numeric or integer");
        
        if (myLim.size() > 2)
            Rcpp::stop("there cannot be more than 2 limitConstraints");
        else if (myLim.size() == 2)
            if (myLim[0] == myLim[1])
                Rcpp::stop("The limitConstraints must be different");
        
        for (std::size_t i = 0; i < myLim.size(); ++i)
            if (Rcpp::NumericVector::is_na(myLim[i]))
                Rcpp::stop("limitConstraints cannot be NA");
        
        std::string mainFun1 = Rcpp::as<std::string >(f1);
        if (mainFun1 != "prod" && mainFun1 != "sum" && mainFun1 != "mean"
                && mainFun1 != "max" && mainFun1 != "min") {
            Rcpp::stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
        }
        
        if (mainFun1 == "mean")
            IsInteger = false;
        
        std::vector<std::string>::const_iterator itComp;
        std::vector<std::string> compFunVec = Rcpp::as<std::vector<std::string> >(f2);
        
        if (compFunVec.size() > 2)
            Rcpp::stop("there cannot be more than 2 comparison operators");
        
        for (std::size_t i = 0; i < compFunVec.size(); ++i) {
            itComp = std::find(compForms.begin(), 
                                compForms.end(), compFunVec[i]);
            
            if (itComp == compForms.end())
                Rcpp::stop("comparison operators must be one of the following: >, >=, <, <=, or ==");
                
            int myIndex = std::distance(compForms.begin(), itComp);
            
            // The first 5 are "standard" whereas the 6th and 7th
            // are written with the equality first. Converting
            // them here makes it easier to deal with later.
            if (myIndex > 4)
                myIndex -= 3;
            
            compFunVec[i] = compForms[myIndex];
        }
        
        if (compFunVec.size() == 2) {
            if (myLim.size() == 1) {
                compFunVec.pop_back();
            } else {
                if (compFunVec[0] == "==" || compFunVec[1] == "==")
                    Rcpp::stop("If comparing against two limitConstraints, the "
                         "equality comparisonFun (i.e. '==') cannot be used. "
                         "Instead, use '>=' or '<='.");
                
                if (compFunVec[0].substr(0, 1) == compFunVec[1].substr(0, 1))
                    Rcpp::stop("Cannot have two 'less than' comparisonFuns or two 'greater than' "
                          "comparisonFuns (E.g. c('<', '<=') is not allowed).");
                
                bool sortNeeded = false;
                
                if (compFunVec[0].substr(0, 1) == ">" && 
                        std::min(myLim[0], myLim[1]) == myLim[0]) {
                    
                    compFunVec[0] = compFunVec[0] + "," + compFunVec[1];
                    sortNeeded = true;
                    
                } else if (compFunVec[0].substr(0, 1) == "<" && 
                    std::max(myLim[0], myLim[1]) == myLim[0]) {
                    
                    compFunVec[0] = compFunVec[1] + "," + compFunVec[0];
                    sortNeeded = true;
                }
                
                if (sortNeeded) {
                    if (std::max(myLim[0], myLim[1]) == myLim[1]) {
                        double temp = myLim[0];
                        myLim[0] = myLim[1];
                        myLim[1] = temp;
                    }
                    compFunVec.pop_back();
                }
            }
        }
        
        bool SpecialCase = false;
        
        // If bLower, the user is looking to test a particular range.
        // Otherwise, the constraint algorithm will simply return (upper
        // - lower) # of combinations/permutations that meet the criteria
        if (bLower) {
            SpecialCase = true;
        } else if (mainFun1 == "prod") {
            for (int i = 0; i < n; ++i) {
                if (vNum[i] < 0) {
                    SpecialCase = true;
                    break;
                }
            }
        }
        
        std::vector<int> limInt(myLim.begin(), myLim.end());
        
        if (SpecialCase) {
            if (IsInteger) {
                return SpecCaseRet<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, nRows, keepRes, startZ,
                                                        mainFun1, IsMultiset, computedRows, compFunVec, limInt,
                                                        IsComb, myReps, freqsExpanded, bLower, permNonTrivial, userNumRows);
            } else {
                return SpecCaseRet<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes, startZ,
                                                        mainFun1, IsMultiset, computedRows, compFunVec, myLim,
                                                        IsComb, myReps, freqsExpanded, bLower, permNonTrivial, userNumRows);
            }
        }
        
        if (keepRes) {
            if (IsInteger) {
                std::vector<int> rowVecInt(m);
                Rcpp::XPtr<funcPtr<int> > xpIntFun1 = putFunPtrInXPtr<int>(mainFun1);
                funcPtr<int> myIntFun = *xpIntFun1;
                Rcpp::IntegerMatrix matRes = CombinatoricsConstraints<Rcpp::IntegerMatrix>(n, m, vNum, IsRepetition,
                                                                                           mainFun1, compFunVec, myLim, nRows,
                                                                                           IsComb, true, myReps, IsMultiset);
                nRows = matRes.nrow();
                
                for (int i = 0; i < nRows; ++i) {
                    for (int j = 0; j < m; ++j)
                        rowVecInt[j] = matRes(i, j);
                    
                    matRes(i, m) = myIntFun(rowVecInt, uM);
                }
                
                return matRes;
            }
            
            Rcpp::XPtr<funcPtr<double> > xpDblFun1 = putFunPtrInXPtr<double>(mainFun1);
            funcPtr<double> myDblFun = *xpDblFun1;
            Rcpp::NumericMatrix matRes = CombinatoricsConstraints<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition,
                                                                                       mainFun1, compFunVec, myLim, nRows,
                                                                                       IsComb, true, myReps, IsMultiset);
            nRows = matRes.nrow();
            
            for (int i = 0; i < nRows; ++i) {
                for (int j = 0; j < m; ++j)
                    rowVec[j] = matRes(i, j);
                
                matRes(i, m) = myDblFun(rowVec, uM);
            }
            
            return matRes;
        }
        
        if (IsInteger) {
            return CombinatoricsConstraints<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, mainFun1,
                                                                 compFunVec, limInt, nRows, IsComb,
                                                                 false, myReps, IsMultiset);
        }
        
        return CombinatoricsConstraints<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, mainFun1,
                                                                 compFunVec, myLim, nRows, IsComb,
                                                                 false, myReps, IsMultiset);
    } else {
        bool applyFun = !Rf_isNull(stdFun) && !IsFactor;

        if (applyFun) {
            if (!Rf_isFunction(stdFun))
                Rcpp::stop("FUN must be a function!");
            
            if (IsCharacter) {
                Rcpp::CharacterVector rcppVChar(vStr.begin(), vStr.end());
                return ApplyFunction(n, m, rcppVChar, IsRepetition, nRows, IsComb, myReps,
                                     freqsExpanded, startZ, permNonTrivial,IsMultiset, stdFun, myEnv);
            } else if (IsLogical || IsInteger) {
                Rcpp::IntegerVector rcppVInt(vInt.begin(), vInt.end());
                return ApplyFunction(n, m, rcppVInt, IsRepetition, nRows, IsComb, myReps,
                                     freqsExpanded, startZ, permNonTrivial,IsMultiset, stdFun, myEnv);
            } else {
                Rcpp::NumericVector rcppVNum(vNum.begin(), vNum.end());
                return ApplyFunction(n, m, rcppVNum, IsRepetition, nRows, IsComb, myReps,
                                     freqsExpanded, startZ, permNonTrivial,IsMultiset, stdFun, myEnv);
            }
        }
        
        if (Rf_isNull(f1) || IsCharacter || IsFactor)
            keepRes = false;
        
        std::string mainFun2;
        funcPtr<double> myFun2;
        
        if (keepRes) {
            mainFun2 = Rcpp::as<std::string >(f1);
            if (mainFun2 != "prod" && mainFun2 != "sum" && mainFun2 != "mean"
                    && mainFun2 != "max" && mainFun2 != "min") {
                Rcpp::stop("contraintFun must be one of the following: prod, sum, mean, max, or min");
            }
            
            Rcpp::XPtr<funcPtr<double> > xpFun2 = putFunPtrInXPtr<double>(mainFun2);
            myFun2 = *xpFun2;
            
            // If our function operating on the rows of our matrix
            // can possibly produce elements greater than INT_MAX.
            // We need a NumericMatrix in this case. We also
            // need to check if our function is the mean that
            // can return not integer values.
            if (IsInteger) {
                if (mainFun2 == "mean") {
                    IsInteger = false;
                } else {
                    int vecMax = *std::max_element(vInt.begin(), vInt.end());
                    
                    for (int i = 0; i < m; ++i)
                        rowVec[i] = (double) vecMax;
                    
                    double testIfInt = myFun2(rowVec, uM);
                    if (testIfInt > INT_MAX)
                        IsInteger = false;
                }
            }
        }
        
        if (IsCharacter) {
            if (IsComb) {
                if (IsMultiset)
                    return Combinations::MultisetCombination<Rcpp::CharacterMatrix>(n, m, vStr, myReps, 
                                                                                    freqsExpanded, nRows, keepRes, startZ);
                else
                    return Combinations::ComboGeneral<Rcpp::CharacterMatrix>(n , m, vStr, IsRepetition, nRows, keepRes, startZ);
            } else {
                if (IsMultiset)
                    return Permutations::MultisetPermutation<Rcpp::CharacterMatrix>(n, m, vStr, nRows, keepRes, startZ);
                else
                    return Permutations::PermuteGeneral<Rcpp::CharacterMatrix>(n, m, vStr, IsRepetition, nRows, keepRes, startZ, permNonTrivial);
            }
        } else if (IsLogical) {
            if (IsComb) {
                if (IsMultiset)
                    return Combinations::MultisetCombination<Rcpp::LogicalMatrix>(n, m, vInt, myReps, 
                                                                                  freqsExpanded, nRows, keepRes, startZ);
                else
                    return Combinations::ComboGeneral<Rcpp::LogicalMatrix>(n , m, vInt, IsRepetition, nRows, keepRes, startZ);
            } else {
                if (IsMultiset)
                    return Permutations::MultisetPermutation<Rcpp::LogicalMatrix>(n, m, vInt, nRows, keepRes, startZ);
                else
                    return Permutations::PermuteGeneral<Rcpp::LogicalMatrix>(n, m, vInt, IsRepetition, nRows, keepRes, startZ, permNonTrivial);
            }
        } else if (IsFactor) {
            Rcpp::IntegerMatrix factorMat;
            Rcpp::IntegerVector testFactor = Rcpp::as<Rcpp::IntegerVector>(Rv);
            Rcpp::CharacterVector myClass = testFactor.attr("class");
            Rcpp::CharacterVector myLevels = testFactor.attr("levels");
            
            if (IsComb) {
                if (IsMultiset)
                    factorMat = Combinations::MultisetCombination<Rcpp::IntegerMatrix>(n, m, vInt, myReps, 
                                                                                       freqsExpanded, nRows, keepRes, startZ);
                else
                    factorMat = Combinations::ComboGeneral<Rcpp::IntegerMatrix>(n , m, vInt, IsRepetition, nRows, keepRes, startZ);
            } else {
                if (IsMultiset)
                    factorMat = Permutations::MultisetPermutation<Rcpp::IntegerMatrix>(n, m, vInt, nRows, keepRes, startZ);
                else
                    factorMat = Permutations::PermuteGeneral<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, nRows, keepRes, startZ, permNonTrivial);
            }
            
            factorMat.attr("class") = myClass;
            factorMat.attr("levels") = myLevels;
            return factorMat;
            
        } else if (IsInteger) {
            Rcpp::IntegerMatrix matResInt;
            
            if (IsComb) {
                if (IsMultiset)
                    matResInt = Combinations::MultisetCombination<Rcpp::IntegerMatrix>(n, m, vInt, myReps, 
                                                                                       freqsExpanded, nRows, keepRes, startZ);
                else
                    matResInt = Combinations::ComboGeneral<Rcpp::IntegerMatrix>(n , m, vInt, IsRepetition, nRows, keepRes, startZ);
            } else {
                if (IsMultiset)
                    matResInt = Permutations::MultisetPermutation<Rcpp::IntegerMatrix>(n, m, vInt, nRows, keepRes, startZ);
                else
                    matResInt = Permutations::PermuteGeneral<Rcpp::IntegerMatrix>(n, m, vInt, IsRepetition, nRows, keepRes, startZ, permNonTrivial);
            }
            
            if (keepRes) {
                funcPtr<int> myFunInt;
                Rcpp::XPtr<funcPtr<int> > xpFunInt2 = putFunPtrInXPtr<int>(mainFun2);
                myFunInt = *xpFunInt2;
                std::vector<int> rowVecInt(m);
                
                for (int i = 0; i < nRows; ++i) {
                    for (int j = 0; j < m; ++j)
                        rowVecInt[j] = matResInt(i, j);
                    
                    matResInt(i, m) = myFunInt(rowVecInt, uM);
                }
            }
            
            return matResInt;
            
        } else {
            Rcpp::NumericMatrix matResDbl;
            
            if (IsComb) {
                if (IsMultiset)
                    matResDbl = Combinations::MultisetCombination<Rcpp::NumericMatrix>(n, m, vNum, myReps, 
                                                                                       freqsExpanded, nRows, keepRes, startZ);
                else
                    matResDbl = Combinations::ComboGeneral<Rcpp::NumericMatrix>(n , m, vNum, IsRepetition, nRows, keepRes, startZ);
            } else {
                if (IsMultiset)
                    matResDbl = Permutations::MultisetPermutation<Rcpp::NumericMatrix>(n, m, vNum, nRows, keepRes, startZ);
                else
                    matResDbl = Permutations::PermuteGeneral<Rcpp::NumericMatrix>(n, m, vNum, IsRepetition, nRows, keepRes, startZ, permNonTrivial);
            }
            
            if (keepRes) {
                funcPtr<double> myFunDbl;
                Rcpp::XPtr<funcPtr<double> > xpFunDbl2 = putFunPtrInXPtr<double>(mainFun2);
                myFunDbl = *xpFunDbl2;
                
                for (int i = 0; i < nRows; ++i) {
                    for (int j = 0; j < m; ++j)
                        rowVec[j] = matResDbl(i, j);
                    
                    matResDbl(i, m) = myFunDbl(rowVec, uM);
                }
            }
            
            return matResDbl;
        }
    }
}
