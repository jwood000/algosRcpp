#include "ComboFunClass.h"

template <int RTYPE>
SEXP ComboFUN::VecFUNRes(const Rcpp::Vector<RTYPE> &v) {

    Rcpp::Vector<RTYPE> vectorPass(m);
    SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));

    for (int j = 0; j < m; ++j)
        vectorPass[j] = v[z[j]];
    
    SETCADR(sexpFun, vectorPass);
    SEXP res = Rf_eval(sexpFun, rho);
    UNPROTECT(1);
    
    return res;
}

template <int RTYPE>
Rcpp::List ComboFUN::ListComboFUN(const Rcpp::Vector<RTYPE> &v, int nRows) {

    Rcpp::List myList(nRows);
    Rcpp::Vector<RTYPE> vectorPass(m);
    SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));
    
    const int loc_n1 = IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1);
    const int lastRow = nRows - 1;
    const int loc_m = m;
    
    // We iterater to the pentultimate row to avoid iterating z one too many times
    for (int count = 0, loc_m1 = m - 1; count < lastRow; ++count) {
        for (int j = 0; j < loc_m; ++j)
            vectorPass[j] = v[z[j]];
        
        SETCADR(sexpFun, vectorPass);
        myList[count] = Rf_eval(sexpFun, rho);
        nextIter(freqs, z, loc_n1, loc_m1);
    }
    
    // Get the last result
    for (int j = 0; j < loc_m; ++j)
        vectorPass[j] = v[z[j]];
    
    SETCADR(sexpFun, vectorPass);
    myList[lastRow] = Rf_eval(sexpFun, rho);
    UNPROTECT(1);
    
    return myList;
}

template <int RTYPE>
Rcpp::List ComboFUN::ListComboFUNReverse(const Rcpp::Vector<RTYPE> &v, int nRows) {

    Rcpp::List myList(nRows);
    Rcpp::Vector<RTYPE> vectorPass(m);
    SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));
    
    const int loc_n1 = IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1);
    const int lastRow = nRows - 1;
    const int loc_m = m;

    // We iterater to the pentultimate row to avoid iterating z one too many times
    for (int count = 0, loc_m1 = m - 1; count < lastRow; ++count) {
        for (int j = 0; j < loc_m; ++j)
            vectorPass[j] = v[z[j]];
        
        SETCADR(sexpFun, vectorPass);
        myList[count] = Rf_eval(sexpFun, rho);
        prevIter(freqs, z, loc_n1, loc_m1);
    }

    // Get the last result
    for (int j = 0; j < loc_m; ++j)
        vectorPass[j] = v[z[j]];
    
    SETCADR(sexpFun, vectorPass);
    myList[lastRow] = Rf_eval(sexpFun, rho);
    UNPROTECT(1);
    
    return myList;
}

template <int RTYPE>
Rcpp::List ComboFUN::SampList(const Rcpp::Vector<RTYPE> &v,
                              const std::vector<double> &mySample,
                              mpz_t *const myBigSamp, std::size_t sampSize) {
    
    Rcpp::List myList(sampSize);
    Rcpp::Vector<RTYPE> vectorPass(m);
    SEXP sexpFun = PROTECT(Rf_lang2(stdFun, R_NilValue));

    if (IsGmp) {
        constexpr double dblDefault = 0;

        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> tempZ = nthResFun(n, m, dblDefault, myBigSamp[i], myReps);

            for (int j = 0; j < m; ++j)
                vectorPass[j] = v[tempZ[j]];
            
            SETCADR(sexpFun, vectorPass);
            myList[i] = Rf_eval(sexpFun, rho);
        }
    } else {
        for (std::size_t i = 0; i < sampSize; ++i) {
            const std::vector<int> tempZ = nthResFun(n, m, mySample[i], mpzTemp, myReps);

            for (int j = 0; j < m; ++j)
                vectorPass[j] = v[tempZ[j]];
            
            SETCADR(sexpFun, vectorPass);
            myList[i] = Rf_eval(sexpFun, rho);
        }
    }
    
    UNPROTECT(1);
    return myList;
}

// The bVec Vector represents IsFac, IsComb, IsMult, IsRep, IsGmp, & IsFull
ComboFUN::ComboFUN(SEXP Rv, int Rm, SEXP RcompRows, Rcpp::LogicalVector bVec,
                   std::vector<std::vector<int>> freqInfo, Rcpp::List funStuff)
                   : n(Rf_length(Rv)), m(Rm), m1(Rm - 1), sexpVec(Rv),
                   IsFactor(bVec[0]), IsComb(bVec[1]), IsMult(bVec[2]),
                   IsRep(bVec[3]), IsGmp(bVec[4]),
                   computedRows(IsGmp ? 0 : Rcpp::as<double>(RcompRows)),
                   freqs(freqInfo[0]), myReps(freqInfo[1]),
                   testFactor(bVec[0] ? Rcpp::as<Rcpp::IntegerVector>(Rv) :
                                     Rcpp::IntegerVector::create(1)),
                   myClass(bVec[0] ? testFactor.attr("class") :
                                  Rcpp::CharacterVector::create("empty")),
                   myLevels(bVec[0] ? testFactor.attr("levels") :
                                   Rcpp::CharacterVector::create("empty")),
                   nthResFun(*putNthResPtrInXPtr(bVec[1], bVec[2], bVec[3], bVec[4])),
                   nextIter(*putNextIterPtrInXPtr(bVec[1], bVec[2], bVec[3], bVec[5])),
                   prevIter(*putPrevIterPtrInXPtr(bVec[1], bVec[2], bVec[3], bVec[5])),
                   n1(IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1)),
                   stdFun(funStuff[0]), rho(funStuff[1]) {

    z.resize(Rm);
    
    // Initialize trivial mpz_t value for functions which require mpz_t
    mpz_init(mpzTemp);
    mpz_init(mpzIndex[0]);
    mpz_init(computedRowsMpz[0]);

    if (IsGmp)
        createMPZArray(RcompRows, computedRowsMpz, 1, "index vector");

    dblIndex = 0;

    SetStartZ(n, m, dblIndex, 0, mpzIndex[0], IsRep,
              IsComb, IsMult, IsGmp, myReps, freqs, z, nthResFun);
}

void ComboFUN::startOver() {
    if (IsGmp) {
        mpz_set_ui(mpzIndex[0], 0u);
    } else {
        dblIndex = 0;
    }

    SetStartZ(n, m, dblIndex, 0, mpzIndex[0], IsRep,
              IsComb, IsMult, IsGmp, myReps, freqs, z, nthResFun);
}

SEXP ComboFUN::nextComb() {
    
    if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 0)) {
        increment(IsGmp, mpzIndex[0], dblIndex);
        RCPP_RETURN_VECTOR(VecFUNRes, sexpVec);
    } else if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        increment(IsGmp, mpzIndex[0], dblIndex);
        nextIter(freqs, z, n1, m1);
        RCPP_RETURN_VECTOR(VecFUNRes, sexpVec);
    } else if (CheckEqInd(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        Rcpp::Rcout << "No more results. To see the last result,"
                       " use the prevIter method(s)\n" << std::endl;
        increment(IsGmp, mpzIndex[0], dblIndex);
        return Rcpp::wrap(false);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP ComboFUN::prevComb() {

    if (CheckIndGrT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        decrement(IsGmp, mpzIndex[0], dblIndex);
        RCPP_RETURN_VECTOR(VecFUNRes, sexpVec);
    } else if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 1)) {
        decrement(IsGmp, mpzIndex[0], dblIndex);
        prevIter(freqs, z, n1, m1);
        RCPP_RETURN_VECTOR(VecFUNRes, sexpVec);
    } else if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 1)) {
        Rcpp::Rcout << "Iterator Initialized. To see the first result,"
                       " use the nextIter method(s)\n" << std::endl;
        decrement(IsGmp, mpzIndex[0], dblIndex);
        return Rcpp::wrap(false);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP ComboFUN::nextNumCombs(SEXP RNum) {

    int num;
    CleanConvert::convertPrimitive(RNum, num, "The number of results");
    
    if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        int nRows = 0;
        int numIncrement = 0;
        
        if (IsGmp) {
            mpz_sub(mpzTemp, computedRowsMpz[0], mpzIndex[0]);
            nRows = mpz_cmp_si(mpzTemp, num) < 0 ? mpz_get_si(mpzTemp) : num;
            numIncrement = mpz_cmp_si(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = computedRows - dblIndex;
            nRows = num > dblTemp ? dblTemp : num;
            numIncrement = num > dblTemp ? (nRows + 1) : nRows;
        }
        
        if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0))
            nextIter(freqs, z, n1, m1);
        
        increment(IsGmp, mpzIndex[0], dblIndex, numIncrement);
        RCPP_RETURN_VECTOR(ListComboFUN, sexpVec, nRows);
    } else if (CheckEqInd(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        Rcpp::Rcout << "No more results. To see the last result,"
                       " use the prevIter method(s)\n" << std::endl;
        increment(IsGmp, mpzIndex[0], dblIndex);
        return Rcpp::wrap(false);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP ComboFUN::prevNumCombs(SEXP RNum) {

    int num;
    CleanConvert::convertPrimitive(RNum, num, "The number of results");
    
    if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 2)) {
        int nRows = 0;
        int numDecrement = 0;
        
        if (IsGmp) {
            mpz_sub_ui(mpzTemp, mpzIndex[0], 1u);
            nRows = mpz_cmp_si(mpzTemp, num) < 0 ? mpz_get_si(mpzTemp) : num;
            numDecrement = mpz_cmp_si(mpzTemp, num) < 0 ? (nRows + 1) : nRows;
        } else {
            dblTemp = dblIndex - 1;
            nRows = num > dblTemp ? dblTemp : num;
            numDecrement = num > dblTemp ? (nRows + 1) : nRows;
        }
        
        if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows, true))
            prevIter(freqs, z, n1, m1);
        
        decrement(IsGmp, mpzIndex[0], dblIndex, numDecrement);
        RCPP_RETURN_VECTOR(ListComboFUNReverse, sexpVec, nRows);
    } else if (CheckEqSi(IsGmp, mpzIndex[0], dblIndex, 2)) {
        Rcpp::Rcout << "No more results. To see the last result,"
                       " use the prevIter method(s)\n" << std::endl;
        decrement(IsGmp, mpzIndex[0], dblIndex);
        return Rcpp::wrap(false);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP ComboFUN::nextGather() {

    if (IsGmp) {
        mpz_sub(mpzTemp, computedRowsMpz[0], mpzIndex[0]);

        if (mpz_cmp_si(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            Rcpp::stop("The number of requested rows is greater than " +
                std::to_string(std::numeric_limits<int>::max()));
        }
    } else {
        dblTemp = computedRows - dblIndex;

        if (dblTemp > std::numeric_limits<int>::max()) {
            Rcpp::stop("The number of requested rows is greater than " +
                std::to_string(std::numeric_limits<int>::max()));
        }
    }

    const int nRows = (IsGmp) ? mpz_get_si(mpzTemp) : dblTemp;

    if (nRows > 0) {
        if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0))
            nextIter(freqs, z, n1, m1);
        
        if (IsGmp) {
            mpz_add_ui(mpzIndex[0], computedRowsMpz[0], 1u);
        } else {
            dblIndex = computedRows + 1;
        }
        
        RCPP_RETURN_VECTOR(ListComboFUN, sexpVec, nRows);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP ComboFUN::prevGather() {

    if (IsGmp) {
        mpz_sub_ui(mpzTemp, mpzIndex[0], 1);

        if (mpz_cmp_si(mpzTemp, std::numeric_limits<int>::max()) > 0) {
            Rcpp::stop("The number of requested rows is greater than " +
                std::to_string(std::numeric_limits<int>::max()));
        }
    } else {
        dblTemp = dblIndex - 1;

        if (dblTemp > std::numeric_limits<int>::max()) {
            Rcpp::stop("The number of requested rows is greater than " +
                std::to_string(std::numeric_limits<int>::max()));
        }
    }

    const int nRows = (IsGmp) ? mpz_get_si(mpzTemp) : dblTemp;

    if (nRows) {
        if (CheckIndLT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows, true))
            prevIter(freqs, z, n1, m1);
        
        if (IsGmp) {
            mpz_set_si(mpzIndex[0], 0u);
        } else {
            dblIndex = 0;
        }
        
        RCPP_RETURN_VECTOR(ListComboFUNReverse, sexpVec, nRows);
    } else {
        return Rcpp::wrap(false);
    }
}

SEXP ComboFUN::currComb() {

    if (CheckIndGrT(IsGmp, mpzIndex[0], dblIndex, computedRowsMpz[0], computedRows)) {
        Rcpp::Rcout << "No more results. To see the last result,"
                       " use the prevIter method(s)\n" << std::endl;
        return Rcpp::wrap(false);
    } else if (CheckGrTSi(IsGmp, mpzIndex[0], dblIndex, 0)) {
        RCPP_RETURN_VECTOR(VecFUNRes, sexpVec);
    } else {
        Rcpp::Rcout << "Iterator Initialized. To see the first result,"
                       " use the nextIter method(s)\n" << std::endl;
        return Rcpp::wrap(false);
    }
}

SEXP ComboFUN::combIndex(SEXP RindexVec) {
    
    std::size_t sampSize;
    std::vector<double> mySample;
    SetIndexVec(RindexVec, mySample, sampSize, IsGmp, computedRows);
    
    const std::size_t bigSampSize = (IsGmp) ? sampSize : 1;
    auto mpzVec = FromCpp14::make_unique<mpz_t[]>(bigSampSize);
    
    for (std::size_t i = 0; i < bigSampSize; ++i)
        mpz_init(mpzVec[i]);
    
    if (IsGmp)
        SetIndexVecMpz(RindexVec, mpzVec.get(), sampSize, computedRowsMpz[0]);
    
    if (sampSize > 1) {
        RCPP_RETURN_VECTOR(SampList, sexpVec, mySample, mpzVec.get(), sampSize);
    } else {
        if (IsGmp) {
            mpz_add_ui(mpzIndex[0], mpzVec[0], 1u);
            mpz_set(mpzTemp, mpzVec[0]);
        } else {
            dblIndex = mySample.front() + 1;
            dblTemp = mySample.front();
        }
        
        z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
        TopOffPartialPerm(z, myReps, n, m, IsComb, IsRep, IsMult);
        RCPP_RETURN_VECTOR(VecFUNRes, sexpVec);
    }
}

SEXP ComboFUN::front() {
    
    if (IsGmp) {
        mpz_set_ui(mpzIndex[0], 1u);
        mpz_set_ui(mpzTemp, 0u);
    } else {
        dblIndex = 1;
        dblTemp = 0;
    }
    
    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    TopOffPartialPerm(z, myReps, n, m, IsComb, IsRep, IsMult);
    RCPP_RETURN_VECTOR(VecFUNRes, sexpVec);
}

SEXP ComboFUN::back() {
    
    if (IsGmp) {
        mpz_set(mpzIndex[0], computedRowsMpz[0]);
        mpz_sub_ui(mpzTemp, computedRowsMpz[0], 1u);
    } else {
        dblIndex = computedRows;
        dblTemp = computedRows - 1;
    }
    
    z = nthResFun(n, m, dblTemp, mpzTemp, myReps);
    TopOffPartialPerm(z, myReps, n, m, IsComb, IsRep, IsMult);
    RCPP_RETURN_VECTOR(VecFUNRes, sexpVec);
}

SEXP ComboFUN::sourceVector() const {
    return sexpVec;
}

Rcpp::List ComboFUN::summary() {
    const std::string CoPerm = (IsComb) ? "Combinations " : "Permutations ";
    const std::string RepStr = (IsRep) ? "with repetition " : "";
    const std::string MultiStr = (IsMult) ? "of a multiset " : "";
    const std::string strDesc = CoPerm + RepStr + MultiStr + "of "
                                + std::to_string(n) + " choose " + std::to_string(m);
    const double dblDiff = (IsGmp) ? 0 : computedRows - dblIndex;

    if (IsGmp)
        mpz_sub(mpzTemp, computedRowsMpz[0], mpzIndex[0]);

    return Rcpp::List::create(Rcpp::Named("description") = Rcpp::wrap(strDesc),
                              Rcpp::Named("currentIndex") = GetCount(IsGmp, mpzIndex[0], dblIndex),
                              Rcpp::Named("totalResults") = GetCount(IsGmp, computedRowsMpz[0], computedRows),
                              Rcpp::Named("totalRemaining") = GetCount(IsGmp, mpzTemp, dblDiff));
}

RCPP_MODULE(ComboFUN) {
    Rcpp::class_<ComboFUN>("ComboFUN")

        .constructor<SEXP, int, SEXP, Rcpp::LogicalVector,
                     std::vector<std::vector<int>>, Rcpp::List>("constructor")
        
        .method("nextIter", &ComboFUN::nextComb, "Get the result of FUN applied to the next combination/permutation")
        .method("nextNIter", &ComboFUN::nextNumCombs, "Get the result of FUN applied to the next n combinations/permutations")
        .method("nextRemaining", &ComboFUN::nextGather, "Get the result of FUN applied to all of the remaining combinations/permutations")
        .method("currIter", &ComboFUN::currComb, "Get the result of FUN applied to the current combination/permutation")
        .method("prevIter", &ComboFUN::prevComb, "Get the result of FUN applied to the previous combination/permutation")
        .method("prevNIter", &ComboFUN::prevNumCombs, "Get the result of FUN applied to the previous n combinations/permutations")
        .method("prevRemaining", &ComboFUN::prevGather, "Get the result of FUN applied to all of the previous combinations/permutations")
        .method("startOver", &ComboFUN::startOver, "Resets the iterator")
        .method("sourceVector", &ComboFUN::sourceVector, "View the source vector")
        .method("summary", &ComboFUN::summary, "See a list of summary information about the iterator")
        .method("front", &ComboFUN::front, "Get the result of FUN applied to the first lexicographical combination/permutation")
        .method("back", &ComboFUN::back, "Get the result of FUN applied to the last lexicographical combination/permutation")
        .method("[[", &ComboFUN::combIndex, "Random access method. Pass a single index or vector of indices")
        ;
}
