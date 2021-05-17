## See src/CheckReturn.cpp for more details
comboGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL,
                         upper = NULL, constraintFun = NULL, comparisonFun = NULL,
                         limitConstraints = NULL, keepResults = NULL, FUN = NULL,
                         Parallel = FALSE, nThreads = NULL, tolerance = NULL) {
    
    IsFactor <- is.factor(v)
    RetValue <- CheckReturn(v, constraintFun, comparisonFun,
                            limitConstraints, IsFactor, keepResults, FUN)

    if (RetValue == 1) {
        CombinatoricsStndrd(v, m, repetition, freqs, lower, upper, TRUE,
                            IsFactor, Parallel, nThreads, pkgEnv$nThreads)
    } else if (RetValue == 2) {
        CombinatoricsApply(v, m, repetition, freqs, 
                           lower, upper, TRUE, FUN, new.env())
    } else {
        CombinatoricsCnstrt(v, m, repetition, freqs, lower, upper, constraintFun,
                            comparisonFun, limitConstraints, TRUE, keepResults,
                            Parallel, nThreads, pkgEnv$nThreads, tolerance)
    }
}

permuteGeneral <- function(v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL,
                           upper = NULL, constraintFun = NULL, comparisonFun = NULL, 
                           limitConstraints = NULL, keepResults = NULL, FUN = NULL,
                           Parallel = FALSE, nThreads = NULL, tolerance = NULL) {

    IsFactor <- is.factor(v)
    RetValue <- CheckReturn(v, constraintFun, comparisonFun,
                            limitConstraints, IsFactor, keepResults, FUN)
    
    if (RetValue == 1) {
        CombinatoricsStndrd(v, m, repetition, freqs, lower, upper, FALSE,
                            IsFactor, Parallel, nThreads, pkgEnv$nThreads)
    } else if (RetValue == 2) {
        CombinatoricsApply(v, m, repetition, freqs, 
                           lower, upper, FALSE, FUN, new.env())
    } else {
        CombinatoricsCnstrt(v, m, repetition, freqs, lower, upper, constraintFun,
                            comparisonFun, limitConstraints, FALSE, keepResults,
                            Parallel, nThreads, pkgEnv$nThreads, tolerance)
    }
}

comboCount <-  function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    CombinatoricsCount(v, m, repetition, freqs, TRUE);
}

permuteCount <- function(v, m = NULL, repetition = FALSE, freqs = NULL) {
    CombinatoricsCount(v, m, repetition, freqs, FALSE);
}

comboSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL, 
                        sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
                        nThreads = NULL, namedSample = FALSE) {
    
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, TRUE, is.factor(v), seed, n,
               sample, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, namedSample)
}

permuteSample <- function(v, m = NULL, repetition = FALSE, freqs = NULL, n = NULL,
                          sampleVec = NULL, seed = NULL, FUN = NULL, Parallel = FALSE,
                          nThreads = NULL, namedSample = FALSE) {
    
    if (!is.null(seed)) {set.seed(seed)}
    SampleRcpp(v, m, repetition, freqs, sampleVec, FALSE, is.factor(v), seed, n,
               sample, FUN, new.env(), Parallel, nThreads, pkgEnv$nThreads, namedSample)
}

comboGroups <- function(v, numGroups, retType = "matrix", lower = NULL,
                       upper = NULL, Parallel = FALSE, nThreads = NULL) {
    
    ComboGroupsRcpp(v, numGroups, retType, lower, upper, is.factor(v),
                    Parallel, nThreads, pkgEnv$nThreads, FALSE, NULL, 
                    NULL, NULL, sample, FALSE)
}

comboGroupsSample <- function(v, numGroups, retType = "matrix", n = NULL,
                              sampleVec = NULL, seed = NULL, Parallel = FALSE, 
                              nThreads = NULL, namedSample = FALSE) {
    
    if (!is.null(seed)) {set.seed(seed)}
    ComboGroupsRcpp(v, numGroups, retType, NULL, NULL, is.factor(v),
                    Parallel, nThreads, pkgEnv$nThreads, TRUE, sampleVec,
                    seed, n, sample, namedSample)
}

comboGroupsCount <- function(v, numGroups) {
    ComboGroupsCountCpp(v, numGroups)
}

