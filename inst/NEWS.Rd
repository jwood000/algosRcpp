\name{NEWS}
\title{News for Package \pkg{RcppAlgos}}

\section{Changes in RcppAlgos version 2.2.1 (Release date: 2018-11-26)}{
    \itemize{
        \item Fixed bug associated with primeSieve for numbers between sqrt(1e9) and 2*3*5*7.
    }
}

\section{Changes in RcppAlgos version 2.2.0 (Release date: 2018-09-15)}{
    \itemize{
        \item Disabled \code{Parallel} argument as it was causing unpredictable errors on certain platforms. Research is ongoing to correct this for use in future versions. The development version will retain this feature.
        \item Corrected UBSAN error that caused by filling a vector of integers with numbers larger than 2^31 - 1.
    }
}

\section{Changes in RcppAlgos version 2.1.0 (Release date: 2018-09-11)}{
    \itemize{
        \item Added argument \code{Parallel} to general and sampling functions for increased gains in efficiency.
            \itemize{
                \item \code{comboGeneral(30, 10, Parallel = TRUE)}
                \item \code{permuteGeneral(12, 7, TRUE, constraintFun = "sum", Parallel = TRUE)}
            }
        \item Logical class is now preserved in combinatorial functions
        \item {Added gmp support to combinatorial functions. Now, one can accurately and quickly work with combinations/permutations of large vectors \emph{E.g.}:
            \itemize{
                \item \code{comboSample(runif(10000), 100, n = 10, seed = 42, Parallel = TRUE)}
                \item \code{permuteGeneral(factor(state.name), 20, 
                            lower = 1e15, upper = 1e15 + 1000)}
            }
        }
        \item {Added \code{FUN} argument to all combinatorial functions. Allows user to pass custom functions to be applied to combinations/permutations.}
    }
}

\section{Changes in RcppAlgos version 2.0.3 (Release date: 2018-06-30)}{
    \itemize{
        \item Corrected clang UBSAN error identified by two different unit tests. In both situations, the problem was occurring as a result of populating a vector of integers with values from a vector of doubles that contained a nan (Not-a-Number). Most information was obtained from Brian J. Knaus's blog titled : "Reproducing a clang-UBSAN issue" (https://knausb.github.io/2017/06/reproducing-a-clang-ubsan-issue/)
    }
}

\section{Changes in RcppAlgos version 2.0.2 (Release date: 2018-06-15)}{
    \itemize{
        \item Corrected divide by zero in divisorsRcpp unit test.
    }
}

\section{Changes in RcppAlgos version 2.0.1 (Release date: 2018-06-10)}{
    \itemize{
        \item Corrected spelling in DESCRIPTION
    }
}

\section{Changes in RcppAlgos version 2.0.0 (Release date: Never Released)}{
    \itemize{
        \item Changed max value and explicitly casted a few values to (int64_t) in PollardRho.cpp for efficiency while still maintaining accuracy from 2^60 to 2^62 (isPrimeRcpp is roughly 10\% faster now).
        \item Updated core permutation algorithm for greater efficiency and generality
        \item Added capability of generating specific combinations/permutations
        \item Changed arguments to comboGeneral/permuteGeneral. rowCap is now upper.
        \item Added comboSample, permuteSample, comboCount, and permuteCount
        \item Fixed bug with numDivisorSieve and divisorSieve when the lower bound was greater than 1 and less than the sqrt of the upper bound. In the previous version, the numbers in this range would have duplicated values/counts.
        \item Increased efficiency of numDivisorSieve by a factor of 2.
        \item Updated unit tests for greater coverage. See the function package_coverage from the package 'covr'.
    }
}

\section{Changes in RcppAlgos version 1.0.1 (Release date: 2018-03-09)}{
    \itemize{
        \item Corrected precision limits in documentation from 2^64 to 2^63.
        \item Changed const type in PollardRho.cpp from (int64_t) to (double) to correct "UndefinedBehaviorSanitizer"
        \item Changed examples in primeFactorizeSieve to reduce check time
        \item Added RcppAlgos-package man file.
    }
}

\section{Changes in RcppAlgos version 1.0.0 (Release date: 2018-03-07)}{
    \itemize{
        \item Added the following functions: primeFactorize (vectorized pollard rho factorization), divisorsRcpp (vectorized factoring (complete)), isPrimeRcpp (vectorized primality testing using Miller-Rabin algorithm), & primeCount (based on the primecount algorithm by Kim Walisch)
        \item Completely revamped the primeSieve function. It is now a segmented sieve of Eratosthenes with wheel factorization based on primesieve by Kim Walisch.
        \item Renamed divisorsList to divisorsSieve (reason for the major version update to 1.0.0)
        \item Renamed primeFactorizationList to primeFactorizeSieve
        \item Made the sieving functions more flexible. They are now able to generate results over a range and can also produce named objects.
        \item All number theoretic functions have been made more efficient. Some make use of the fast integer division library 'libdivide' by ridiculousfish.
    }
}

\section{Changes in RcppAlgos version 0.2.5 (Release date: 2018-01-04)}{
    \itemize{
        \item Added unit tests.
        \item Removed unnecessary files.
        \item Fixed bug in primeSieve that occurred when a number with a decimal was passed (e.g. 2.01).
        \item Adjusted accepted lower bound for numDivisorSieve, eulerPhiSieve, divisorsList, primeSieve, and primeFactorizationList.
        \item Fixed bug when non-unique elements are present with factors.
    }
}

\section{Changes in RcppAlgos version 0.2.4 (Release date: 2017-12-18)}{
    \itemize{
        \item Fixed bug that occurs when non-unique elements are present for combinations with replacement.
    }
}

\section{Changes in RcppAlgos version 0.2.3 (Release date: 2017-12-18)}{
    \itemize{
        \item Fixed segmentation fault error highlighted by valgrind check in version 0.2.2.
        \item Updated DESCRIPTION file.
    }
}

\section{Changes in RcppAlgos version 0.2.2 (Release date: 2017-12-15)}{
    \itemize{
        \item Fixed bug in constraint functions that occurred when m = 1 and the constraint limit was equal to the last element in v. It was returning a 2x1 matrix with the same value twice.  It is now correctly returning a 1x1 matrix with the correct value 1 time.
        \item Reorganized source code such that all utility functions for the combinatoric functions are now in their own file. Additionally added header for this file.
        \item All combinatoric functions can now utilize the rowCap argument. Before, rowCap only applied to the combinatorial functions with constraints.
        \item comboGeneral can now find all combinations of multisets.
        \item Both comboGeneral and permuteGeneral can utilize the argument m when dealing with multisets.  Before, permuteGeneral would simply return all permutations of a multiset. Now you can specify the lengths of the output.
    }
}

\section{Changes in RcppAlgos version 0.2.1 (Release date: 2017-11-29)}{
    \itemize{
        \item Fixed bug that would occur in two edge cases involving the constraint functions.
        \item Slightly modified formatting for primeSieve.Rd
    }
}

\section{Changes in RcppAlgos version 0.2.0 (Release date: 2017-11-28)}{
    \itemize{
        \item Updated combination algorithms. They are now more than twice as fast.
        \item Updated constraint functions so that memory access is always within container bounds
        \item Consolidated redundant code for outputting different Rcpp types (e.g. IntegerMatrix, CharacterMatrix, etc.) via a templated approach.
        \item Added the function permuteGeneral that is analogous to comboGeneral only instead of combinations, it gives all permutations. It has an additional argument (i.e. 'freqs') that is used to generate permutations of multisets.
        \item All combinatoric functions now support factor types.
    }
}

\section{Changes in RcppAlgos version 0.1.2 (Release date: 2017-11-03)}{
    \itemize{
        \item Corrected minor typo in README file.
        \item Fixed minor error regarding explicitly comparing variables to large numbers that are typed out. Simply adding a decimal along with a zero remedies the situation.
    }
}

\section{Changes in RcppAlgos version 0.1.1 (Release date: 2017-11-03)}{
    \itemize{
        \item Improved ComboConstraint function by removing unnecessary subsetting.
        \item Improved PrimeSieve internal C++ algorithm.
        \item Corrected the errors with respect to the math functions in C++. Explicitly overloaded the parameters of these functions by casting them to the double type.
    }
}

\section{RcppAlgos version 0.1.0 (Release date: 2017-10-26)}{
    \itemize{
        \item Initial Release
    }
}


