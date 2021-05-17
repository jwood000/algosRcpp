\name{comboGeneral}
\alias{comboGeneral}
\alias{permuteGeneral}
\title{
Generate Combinations and Permutations of a Vector with/without Constraints
}
\description{
    \itemize{
        \item Generate combinations or permutations of a vector with or without constraints.
        \item Produce results in parallel using the \code{Parallel} or \code{nThreads} arguments. You can also apply each of the five compiled functions given by the argument \code{constraintFun} in parallel.
        \item Alternatively, the arguments \code{lower} and \code{upper} make it possible to generate combinations/permutations in chunks allowing for parallelization via the parallel package. This is convenient when you want to apply a custom function to the output in parallel.
        \item Attack integer partition and general subset sum problems.
        \item GMP support allows for exploration of combinations/permutations of vectors with many elements.
        \item The output is in lexicographical order.
    }
}
\usage{
comboGeneral(v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL, 
              upper = NULL, constraintFun = NULL, comparisonFun = NULL,
               limitConstraints = NULL, keepResults = NULL, FUN = NULL,
                Parallel = FALSE, nThreads = NULL, tolerance = NULL)

permuteGeneral(v, m = NULL, repetition = FALSE, freqs = NULL, lower = NULL, 
                upper = NULL, constraintFun = NULL, comparisonFun = NULL,
                 limitConstraints = NULL, keepResults = NULL, FUN = NULL,
                  Parallel = FALSE, nThreads = NULL, tolerance = NULL)
}
\arguments{
  \item{v}{Source vector. If \code{v} is an integer (including nonpositive integers), it will be converted to the sequence \code{1:v}. All atomic types are supported (See \code{?is.atomic}).}
  \item{m}{Number of elements to choose.  If \code{repetition = TRUE} or \code{freqs} is utilized, \code{m} can exceed the length of \code{v}. If \code{m = NULL}, the length will default to \code{length(v)} or \code{sum(freqs)}.}
  \item{repetition}{Logical value indicating whether combinations/permutations should be with or without repetition. The default is \code{FALSE}.}
  \item{freqs}{A vector of frequencies used for producing all combinations/permutations of a multiset of \code{v}. Each element of \code{freqs} represents how many times each element of the source vector, \code{v}, is repeated. It is analogous to the \code{times} argument in \code{\link[base]{rep}}. The default value is \code{NULL}.}
  \item{lower}{The lower bound. Combinations/permutations are generated lexicographically, thus utilizing this argument will determine which specific combination/permutation to start generating from (\emph{e.g.} \code{comboGeneral(5, 3, lower = 6)} is equivalent to \code{comboGeneral(5, 3)[6:choose(5, 3), ]}). This argument along with \code{upper} is very useful for generating combinations/permutations in chunks allowing for easy parallelization.}
  \item{upper}{The upper bound. Similar to \code{lower}, however this parameter allows the user to \emph{stop} generation at a specific combination/permutation (\emph{e.g.} \code{comboGeneral(5, 3, upper = 5)} is equivalent to \code{comboGeneral(5, 3)[1:5, ]})
  
  If the output is constrained and \code{lower} isn't supplied, \code{upper} serves as a cap for how many results will be returned that meet the criteria (\emph{e.g.} setting \code{upper = 100} alone will return the first 100 results that meet the criteria, while setting \code{lower = 1} and \code{upper = 100} will test the first 100 results against the criteria).
  
  In addition to the benefits listed for \code{lower}, this parameter is useful when the total number of combinations/permutations without constraint is large and you expect/need a small number of combinations/permutations that meet a certain criteria. Using \code{upper} can improve run time if used judiciously as we call the member function \href{https://en.cppreference.com/w/cpp/container/vector/reserve}{reserve} of \href{https://en.cppreference.com/w/cpp/container/vector}{std::vector}. See examples below.}
  \item{constraintFun}{Function to be applied to the elements of \code{v} that should be passed as a string (\emph{e.g.} \code{constraintFun = "sum"}). The possible constraint functions are: \code{"sum"}, \code{"prod"}, \code{"mean"}, \code{"max"}, & \code{"min"}. The default is \code{NULL}, meaning no function is applied.}
  \item{comparisonFun}{Comparison operator that will be used to compare \code{limitConstraints} with the result of \code{constraintFun} applied to \code{v}. It should be passed as a string or a vector of two strings (\emph{e.g.} \code{comparisonFun = "<="} or  \code{comparisonFun = c(">","<")}). The possible comparison operators are: \code{"<"}, \code{">"}, \code{"<="}, \code{">="}, \code{"=="}. The default is \code{NULL}.
  
  When \code{comparisonFun} is a vector of two comparison strings, \emph{e.g} \code{comparisonFun = c(comp1, comp2)}, and \code{limitConstraints} is a vector of two numerical values, \emph{e.g} \code{limitConstraints = c(x1, x2)}, the combinations/permutations will be filtered in one of the following two ways:
  \enumerate{
      \item{When \code{comp1} is one of the 'greater-than' operators (\emph{i.e.} ">=" or ">"), \code{comp2} is one of the 'less-than' operators (\emph{i.e.} "<=" or "<"), and \code{x1 < x2}, the combinations/permutations that are returned will have a value (after \code{constraintFun} has been applied) between \code{x1} and \code{x2}.}
      \item{When \code{comp1} and \code{comp2} are defined as in #1 and \code{x1 > x2}, the combinations/permutations that are returned will have a value outside the range of \code{x1} and \code{x2}. See the examples below.}
        }
  In other words, the first comparison operator is applied to the first limit and the second operator is applied to the second limit.}
  \item{limitConstraints}{This is the value(s) that will be used for comparison. Can be passed as a single value or a vector of two numerical values. The default is \code{NULL}. See the definition of \code{comparisonFun} as well as the examples below for more information.}
  \item{keepResults}{A logical flag indicating if the result of \code{constraintFun} applied to \code{v} should be displayed; if \code{TRUE}, an additional column of results will be added to the resulting matrix. The default is \code{FALSE}. If user is only applying \code{constraintFun}, \code{keepResults} will default to \code{TRUE}.
  
      \emph{E.g}. The following are equivalent and will produce a \eqn{4^{th}}{4th} column of row sums:
      \itemize{
            \item \code{comboGeneral(5, 3 constraintFun = "sum", keepResults = TRUE)}
            \item \code{comboGeneral(5, 3 constraintFun = "sum")}
        }
    }
  \item{FUN}{Function to be applied to each combination/permutation. The default is \code{NULL}.}
  \item{Parallel}{Logical value indicating whether combinations/permutations should be generated in parallel using \eqn{n - 1} threads, where \eqn{n} is the maximum number of threads. The default is \code{FALSE}. If \code{nThreads} is not \code{NULL}, it will be given preference (\emph{e.g.} if user has 8 threads with \code{Parallel = TRUE} and \code{nThreads = 4}, only 4 threads will be spawned). If your system is single-threaded, the arguments \code{Parallel} and \code{nThreads} are ignored.}
  \item{nThreads}{Specific number of threads to be used. The default is \code{NULL}. See \code{Parallel}.}
  \item{tolerance}{A numeric value greater than or equal to zero. This parameter is utilized when a constraint is applied on a numeric vector. The default value is 0 when it can be determined that whole values are being utilized, otherwise it is \code{sqrt(.Machine$double.eps)} which is approximately \eqn{1.5e-8}. N.B. If the input vector is of type integer, this parameter will be ignored and strict equality will be enforced.}
}
\value{
In general, a matrix is returned with each row containing a vector of length \eqn{m} or \eqn{m + 1} depending on the value of \code{keepResults}. If \code{FUN} is utilized, a list is returned.
}
\details{
Finding all combinations/permutations with constraints is optimized by organizing them in such a way that when \code{constraintFun} is applied, a \emph{partially} monotonic sequence is produced. Combinations/permutations are added successively, until a particular combination exceeds the given constraint value for a given constraint/comparison function combo. After this point, we can safely skip several combinations knowing that they will exceed the given constraint value.

When there are any negative values in \code{v} and \code{constraintFun = "prod"}, producing a monotonic set is non-trivial for the general case. As a result, performance will suffer as all combinations/permutations must be tested against the constraint criteria. Additionally, \code{upper} will not have its normal effectiveness (\emph{i.e.} it will only limit the number of rows after producing all combinations/permutations).
}
\note{
    \itemize{
        \item {\code{Parallel} and \code{nThreads} will be ignored in the following cases:
            \itemize{
                \item When the output is constrained.
                \item If the class of the vector passed is \code{character}, \code{raw}, and \code{complex} (N.B. \code{Rcpp::CharacterMatrix} is not thread safe). Alternatively, you can generate an indexing matrix in parallel.
                \item If \code{FUN} is utilized.
            }
        }
        \item If either \code{constraintFun}, \code{comparisonFun} or \code{limitConstraints} is \code{NULL} --or-- if the class of the vector passed is \code{logical}, \code{character}, \code{raw}, and \code{complex}, the constraint check will not be carried out. This is equivalent to simply finding all combinations/permutations of \eqn{v} choose \eqn{m}.
        \item The maximum number of combinations/permutations that can be generated at one time is \eqn{2^{31} - 1}{2^31 - 1}. Utilizing \code{lower} and \code{upper} makes it possible to generate additional combinations/permutations.
        \item Factor vectors are accepted. Class and level attributes are preserved.
        \item Lexicographical ordering isn't guaranteed for permutations if \code{lower} isn't supplied and the output is constrained.
        \item If \code{lower} is supplied and the output is constrained, the combinations/permutations that will be tested will be in the lexicographical range \code{lower} to \code{upper} or up to the total possible number of results if \code{upper} is not given. See the second paragraph for the definition of \code{upper}.
        \item \code{FUN} will be ignored if the constraint check is satisfied.
    }
}
\author{
Joseph Wood
}
\references{
    \itemize{
        \item{\href{http://gallery.rcpp.org/articles/passing-cpp-function-pointers/}{Passing user-supplied C++ functions}}
        \item{\href{https://en.wikipedia.org/wiki/Monotonic_function}{Monotonic Sequence}}
        \item{\href{https://en.wikipedia.org/wiki/Multiset}{Multiset}}
        \item{\href{https://en.wikipedia.org/wiki/Lexicographical_order}{Lexicographical Order}}
        \item{\href{https://en.wikipedia.org/wiki/Subset_sum_problem}{Subset Sum Problem}}
        \item{\href{https://en.wikipedia.org/wiki/Partition_(number_theory)}{Partition (Number Theory)}}
    }
}
\examples{
system.time(comboGeneral(17, 8))
system.time(permuteGeneral(13, 6))

system.time(comboGeneral(13,10,repetition = TRUE))
system.time(permuteGeneral(factor(letters[1:9]),6,TRUE))

## permutations of the multiset (with or w/o setting m) :
## c(1,1,1,1,2,2,2,3,3,4)
system.time(permuteGeneral(4, freqs = c(4,3,2,1)))

#### Examples using "upper" and "lower":
## Generate some random data
set.seed(1009)
mySamp = sort(rnorm(75, 997, 23))

permuteCount(75, 10, freqs = rep(1:3, 25))
# >Big Integer ('bigz') :
# >[1] 4606842576291495952

## See specific range of permutations
permuteGeneral(75, 10, freqs = rep(1:3, 25), 
               lower = 1e12, upper = 1e12 + 10)

## Researcher only needs 1000 7-tuples of mySamp
## such that the sum is greater than 7200.
system.time(comboGeneral(mySamp, 7, FALSE, constraintFun = "sum",
             comparisonFun = ">", limitConstraints = 7200, upper = 1000))

## Similarly, you can use "lower" to obtain the last rows.
## Generate the last 10 rows
system.time(comboGeneral(mySamp, 7, lower = choose(75, 7) - 9))

## Or if you would like to generate a specific chunk,
## use both "lower" and "upper". E.g. Generate one
## million combinations starting with the 900,000,001
## lexicographic combination.
t1 = comboGeneral(mySamp, 7, lower = 9*10^8 + 1, 
                        upper = 9*10^8 + 10^6)

## class of the source vector is preserved
class(comboGeneral(5,3)[1,]) == class(1:5)
class(comboGeneral(c(1,2:5),3)[1,]) == class(c(1,2:5))
class(comboGeneral(factor(month.name),3)[1,]) == class(factor(month.name))

## Using keepResults will add a column of results
t2 = permuteGeneral(-3,6,TRUE,constraintFun = "prod", 
                         keepResults = TRUE)

t3 = comboGeneral(-3,6,TRUE,constraintFun = "sum", 
                        comparisonFun = "==",
                        limitConstraints = -8, 
                        keepResults = TRUE)
                        
## Using multiple constraints:

## Get combinations such that the product
## is between 3000 and 4000 inclusive
comboGeneral(5, 7, TRUE, constraintFun = "prod",
             comparisonFun = c(">=","<="),
             limitConstraints = c(3000, 4000),
             keepResults = TRUE)
             
## Or, get the combinations such that the
## product is less than or equal to 10 or
## greater than or equal to 40000
comboGeneral(5, 7, TRUE, constraintFun = "prod",
             comparisonFun = c("<=",">="),
             limitConstraints = c(10, 40000),
             keepResults = TRUE)
             
#### General subset sum problem
set.seed(516781810)
comboGeneral(runif(100, 0, 42), 5, constraintFun = "mean",
             comparisonFun = "==", limitConstraints = 30,
             tolerance = 0.0000002)

             
#### Integer Partitions
comboGeneral(0:5, 5, TRUE, constraintFun = "sum",
             comparisonFun = "==", limitConstraints = 5)
             
             
## Using FUN
comboGeneral(10000, 5, lower = 20, upper = 22, 
             FUN = function(x) {
                 which(cummax(x) \%\% 2 == 1)
             })

\dontrun{
## Parallel example generating more than 2^31 - 1 combinations.
library(parallel)
numCores = detectCores() - 1

system.time(mclapply(seq(1, comboCount(35, 15), 10086780), function(x) {
    a = comboGeneral(35, 15, lower = x, upper = x + 10086779)
    ## do something
    x
}, mc.cores = numCores))


## Find 13-tuple combinations of 1:25 such
## that the mean is less than 10
system.time(myComb <- comboGeneral(25, 13, FALSE,
                                   constraintFun = "mean",
                                   comparisonFun = "<",
                                   limitConstraints = 10))

## Alternatively, you must generate all combinations and subsequently
## subset to obtain the combinations that meet the criteria
system.time(myComb2 <- combn(25, 13))
system.time(myCols <- which(colMeans(myComb2) < 10))
system.time(myComb2 <- myComb2[, myCols])

## Any variation is much slower
system.time(myComb2 <- combn(25, 13)[,combn(25, 13, mean) < 10])

## Test equality with myComb above
all.equal(myComb, t(myComb2))


## much faster using Parallel = TRUE
system.time(permuteGeneral(15, 8, lower = 250000000, Parallel = TRUE))
system.time(permuteGeneral(15, 8, lower = 250000000))

system.time(comboGeneral(30, 10, Parallel = TRUE))
system.time(comboGeneral(30, 10))

## Parallel also works when applying constraintFun solely
system.time(comboGeneral(30, 10, Parallel = TRUE, constraintFun = "sum"))
system.time(comboGeneral(30, 10, constraintFun = "sum"))

## Using Parallel argument only uses 1 minus the maximal # of threads.
## The below should run a bit faster as we are utilizing all threads.
maxThreads = RcppAlgos::stdThreadMax()
system.time(comboGeneral(30, 10, nThreads = maxThreads, constraintFun = "sum"))

## Depending on # of cores available, using Parallel with
## constraintFun is faster than rowSums or rowMeans alone
combs = comboGeneral(30, 10)
system.time(rowSums(combs))


## Fun example... see stackoverflow:
## https://stackoverflow.com/q/22218640/4408538
system.time(permuteGeneral(seq(0L,100L,10L),8,TRUE,
                                constraintFun = "sum", 
                                comparisonFun = "==",
                                limitConstraints = 100))
}
}