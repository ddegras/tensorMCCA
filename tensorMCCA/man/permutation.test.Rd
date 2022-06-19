\name{permutation.test}
\alias{permutation.test}

\title{
Permutation Testing for Canonical Components 
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
permutation.test(x, obj, nperm = 100, parallel = TRUE, ncores = NULL)
}

\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}

\item{obj}{The result of a call to \code{\link{MCCA.cor}} or \code{\link{MCCA.cov}}.}

\item{nperm}{Number of random permutations to use for the test.}

\item{parallel}{Should calculations be conducted in parallel?}

\item{ncores}{Number of cores for parallel computation.}
}

\details{
For large datasets, it may be useful to first carry out the permutation test with a small number of permutations. This may help calibrate the number of canonical components to test with a larger number of permutations - and reduce computation time. 
}

\value{
A vector of p-values, one for each canonical component.
}


\seealso{
\code{\link{MCCA.cov}}, \code{\link{MCCA.cor}}
}

\examples{
}
