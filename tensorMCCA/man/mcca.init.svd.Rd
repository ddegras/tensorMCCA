\name{mcca.init.svd}
\alias{mcca.init.svd}

\title{Initialize TMCCA Optimization by SVD}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
mcca.init.svd(x, w = NULL, objective = c("cov", "cor"), scope = c("block", 
    "global"), k = NULL, maxit = 100, tol = 1e-04) 
}

\arguments{
\item{x}{List of data arrays. All arrays must have the same value for their last dimension (= individual or object).}
\item{w}{Matrix of weights for pairs of datasets in objective function.}
\item{objective}{Maximize the sum of canonical covariances or correlations?}
\item{scope}{Apply scaling constraints to canonical weights separately for each dataset (\code{block}) or across datasets (\code{global})?}
\item{k}{Rank of truncated SVD to approximate each dataset. Default = \code{NULL} (full SVD).}
\item{maxit}{Maximum number of iterations in rank-1 tensor approximation.}
\item{tol}{Numerical tolerance in rank-1 tensor approximation.}

}

\details{
This initialization consists in first flattening the data tensors along their common mode (= last mode), then concatenating them in one large matrix, and finally obtain canonical tensors from successive rank-1 SVDs of the large matrix.}

\value{
A list of \eqn{m} sublists of canonical vectors, where \eqn{m} is the length of \code{x}. The number of vectors in each sublist equals the number of modes in the corresponding data array minus 1 (not counting the last mode which corresponds to individuals or objects). The lengths of these vectors are equal to the numbers of dimensions in the modes of the corresponding data array.}

\seealso{
\code{\link{mcca.init.cca}}, \code{\link{mcca.init.random}}
}

\examples{
}

