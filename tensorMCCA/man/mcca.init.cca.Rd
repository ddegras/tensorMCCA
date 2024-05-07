\name{mcca.init.cca}
\alias{mcca.init.cca}

\title{Initialize Multiple CCA With Standard CCA}

\description{
This function calculates a starting point for multiple CCA by: 
(i) matricizing the data (ignore the tensor structure), 
(ii) calculating candidate canonical vectors for each pair of datasets by standard CCA, 
(iii) recovering best rank-1 tensor approximations to each candidate, 
and (iv) solving exactly or approximately a multidimensional assignment problem to find the best candidate for each dataset. 
}

\usage{
mcca.init.cca(x, w = NULL, objective = c("cov", "cor"), scope = c("block", 
    "global"), k = NULL, optim = NULL, maxit = 100, tol = 1e-04) 
}

\arguments{
\item{x}{list of data arrays. All arrays must have the same number of dimensions in their last mode/dimension.}
\item{w}{matrix of weights for pairs of datasets in objective function.}
\item{objective}{character: maximize the sum of canonical covariances (\code{cov}) or canonical correlations (\code{cor})?
}
\item{scope}{Apply scaling constraints to canonical weights separately for each dataset (\code{block}) or across datasets (\code{global})?}

\item{k}{Rank of truncated SVD to approximate each dataset. Default = \code{NULL} (full SVD).}
\item{optim}{
Type of search in a (possibly large) space of candidate solutions: \code{exhaustive} or \code{approximate}.}
\item{maxit}{Maximum number of iterations in rank-1 tensor approximation.}
\item{tol}{Numerical tolerance in rank-1 tensor approximation.}
}

\details{
It can be useful to set \code{k} to a small-to-moderate value if one or several datasets are very big.

By default, \code{optim} is set to \code{exhaustive} if \eqn{m \le 5} and to 
\code{approximate} otherwise, where \eqn{m} is the length of \code{x}, i.e. the number of datasets.
}


\value{
A LIST of rank-1 canonical tensors: one per dataset. Each rank-1 tensor is specified as a list of canonical vectors. 
}

\references{
Degras, Shen, and Chu (2022). Tensor Multiple Canonical Correlation Analysis.
}



\seealso{
\code{\link{mcca.init.random}}, \code{\link{mcca.init.svd}}
}

\examples{
}

