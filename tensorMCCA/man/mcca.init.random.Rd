\name{mcca.init.random}
\alias{mcca.init.random}

\title{
Initialize Canonical Tensors Randomly
}

\description{
Generate random rank-1 tensors of given size under scaling constraints and optional orthogonality constraints.}

\usage{
mcca.init.random(x, r = 1L, objective = c("cov", "cor"), ortho.mode = NULL)
}

\arguments{
\item{x}{list of data arrays. All arrays must have the same last dimension (= individual or object).}
\item{r}{integer number of sets of canonical tensors to calculate.}
\item{objective}{character indicating the objective function: sum of covariances (\code{cov}) or sum of correlations (\code{cor}).}
\item{ortho.mode}{integer or list indicating in which mode(s) -if any- to orthogonalize the canonical tensors of each dataset. This argument is recycled as needed to a list of same length as \code{x}.}

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
A vector or matrix of lists representing canonical tensors. If \eqn{r = 1}, the result is a list of length \eqn{m}, the number of datasets in \code{x}. If \eqn{r > 1}, the result is an \eqn{m} by \eqn{r} matrix. Each list component is itself a list of canonical vectors (as many as the are modes in the corresponding dataset, not counting the last mode which corresponds to replications). The canonical vectors have length equal to the numbers of dimensions in the modes of the corresponding dataset.}

\seealso{
\code{\link{mcca.init.cca}}, \code{\link{mcca.init.svd}}
}

\examples{
}
