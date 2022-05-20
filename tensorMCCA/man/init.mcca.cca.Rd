\name{init.mcca.cca}
\alias{init.mcca.cca}

\title{Initialize Multiple CCA With Standard CCA}

\description{
This function calculates a starting point for multiple CCA by: 
(i) matricizing the data (ignore the tensor structure), 
(ii) calculating (long) candidate canonical vectors for each pair of datasets by standard CCA, 
(iii) recovering best rank-1 tensor approximations to each candidate, 
and (iv) solving exactly or approximately a multidimensional assignment problem to find the best candidate for each dataset. 
}

\usage{
init.mcca.cca(x, k = NULL, c = 1, objective = c("covariance", "correlation"), cnstr = c("block", "global"), 
center = TRUE, search = c("exhaustive", "approximate"))
}

\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
\item{k}{Order of (truncated) SVD to approximate each dataset.}
\item{c}{Weight matrix for pairs of datasets.}
\item{objective}{
Objective to maximize in MCCA: \code{covariance} 
or \code{correlation}. 
}
\item{cnstr}{Norm constraints for canonical tensors.}
\item{center}{
Logical: should the data be centered before analysis? Default = \code{TRUE}.
}
\item{search}{
Type of search in a (possibly large) space of candidate solutions: \code{exhaustive} or \code{approximate}.}
}

\details{
}


\value{
A LIST with components
\describe{
\item{\code{v}}{List of canonical tensors with dimensions \eqn{m}-by-\eqn{r} (dataset x component)}
}
}

\references{
Degras, Shen, and Chu (2022). Tensor Multiple Canonical Correlation Analysis.
}



\seealso{
\code{\link{init.mcca.random}}, \code{\link{init.mcca.svd}}
}

\examples{
}

