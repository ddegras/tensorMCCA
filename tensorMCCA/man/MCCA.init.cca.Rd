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
mcca.init.cca(x, k = NULL, c = 1, objective = c("cov", "cor"), norm = c("block", "global"), 
center = TRUE, search = c("exhaustive", "approximate"))
}

\arguments{
\item{x}{list of data arrays. All arrays must have the same number of dimensions in their last mode/dimension.}
\item{k}{integer: rank of (truncated) SVD to approximate each dataset. Default = full SVD.}
\item{w}{matrix of weights for pairs of datasets in objective function.}
\item{objective}{character: maximize sum of covariances (\code{cov}) or sum of  
or \code{cor}?
}
\item{norm}{character: should norm constraints for canonical tensors apply at the \code{block} level, that is, separately for each dataset, or at a \code{global} level, that is, on average across datasets?}
\item{center}{
Logical: should the data be centered before analysis? Default = \code{TRUE}.
}
\item{search}{
Type of search in a (possibly large) space of candidate solutions: \code{exhaustive} or \code{approximate}.}
}

\details{
}


\value{
A LIST of canonical tensors with dimensions \eqn{m}-by-\eqn{r} (dataset x component). Each tensor is specified as a list of canonical vectors. 
}

\references{
Degras, Shen, and Chu (2022). Tensor Multiple Canonical Correlation Analysis.
}



\seealso{
\code{\link{mcca.init.random}}, \code{\link{mcca.init.svd}}
}

\examples{
}

