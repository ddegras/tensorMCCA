\name{mcca.init.svd}
\alias{mcca.init.svd}

\title{Initialize Canonical Tensors by SVD (Quick Rank 1)}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
mcca.init.svd(x, objective = c("cov", "cor"), cnstr = c("block", "global"), 
	center = TRUE, balance = TRUE)
}

\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
\item{objective}{Maximize the sum of covariances or sum of correlations?}
\item{cnstr}{Scaling of canonical tensors: block constraints (i.e. for each dataset) or global constraints?}
\item{center}{Center the data?}
\item{balance}{Should the vectors making up the rank-1 tensors have equal norm?}
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

