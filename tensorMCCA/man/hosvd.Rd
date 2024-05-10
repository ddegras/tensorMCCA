\name{hosvd}
\alias{hosvd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Higher-Order Singular Value Decomposition
}
\description{
Full or truncated HOSVD
}

\usage{
hosvd(x, r = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Tensor (array). 
}
  \item{r}{
Multilinear rank. Can be specified as an integer vector or a single number. 
}
}

\details{
If \code{r} is not specified, a full-rank HOSVD (aka Tucker decomposition) is performed. In other cases, \code{r} is recycled to match the length of \code{dim(x)} as necessary and the (rank) components of this vector are adjusted to not exceed the corresponding dimension in \code{x}. 

}

\value{
A LIST with components
\item{u}{List of factors/singular vectors, one matrix for each mode.}
\item{core}{Core tensor.}
}



\seealso{
\code{\link{tnsr.mat.prod}}
}

\examples{
set.seed(2048)
x <- array(runif(24), 2:4)

## Full-rank HOSVD 
fullsvd <- hosvd(x)  

## Check that full-rank HOSVD exactly decomposes the tensor 
xx <- tnsr.mat.prod(fullsvd$core, fullsvd$factors)
identical(dim(x), dim(xx)) # TRUE
isTRUE(all.equal(x, xx)) # TRUE
range(x - xx) 

## Check orthogonality of factors
lapply(fullsvd$factors, crossprod)

## Truncated HOSVD 
hosvd(x, 1) # HOSVD of multilinear rank (1,1,1)
hosvd(x, c(1,2,2)) # HOSVD of multilinear rank (1,2,2)
}

