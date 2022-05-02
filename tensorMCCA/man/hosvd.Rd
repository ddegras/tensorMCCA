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
If \code{r} is not specified, a full HOSVD is performed. In other cases, \code{r} is recycled to match the length of \code{dim(x)} as necessary and the (rank) components of this vector are adjusted to not exceed the corresponding dimension in \code{x}. 

}

\value{
A LIST with components
\item{u}{List of factors/singular vectors, one matrix for each mode.}
\item{core}{Core tensor}
}



%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
x <- array(runif(24), 2:4)
hosvd(x, 1)
hosvd(x, c(1,2,2))
}

