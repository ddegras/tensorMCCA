\name{canon.scores}
\alias{canon.scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Canonical Scores
}
\description{
Calculate the canonical scores of observed individuals/objects for given datasets and canonical tensors. 
}
\usage{
canon.scores(x, v)
}

\arguments{
  \item{x}{list of data arrays.}
  \item{v}{list of lists. Each sublist is a list of single vectors that specify a canonical tensor.}
}

\details{
The list \code{x} contains data arrays such that the last mode/dimension of each array indexes the same individuals or objects. Calling \eqn{m} the number of datasets (i.e., the length of \code{x}) and \eqn{n} the common number of individuals/objects, each data array \code{x[[i]]} has dimensions of the form 
\eqn{(p[i1], ..., p[id[i]], n)} for \eqn{i=1,...,m}.

The list \code{v} contains sublists of vectors such that each sublist specifies rank-1 canonical tensor. The list \code{v} can either have length \eqn{m} (no dimension attributes, single set of canonical tensors) or have dimensions \eqn{m}-by-\eqn{r}, where \eqn{r} is the number of sets of canonical tensors to consider. 
The length of each vector in \code{v} must match the corresponding dimension in \code{x}. More precisely, \code{v[[i]]} or \code{v[i,l]} (\eqn{i = 1,...,m, l = 1,...,r}) must be a list of \eqn{d[i]} vectors of respective lengths \eqn{p[i1], ..., p[id[i]]}. 
 
For each dataset (\eqn{i=1,...,m}), individual (\eqn{t=1,...,n}), and canonical component (\eqn{l=1,...,r}), the function \code{canon.score} calculates the inner (tensor) product \eqn{<X[it],v[il]>} and stores this score in a matrix (if \eqn{r=1}) or in an array (if \eqn{r>1}). 
}

\value{
If \eqn{r = 1} (\code{NCOL(v) == 1}), a matrix of size \eqn{n}-by-\eqn{m} with individuals in rows and datasets in columns. \cr
If \eqn{r > 1}, (\code{NCOL(v) > 1}), an array of dimensions \eqn{(n, m, r)}. 
}

\examples{
## Simulate data
x <- list(matrix(runif(50), 5, 10), array(runif(60), c(2, 3, 10)), 
	array(runif(80), c(2, 2, 2, 10))) 
	
## Example with single set of canonical vectors
v1 <- list(list(runif(5)), list(runif(2), runif(3)), 
	list(runif(2), runif(2), runif(2))) 
stopifnot(length(x) == length(v1))
canon.scores(x, v1) # 10-by-3 matrix

## Example with multiple sets of canonical vectors
m <- 3 # number of datasets
r <- 2 # number of sets of canonical vectors
p <- list(list(5), list(2, 3), list(2, 2, 2)) # data dimensions
v2 <- vector("list", m * r)
dim(v2) <- c(m, r)
for (i in 1:m) {
	for (l in 1:r) {
	v2[[i,l]] <- lapply(p[[i]], runif)
	}
}
canon.scores(x, v2) # 10-by-3-by-2 array
}

