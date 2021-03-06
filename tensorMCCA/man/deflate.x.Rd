\name{deflate.x}
\alias{deflate.x}

\title{
Deflate Data Tensors
}

\description{
Deflate data tensors with respect to canonical vectors or canonical scores.
}

\usage{
deflate.x <- function(x, v = NULL, score = NULL, check.args = TRUE)

}

\arguments{
\item{x}{
List of data tensors (matrices or arrays).}

\item{v}{Optional list of sublists of canonical vectors.   }

\item{score}{Optional vector or matrix of canonical scores.}

\item{check.args}{
Should the structure, dimensions, and values of \code{x} be checked? 
}
}

\details{
The arrays \code{x[[1]], x[[2]], ...} may have different orders 
(order = number of array indices) but their last dimension must be the same, 
say \deqn{n}. 

The deflation consists in applying matrix-tensor multiplications 
to \code{x} along certain modes (i.e. dimensions), with the matrices 
equal to orthogonal projectors on vector spaces defined by \code{v} or \code{score}. 
More precisely, if \code{v} is nonnull, for each \deqn{(i,k)} such that  
\code{v[[i]][[k]]} is not empty, the data array \code{x[[i]]} will be multiplied 
 in its mode/dimension \deqn{k} by the projector matrix onto the orthogonal space 
 of \code{v[[i]][[k]]}. 
Similarly, if \code{score} is nonnull, each data array \code{x[[i]]}
will be multiplied in its last mode by the projector onto the orthogonal
space of \code{score[,i]}. (If \code{score} is a vector, 
it will be replicated to form a matrix with \deqn{m} columns, 
where \deqn{m} is the number of data arrays.)


*Only one* of the arguments \code{v} and \code{score} must be specified. 
The other one should be left to its default value \code{NULL}.
 
If \code{v} is specified, this list must have the same length as \code{x}. 
Each sublist \code{v[[i]]} must have its length equal to the order 
(i.e. number of dimensions) of \code{x[[i]]} minus 1. 
Each component \code{v[[i]][[k]]} can either be an empty list 
(no deflation for the \deqn{i}th dataset in the \deqn{k}th mode), 
a vector with length equal to the number \deqn{p_{ik}} of dimensions 
of \code{x[[i]]} in its \code{k}-th mode, or a matrix with \deqn{p_{ik}} rows. 

If \code{score} is specified, it must be a vector of length \deqn{n} 
(to be recycled) or a matrix with \deqn{n} rows and \deqn{m} columns. 



}

\value{
A list of deflated data arrays of the same size as \code{x}.
}



\examples{
## Create list of random arrays
x <- list(matrix(runif(8), 2, 4), array(runif(24), dim = c(2,3,4)), 
	array(runif(48), dim = c(2,3,2,4)))
# The 3 components of x are compatible because their last dimensions are equal (to 4).

## Create compatible list of random (canonical) vectors
v <- vector("list", 3)
v[[1]] <- list(runif(2)) 
# To deflate x[[1]] w.r.t. single vector of length nrow(x[[1]]) = 2
v[[2]] <- list(list(), list(runif(2)), list(list() ) )

## Deflate x with respect to v
x1 <- deflate(x, v = v)


## Deflate x with respect to score
x2 <- deflate(x, score = score)
 



}
