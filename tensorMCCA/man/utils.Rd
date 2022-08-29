\name{tnsr.rk1.utils}
\alias{tnsr.rk1.cp}
\alias{tnsr.rk1.expand}
\alias{tnsr.rk1.nrm}

\title{
Utility functions for Rank-1 Tensors}

\description{
Utility functions rank-1 tensors specified as lists of vectors. 
Functions include calculating inner products and norms of tensors 
as well as expanding tensors from list format to array.  
}

\usage{
tnsr.rk1.cp(v, v2 = NULL)
tnsr.rk1.expand(v)
tnsr.rk1.nrm(v, norm = c("block", "global"))
}

\arguments{
  \item{v}{
list of lists in vector or matrix format. Each sublist must contain numeric vectors whose outer product form a rank-1 tensor.   
}
\item{v2}{optional list. If specified, must have the same dimensions as \code{v}.}

\item{norm}{
character: type of norm to be calculated. For \code{block} 
}
}
\details{
In function \code{tnsr.rk1.cp}, if \code{v} is a (list) matrix, 
all sublists of a given row must specify tensors of the same dimensions. 
That is, for each row \eqn{i}, the lists \code{x[[i,1]]}, \code{x[[i,2]]}, etc, 
must contain the same number of vectors and these vectors must have the same length. 
If \code{v2} is \code{NULL}) and \code{v} is a (list) vector, \code{tnsr.rk1.cp} returns the squared norms of the tensors specified by \code{v}. 
If \code{v2} is \code{NULL}) and \code{v} is a (list) matrix of dimensions \eqn{m}-by-\eqn{r}, \code{tnsr.rk1.cp} returns an array of dimensions \eqn{r}-by-\eqn{r}-by-\eqn{m} whose entry \eqn{(k,l,i)} is the inner product between tensors \code{v[[i,k]]} and \code{v[[i,l]]}. If \code{v2} is not \code{NULL}, its dimensions must match those of \code{v}. In this case, \code{tnsr.rk1.cp} returns an array of same dimensions as \code{v} containing the inner products between \code{v[[i]]} and \code{v2[[i]]}. 

Function \code{tnsr.rk1.expand} expands each sublist in \code{v} into a rank-1 tensor by calculating the outer product of its elements. For example, if \code{v[[i]]} is a list of 3 numerical vectors of respective lengths 2, 3, and 4, the associated output will be an array of dimensions 2-by-3-by-4. 

Function \code{tnsr.rk1.nrm} calculates the norm of each tensor specified by \code{v} if \code{norm} is \code{block}. If \code{norm} is \code{global}, \code{v} is treated as a matrix and \code{tnsr.rk1.nrm} returns the square root of the average of the squared tensor norms for each row in \code{v}. That is, for each row \eqn{i}, the squared norms of the rank-1 tensors (specified by) \code{v[[i,1]], ..., v[[i,NCOL(v)]]} are calculated, averaged, and finally square-rooted. 

}

\value{
\code{tnsr.rk1.cp} returns a numerical vector or array depending on the dimensions of \code{v}.  

\code{tnsr.rk1.expand} returns a list of same length and dimensions as \code{v}. The elements of that list are arrays whose dimensions match the length of the associated vectors in \code{v}.

\code{tnsr.rk1.norm} returns a numerical array of tensor norms eith same dimensions as \code{v} if \code{norm} equals \code{block}; otherwise it returns a numerical vector of average norms (one for each row in \code{v}). 
}



\examples{
## Generate matrix of lists of random vectors
v <- vector("list", 6)
dim(v) <- c(2, 3)
v[1,] <- replicate(3, list(lapply(4:5, runif)))
v[2,] <- replicate(3, list(lapply(4:6, runif)))
# The elements in the rows of v have same dimensions

## Calculate the inner products between these rank-1 tensors
cp <- tnsr.rk1.cp(v) # 3-by-3-by-2 array

## Compare to direct calculation through tensor expansion
cp2 <- array(dim = c(3, 3, 2))
for (i in 1:2) {
	for (k in 1:3) {
		for (l in 1:k) {
			vik <- tnsr.rk1.expand(v[[i,k]])
			vil <- tnsr.rk1.expand(v[[i,l]])
			cp2[k,l,i] <- sum(vik * vil)
			if (l < k) cp2[l,k,i] <- cp2[k,l,i]
		}
	}
}
stopifnot(isTRUE(all.equal(cp, cp2)))

## Note: all tensor expansions can be calculated in one go as
tnsr <- tnsr.rk1.expand(v) # 2-by-3 matric list of arrays
tnsr[[1,2]]

## Calculate tensor inner products with 2 arguments 
cp <- tnsr.rk1.cp(v[,1], v[,2])
cp2 <- c(sum(tnsr[[1,1]] * tnsr[[1,2]]), 
	sum(tnsr[[2,1]] * tnsr[[2,2]]))
stopifnot(isTRUE(all.equal(cp, cp2)))

## Calculate tensor norms
nrm <- tnsr.rk1.nrm(v, "block")
nrm2 <- matrix(sapply(tnsr, function(x) sqrt(sum(x^2))), 2, 3)
stopifnot(isTRUE(all.equal(nrm, nrm2)))
global.nrm <- tnsr.rk1.nrm(v, "global")
stopifnot(isTRUE(all.equal(global.nrm, sqrt(rowMeans(nrm^2)))))
}

