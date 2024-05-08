\name{tnsr.mat.prod}
\alias{tnsr.mat.prod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tensor-Matrix Product
}
\description{
Multiply a tensor by one or several matrices in given modes/dimensions. 
}
\usage{
tnsr.mat.prod(x, mat, modes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Tensor (array object). 
}
  \item{mat}{
Single matrix or list of matrices to multiply the tensor with.
}
\item{modes}{
Integer vector of modes in which to perform the tensor-matrix multiplications. 
}
}
\details{
The argument \code{mat} can be provided as a single matrix or a list containing one or several matrices. 

The argument \code{modes} should be of same length as \code{mat}. Or, if \code{mat} has the same length as the number of dimensions of \code{x} (tensor order), the \code{NULL} default will be changed to \code{1:length(mat)}.  


The numbers of columns of the matrices in \code{mat} should match the corresponding dimensions of \code{x}. That is, \code{dim(x)[modes[1]] == ncol(mat[[1]])}, \code{dim(x)[modes[2]] == ncol(mat[[2]])}, ... should all evaluate to \code{TRUE}. 
After the tensor-matrix multiplication, the numbers of rows of the matrices are equal tto the corresponding dimensions of the new tensor.  



}
\value{
Array object.
}

\seealso{
\code{\link{hosvd}}
}

\examples{
x <- array(runif(24), dim = 2:4)
mat <- list(matrix(runif(8),4,2), matrix(runif(8),2,4))
modes <- c(1,3)

## Check that numbers of matrix columns match tensor
## dimensions in associated modes before multiplication
ncolmat <- sapply(mat, ncol) # 2 4
all(dim(x)[modes] == ncolmat) # TRUE

## Perform tensor-matrix multiplication
tmprod <- tnsr.mat.prod(x, mat, modes)

dim(tmprod) # 4 3 2

## Check that numbers of matrix rows match tensor
## dimensions in associated modes after multiplication
nrowmat <- sapply(mat, nrow) # 4 2
all(dim(tmprod)[modes] == nrowmat) # TRUE
}

