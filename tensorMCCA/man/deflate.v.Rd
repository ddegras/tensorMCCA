\name{deflate.v}
\alias{deflate.v}

\title{Deflation of Rank-1 Tensors}

\description{
Deflate (i.e. project orthogonally) rank-1 tensors with respect to other rank-1 tensors. The deflation can apply to the first mode or to all modes of the tensor.  
}

\usage{
deflate.v(v, vprev, ortho = c("canon.t.1", "canon.t.all"))
}

\arguments{
  \item{v}{
List of rank-1 tensors to deflate. 
}
\item{vprev}{
List of rank-1 tensors to deflate against. 
}
\item{ortho}{
Type of deflation: orthogonal projection in first or in all modes. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of the same length as \code{v} 
and with the same component dimensions. 
}




\seealso{
\code{\link{deflate.x}}
}
\examples{
## Tensors to deflate: 2 rank-1 tensors v1 and v2
## with respective dimensions 2 x 5 and 2 x 3 x 5
v <- list(
	v1 = list(mode1 = runif(3), mode2 = runif(5)), 
	v2 = list(mode1 = runif(3), mode2 = runif(4), mode3 = runif(5)))
	
## Tensors to deflate against: 2 rank-1 tensors for each of v1 and v2
vprev <- list(
	v1 = list(mode1 = matrix(runif(6), 3, 2), mode2 = matrix(runif(10), 5, 2)),
 	v2 = list(mode1 = matrix(runif(6), 3, 2), mode2 = matrix(runif(8), 4, 2), 
 		mode3 = matrix(runif(10), 5, 2)))

# Note: list component names are given here for clarity but are not necessary 

## Deflate v1 and v2 in first mode only
vnew <- deflate.v(v, vprev, "canon.t.1")

## Check that orthogonality only applies in first mode
mapply(crossprod, vprev[[1]], vnew[[1]], SIMPLIFY = FALSE)
mapply(crossprod, vprev[[2]], vnew[[2]], SIMPLIFY = FALSE)

## Deflate v1 and v2 in all modes and check orthogonality
vnew <- deflate.v(v, vprev, "canon.t.all")

## Check that orthogonality only applies in first mode
mapply(crossprod, vprev[[1]], vnew[[1]], SIMPLIFY = FALSE)
mapply(crossprod, vprev[[2]], vnew[[2]], SIMPLIFY = FALSE)
}

