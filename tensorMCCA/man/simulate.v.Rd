\name{simulate.v}
\alias{simulate.v}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Random Sets of Rank-1 Tensors
}
\description{
This function randomly generates sets of rank-1 tensors. The tensors of a given set have the same dimensions, but tensor dimensions may vary across sets. Several options are available to control the scale (Euclidean norm) of tensors and the orthogonality constraints between them in each set. 
}

\usage{
simulate.v(p, r, scale = 1, 
	ortho = c("cyclical", "random", "alldim", "maxdim", "none"))
}

\arguments{
  \item{p}{
List of tensor dimensions as integer vectors. 
}
  \item{r}{
Number of tensors to generate for each set of dimensions. 
}
  \item{scale}{
Euclidean norm of tensor for each dataset. Numerical vector of same length as \code{p}. 
}
  \item{ortho}{
Modes in which tensors of the same set are orthogonal to each other. Can be specified as a string or as an integer vector of same length as \code{p}.
}
}

\details{
Let's call \eqn{m} the number of tensor sets (=length of \code{p}), 
\eqn{d_i}{di} the order of the tensors in the \eqn{i}-th set (=length of \code{p[[i]]}), and \eqn{p_{i1},...,p_{id_i}}{p(i,1),...,p(i,di)} the tensor dimensions for the \eqn{i}th set (= entries of \code{p[[i]]}). 
Each generated rank-1 tensor \eqn{v_i^{(l)}}{v(i,l)}, \eqn{i=1,...m}, \eqn{l=1,...,r}, is represented as a list of vectors \eqn{v_{i1}^{(l)}, ..., v_{id_i}^{(l)}}{v(i,1,l), ..., v(i,di,l)} of respective lengths \eqn{p_{i1}, ..., p_{id_i}}{p(i,1),...,p(i,di)}. The actual \eqn{v_i^{(l)}}{v(i,l)} is the outer products of these vectors. 
	
For \eqn{i = 1, ..., m}, the rank-1 tensors \eqn{v_i^{(1)}}{v(i,1)}, ..., \eqn{v_i^{(r)}}{v(i,r)} have their Euclidean norms equal to \code{scale[i]} and they 
are orthogonal to each other in one or several modes \eqn{k} in \eqn{\{1, ..., d_i\}}{\{1, ..., di\}}. Recall that orthogonality between rank-1 tensors \eqn{v_i^{(l_1)}}{v(i,l1)} and \eqn{v_i^{(l_2)}}{v(i,l2)} in mode \eqn{k} means that the associated vectors \eqn{v_{ik}^{(l_1)}}{v(i,k,l1)} and \eqn{v_{ik}^{(l_2)}}{v(i,k,l2)} are orthogonal. Now the modes \eqn{k = k(i,l_1,l_2)}{k = k(i,l1,l2)} (\eqn{l_1 < l_2}{l1 < l2}) in which orthogonality holds are determined by \code{ortho} as follows: 

- \code{cyclical}: \eqn{k = k(i,l_1,l_2)}{k = k(i,l1,l2)} runs cyclically through \eqn{1, ..., d_i, 1, ..., d_i, ...}{1, ..., di, 1, ..., di, ...} as \eqn{(l_1,l_2)}{(l1,l2)} runs through \eqn{(1,2), ..., (1,r)}, \eqn{(2,3), ..., (2,r)},  \eqn{..., (r-1,r)}. 

- \code{random}: \eqn{k=k(i,l_1,l_2)}{k=k(i,l1,l2)} is randomly generated in \eqn{\{1, ..., d_i\}}{\{1, ..., d[i]\}}.  

- \code{alldim}: \eqn{v_{ik}^{(l_1)}}{v(i,k,l1)} and \eqn{v_{ik}^{(l_2)}}{v(i,k,l2)} are orthogonal for all \eqn{k=1, ..., d_i}{k=1, ..., di}.  

- \code{maxdim}: \eqn{k=k(i,l_1,l_2)}{k=k(i,l1,l2)} is the mode associated with the largest dimension \eqn{p_{ik}}{pik} among \eqn{p_{i1},...,p_{id_i}}{pi1,...,pidi}. 

- \code{none}: \eqn{v_i^{(1)}}{v(i,1)}, ..., \eqn{v_i^{(r)}}{v(i,r)} are generated without orthogonality constraints. 
  
}


\value{
A list \code{v} of dimensions \eqn{m \times r}. 
Each component \code{v[[i,l]]} is a list of vectors of respective lengths 
\code{p[[i]][1]}, ..., and \code{p[[i]][length(p[[i]])]}, 
representing a rank-1 tensor.  
}

\references{
%% ~put references to the literature/web site here ~
}


\seealso{
\code{\link{simulate.factor.model}}
}

\examples{
p <- list(10, c(5,10), c(6,8,10)) # data dimensions
r <- 4 # number of factors

## Example 1: simulate rank-1 tensors with unit norm 
## and orthogonal in all dimensions
v <- simulate.v(p, r, ortho = "alldim")
tnsr.rk1.nrm(v) # check norms
tnsr.rk1.cp(v) # check inner products

## Example 2: simulate rank-1 tensors with specified norms 
## and orthogonal in the first dimension
v <- simulate.v(p, r, scale = c(1, 1.25, 1.5), ortho = 1)
tnsr.rk1.nrm(v) # check norms
tnsr.rk1.cp(v) # check inner products (full tensor)
v1 <- matrix(lapply(v, "[", 1), nrow(v), ncol(v)) # extract vectors from 1st dimension
tnsr.rk1.cp(v1) # check that scalar products for these vectors are zero 
}


