\name{mcca.cov}
\alias{mcca.cov}

\title{Multiple CCA: Maximize Sum of Covariances}

\description{Perform multiple CCA by maximizing sum of covariances between canonical scores.}

\usage{
mcca.cov(x, r, w = 1, cnstr = c("block", "global"), ortho = c("block.score", "global.score", "canon.tnsr"), 
	init = c("cca", "svd", "random"), maxit = 1000, tol = 1e-6, sweep = c("cyclical", "random"), 
	control = list(), verbose = FALSE)
}

\arguments{
\item{x}{list of data arrays. All arrays must have the same size in their last dimension (= common instances, cases, or subjects).}
\item{r}{number of sets of canonical components to calculate.}
\item{w}{matrix of weights in objective function. Can be specified as a single number or as a square matrix with numbers of rows and columns equal to the length of \code{x}.}
\item{cnstr}{character: should the norm constraints on canonical tensors apply to each individual dataset (\code{block}) or on average across datasets (\code{global})?}
\item{ortho}{character: type of orthogonality constraints for higher order canonical components.}
\item{init}{character (initialization method) or list of initial values for the canonical tensors.}
\item{maxit}{maximum number of algorithm iterations.}
\item{tol}{numerical tolerance for algorithm convergence.}
\item{sweep}{character: sweeping method in block optimization.}
\item{verbose}{display algorithm progress?}
}

\details{
The objective weights \code{w} can be specified as a single positive number (to be recycled) or as an \eqn{m}-by-\eqn{m} matrix, where \eqn{m} is the length of \code{x} (number of datasets). These weights must be nonnegative and have a positive sum. Within the function \code{mcca.cov}, \code{w} is made symmetric and scaled so its elements add up to 1.

For each dataset, the number \code{r} of canonical components cannot exceed the number of cases, say \eqn{n}, or the number of elements in each image (the dataset containing \eqn{n} images). If needed, \code{r} will be reduced to meet these limits. 

The orthogonality constraints for estimating higher order canonical components (\eqn{2,...,r}) can apply either to canonical scores (defined as \eqn{m} data block-specific scores or as a single global score for each case) or to canonical tensors. Either way, the constraints are enforced by deflating the data arrays/tensors.  

The norm constraint \code{cnstr = "block"} requires that each canonical tensor has norm equal to 1. The weaker constraint \code{cnstr = "global"} requires that on average across datasets, the canonical tensors have norm equal to 1 (at each stage \eqn{1,...,r} of the MCCA). 
	
The optional list \code{init} can take as arguments: 
\describe{
\item{\code{method}}{initialization method: \code{cca} (standard CCA between pairs of matricized datasets), \code{svd} (higher order SVD), or \code{random}. Refer to Degras et al. (2022) for details.}

\item{\code{value}}{List of starting values for the mcca. If \code{value} is provided, the argument \code{method} will be ignored. \code{value} can either have the same length as \code{x}, say \eqn{m}, or have dimensions \eqn{m}-by-\code{r}. In the first case, the same starting point is used for all \code{r} stages of the mcca. In the second, each column of \code{value} will be used as separate starting point. The \eqn{i}-th field of \code{value} (or all fields of its \eqn{i}-th row) must match the dimensions of the \eqn{i}-th field of \code{x}. More precisely, the \eqn{i}-th field of \code{value} must be a list of length equal to that of \code{x[[i]]} minus 1. If for example \code{x[[i]]} is a 4D array of 3D images of dimensions 10-by-15-by-20, then \code{value[[i]]} must be a list of length 3 with each component a single vector of length 10, 15, and 20, respectively.}

\item{\code{search}}{Only used for \code{method = cca}. 
Type of search in the space of possible solutions (size \eqn{m^m} where \eqn{m} is the number of datasets): \code{exhaustive} or \code{approximate}. By default, \code{search} is set to \code{exhaustive} if \eqn{m \le 5} and to 
\code{approximate} otherwise.} 

\item{\code{k}}{Only used for \code{method = cca}. Number of singular values used to approximate each dataset. It can be useful to set \code{k} to a small-to-moderate value if one or several datasets are very big.}
}
}

\value{
A LIST with components
\describe{
\item{\code{v}}{canonical tensors}
\item{\code{block.score}}{canonical scores: array of dimensions individuals x datasets x canonical components}
\item{\code{global.score}}{weighted average of the canonical scores with weights given by the column sums of \code{c} (individuals x canonical components)}
\item{\code{objective}}{objective value (sum of correlations) for each canonical component (vector of length \code{r})}
\item{\code{iters}}{Number of algorithm iterations for each canonical component (length \code{r})}
\item{\code{input}}{List of input arguments}
}
}



\seealso{
\code{\link{mcca.cor}}
}

% \examples{
% }
