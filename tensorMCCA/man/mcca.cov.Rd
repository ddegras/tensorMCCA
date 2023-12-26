\name{mcca.cov}
\alias{mcca.cov}

\title{Multiple CCA: Maximize Sum of Covariances}

\description{Perform multiple CCA by maximizing sum of covariances between canonical scores.}

\usage{
mcca.cov(x, r = 1, w = NULL, scale = c("block", "global"), 
	ortho = c("score", "weight"), optim = c("bca", "grad.scale", 
	"grad.rotate"), init = c("cca", "svd", "random"), maxit = 1000, 
	tol = 1e-6, sweep = c("cyclical", "random"), control = list(), 
	verbose = FALSE)
}

\arguments{
\item{x}{list of data arrays. All arrays must have the same size in their last dimension (= common instances, cases, or subjects).}
\item{r}{number of sets of canonical components to calculate.}
\item{w}{matrix of weights in objective function. Can be specified as a single number or as a square matrix with numbers of rows and columns equal to the length of \code{x}.}
\item{scale}{character: should the norm constraints on canonical tensors apply to each individual dataset (\code{block}) or on average across datasets (\code{global})?}
\item{ortho}{character: type of orthogonality constraints for higher order canonical components.}
\item{optim}{character: optimization method.}
\item{init}{character (initialization method) or list of initial values for the canonical tensors.}
\item{maxit}{maximum number of algorithm iterations.}
\item{tol}{numerical tolerance for algorithm convergence.}
\item{sweep}{character: sweeping method in block optimization.}
\item{verbose}{display algorithm progress?}
}


\details{
The objective weights \code{w} can be specified as a single positive number (to be recycled) or as an \eqn{m}-by-\eqn{m} matrix, where \eqn{m} is the length of \code{x} (number of datasets). These weights must be nonnegative and have a positive sum. By default (\code{w=NULL}), \code{w} is set to a mmatrix with diagonal coefficients zero and off-diagonal coefficients equal to a positive number. Within the function \code{mcca.cov}, \code{w} is made symmetric and scaled so its elements add up to 1.

For each data array, the number \code{r} of canonical components cannot exceed the number of cases, say \eqn{n}, or the number of elements in each image (the dataset containing \eqn{n} images). If needed, \code{r} will be reduced to meet these limits. 

The scaling constraint \code{scale} can take the value \code{block} (each canonical tensor must have norm equal to 1) or \code{global} (on average across datasets, the canonical tensors have squared norm equal to 1). (at each stage \eqn{1,...,r} of the MCCA). 

The orthogonality constraints \code{ortho} for estimating higher order canonical components (\eqn{2,...,r}) can apply either to canonical scores (\code{score}, defined as \eqn{m} data block-specific scores for each case if \code{scale=block} or as a single global score if \code{scale=global}) or to canonical tensors (\code{weight}). Either way, the constraints are enforced by deflating the data arrays/tensors.  

The optimization method \code{optim} can be set to \code{bca}: block coordinate ascent, \code{grade.scale}: gradient-based with retraction by scaling, or 	\code{grade.rotate}: gradient-based with retraction by rotation.
	
The initialization method \code{init} can be set to \code{cca}: standard CCA between pairs of matricized datasets, \code{svd}: higher order SVD, or \code{random}: random starting point. Alternatively, \code{init} can be a list of length \eqn{m} (the number of datasets), with each component a sublist of vectors whose lengths match the dimensions of the corresponding dataset. For example, if \code{x[[1]]} is an array of dimensions 2-by-3-by-4-by-5, init[[1]] can be specified as a list of 3 numerical vectors of respective lengths 2, 3, and 4. The argument \code{init} can also be a matrix of lists, with the rows indicating datasets \eqn{i=1,...,m} and the columns indicating the optimization stage \eqn{l=1,...,r}. Refer to Degras et al. (2022) for details.

The optional argument \code{control} is a list that provides additional arguments for the optimization. This list can contain any of the named sublists \code{init}, \code{ortho}. The sublist \code{control$init} passes additional arguments to the selected initialization function. (At the moment, only \code{mcca.init.cca} can receive additional arguments; see the help page for this function.)
}

\value{
A LIST with components
\describe{
\item{\code{v}}{canonical tensors. Matrix of lists with \eqn{m} rows and \eqn{r} columns (vector of lists if \code{r=1}). \code{v[[i,l]]} is a list that contains the vectors whose outer product form the canonical tensor solution for the \eqn{i}th dataset at the \eqn{l}th stage of optimization.}
\item{\code{block.score}}{canonical scores: array of dimensions \eqn{n\times m\times r} (individuals x datasets x canonical components).}
\item{\code{global.score}}{Matrix or vector of average  canonical scores across datasets. Dimensions \eqn{n \times r} (individuals x canonical components).}
\item{\code{objective}}{objective value (weighted sum of covariances) for each canonical component (vector of length \code{r})}
\item{\code{iters}}{Number of algorithm iterations for each canonical component (length \code{r})}
\item{\code{input}}{List of input arguments (all but \code{x}).}
}
}



\seealso{
\code{\link{mcca.cor}}, \code{\link{mcca.init.cca}}, \code{\link{mcca.init.svd}}, \code{\link{mcca.init.random}}
}

\examples{
## Simulate three tensor datasets measured on 100 cases with 4 common factors
m <- 3 
n <- 100 
r <- 4 
dimx <- list(c(10,n), c(5,10,n), c(6,8,10,n)) 
sim <- simulate.factor.model(dimx, r, score.cov = r:1, noise.cov = .5, ortho.v = 1)
x <- sim$x

fit <- mcca.cov(x, r, verbose = TRUE)
}
