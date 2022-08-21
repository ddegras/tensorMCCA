\name{mcca.cor}
\alias{mcca.cor}

\title{Multiple CCA: Maximize Sum of Correlations}

\description{Perform mCCA by maximizing sum of correlations between canonical scores. Optimization method: block coordinate ascent.}

\usage{
mcca.cor(x, r, w = 1, norm = c("block", "global"), 
	ortho = c("score", "weight"), init = c("cca", "svd", "random"), 
	maxit = 1000, tol = 1e-6, sweep = c("cyclical", "random"), 
	control = list(), verbose = FALSE)
}

\arguments{
\item{x}{list of data arrays. All arrays must have the same last dimension (= common measured individuals or objects).}
\item{r}{integer number of sets of canonical components to calculate.}
\item{w}{weight matrix for pairs of datasets. Can be specified as a single number (recycled) or as a square matrix with numbers of rows and columns equal to the length of \code{x}.}
\item{norm}{character specifying the scope of orthogonality constraints for higher-order canonical components: within datasets (\code{block}) or across datasets (\code{global}).}
\item{ortho}{character specifying whether orthogonality constraints for higher-order canonical components apply to canonical scores (\code{score}) or canonical tensors (\code{weight}).}
\item{init}{character describing the initialization method or list of starting points for the optimization.}
\item{maxit}{maximum number of algorithm iterations.}
\item{tol}{numerical (relative and absolute) tolerance for algorithm convergence.}
\item{control}{optional list of arguments for initialization and orthogonality constraints.}
\item{verbose}{logical: display algorithm progress?}
}

\details{
The arrays in \code{x} may have different numbers of dimensions. For example, \code{x} may contain a 2D matrix, a 3D array, and 4D array. The last dimension of these arrays represents the measured individuals or objects, which are the same for all data arrays. Therefore this dimension must be the same in all arrays: \code{tail(dim(x[[1]],1))},  \code{tail(dim(x[[2]],1))}, ..., \code{tail(dim(x[[length(x)]]),1)} should all be equal. 
The arrays in \code{x} are automatically centered along the last dimension prior to optimization. 

If the required number \code{r} of sets of canonical components to calculate is too large,   
it will be reduced before the optimization to the largest theoretical number of possible sets of nonzero canonical components. This theoretical number depends on the data dimensions and the orthogonality constraints. During the optimization, the algorithm stops when the canonical tensors obtained at the end of a stage yield a null objective value. Only solutions associated to positive objective values are retained. 

The objective weight matrix \code{w} must contain only nonnegative values (with at least one positive). It will be made symmetric and scaled so its elements add up to 1.

The arguments \code{norm} and \code{ortho} respectively specify the scope and type of orthogonality constraints used to estimate higher-order canonical components (\eqn{2,...,r}). Details on these constraints can be found in the reference below. 
The constraints are enforced by deflating/reducing the data arrays.  

The argument \code{init} can either be a character specifying an initialization method or  
or a list directly providing starting points for the optimization. In the former case, possible character values are \code{cca} (standard CCA between pairs of matricized datasets), \code{svd} (higher order SVD), or \code{random}. Refer to Degras et al. (2022) for details. If \code{init} is a list, it should either have vector type (no dimensions) and same length as \code{x} OR have matrix type with number of rows equal to the length of \code{x}. In all cases, \code{init[[i]]} or \code{init[[i,l]]} specifies the starting point for the \eqn{i}th dataset and \eqn{l}th round of optimization. Eachof these values should be a list of vectors whose lengths match (all but the last) dimensions of \code{x[[i]]}. If the dimensions of the list \code{init} are less than the effective value of \code{r}, then the values in \code{init} are recycled as necessary.   

 


% % \item{\code{value}}{List of starting values for the MCCA. If \code{value} is provided, the argument \code{method} will be ignored. \code{value} can either have the same length as \code{x}, say \eqn{m}, or have dimensions \eqn{m}-by-\code{r}. In the first case, the same starting point is used for all \code{r} stages of the MCCA. In the second, each column of \code{value} will be used as separate starting point. The \eqn{i}-th field of \code{value} (or all fields of its \eqn{i}-th row) must match the dimensions of the \eqn{i}-th field of \code{x}. More precisely, the \eqn{i}-th field of \code{value} must be a list of length equal to that of \code{x[[i]]} minus 1. If for example \code{x[[i]]} is a 4D array of 3D images of dimensions 10-by-15-by-20, then \code{value[[i]]} must be a list of length 3 with each component a single vector of length 10, 15, and 20, respectively.}

% \item{\code{k}}{Only used for \code{method = cca}. Number of singular values used to approximate each dataset. It can be useful to set \code{k} to a small-to-moderate value if one or several datasets are very big.}

% \item{\code{search}}{Only used for \code{method = cca}. 
% Type of search in the space of possible solutions (size \eqn{m^m} where \eqn{m} is the number of datasets): \code{exhaustive} or \code{approximate}. By default, \code{search} is set to \code{exhaustive} if \eqn{m \le 5} and to 
% \code{approximate} otherwise.} 


}

\value{
A LIST with components
\describe{
\item{\code{v}}{list of canonical tensors with dimensions \eqn{m}-by-\eqn{r} (dataset x component)}
\item{\code{block.score}}{canonical scores: array of dimensions individuals x datasets x canonical components}
\item{\code{global.score}}{average of the canonical scores across datasets: matrix of dimensions individuals x canonical components}
\item{\code{objective}}{objective values (sum of correlations) for each canonical component (vector of length \code{r})}
\item{\code{iters}}{number of algorithm iterations for each canonical component (length \code{r})}
\item{\code{input}}{list of input arguments}
}
}



\seealso{
\code{\link{mcca.cov}}, \code{\link{mcca.init.cca}}, \code{\link{mcca.init.random}}, \code{\link{mcca.init.svd}}
}

\references{
D. Degras, Y. Shen, X. Chu (2022). Multiple canonical correlation analysis for mixed tensor data. \emph{Manuscript in progress.}	
}

% \examples{
% }
