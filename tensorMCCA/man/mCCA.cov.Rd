\name{mCCA.cov}
\alias{mCCA.cov}

\title{Multiple CCA: Maximize Sum of Covariances}

\description{Perform multiple CCA by maximizing sum of covariances between canonical components.}

\usage{
mCCA.cov(x, r, c = 1, init.type = c("svd", "ones", "random"), init.value = NULL,
	cnstr = c("block", "global"), ortho = c("block.score", "global.score", "canon.t.1", "canon.t.all"), 
	maxit = 1000, tol = 1e-6, sweep = c("cyclical", "random"), verbose = FALSE)
}

\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
\item{r}{Number of canonical components to calculate.}
\item{c}{Weight matrix for pairs of datasets.}
\item{init.type}{Initialization method.}
\item{init.value}{Optional starting point for optimization.}
\item{cnstr}{Norm constraints for canonical tensors.}
\item{ortho}{Orthogonality constraints for higher order canonical components.}
\item{maxit}{Maximum number of algorihm iterations.}
\item{tol}{Tolerance for algorithm convergence.}
\item{sweep}{Sweeping method}
\item{verbose}{Display algorithm progress?}
}

% \details{
% }

\value{
A LIST with components
\describe{
\item{\code{v}}{canonical tensors}
\item{\code{block.score}}{canonical scores: array of dimensions individuals x datasets x canonical components}
\item{\code{global.score}}{weighted average of the canonical scores with weights given by the column sums of \code{c} (individuals x canonical components)}
\item{\code{objective}}{objective value (sum of correlations) for each canonical component (vector of length \code{r})}
\item{\code{iters}}{Number of algorithm iterations for each canonical component (length \code{r})}
\item{\code{c}}{Input argument \code{c}}
\item{\code{maxit}}{Input argument \code{maxit}}
\item{\code{tol}}{Input argument \code{tol}}
\item{\code{init.type}}{Input argument \code{init.type}}
\item{\code{init.value}}{Input argument \code{init.value}}
\item{\code{ortho}}{Input argument \code{ortho}}
}
}



\seealso{
\code{\link{mCCA.cor}}
}

% \examples{
% }
