\name{mCCA.cor}
\alias{mCCA.cor}

\title{Multiple CCA: Maximize Sum of Correlations}

\description{Perform multiple CCA by maximizing sum of correlations between canonical components.}

\usage{
mCCA.cor(x, r, c = 1, maxit = 1000, tol = 1e-06, init.type = c("svd", "ones", "random"), init.value = NULL, ortho = c("block.score", "global.score", "canon.t.1", "canon.t.all"), balance = TRUE, verbose = FALSE)
}

\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
\item{r}{Number of canonical components to calculate.}
\item{c}{Weight matrix for pairs of datasets.}
\item{maxit}{Maximum number of algorihm iterations.}
\item{tol}{Tolerance for algorithm convergence.}
\item{init.type}{Initialization method.}
\item{init.value}{Optional starting point for optimization.}
\item{ortho}{Orthogonality constraints for canonical components.}
\item{balance}{Should canonical vectors be rescaled to have equal norms within each dataset?}
\item{verbose}{Display algorithm progress?}
}

% \details{
% }

\value{
A list with components
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% \examples{
% }
