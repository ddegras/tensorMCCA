\name{init.mcca.random}
\alias{init.mcca.random}

\title{
Initialize Canonical Tensors Randomly
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
init.mcca.random(x, objective = c("cov", "cor"), cnstr = c("block", "global"), 
	center = TRUE, balance = TRUE)
}

\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
\item{objective}{Objective function: sum of covariances or sum of correlations?}
\item{scope}{Scaling of canonical tensors: block constraints (i.e. for each dataset) or global constraints?}
\item{balance}{Should the canonical vectors have equal norm?}

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
A list representing canonical tensors. Each sublist contains as many canonical vectors as there are modes in the corresponding dataset (not counting the last mode which corresponds to individuals or objects). The vectors of a sublist have length equal to the numbers of dimensions in the modes of the corresponding dataset.}

\seealso{
\code{\link{init.mcca.cca}}, \code{\link{init.mcca.svd}}
}

\examples{
}
