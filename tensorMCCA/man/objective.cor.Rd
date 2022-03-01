\name{objective.cor}
\alias{objective.cor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the Sum of Correlations Between Canonical Scores
}

\description{
Calculate the sum of correlations between canonical scores for given datasets, canonical tensors, and weight matrices.}

\usage{
objective.cor(x, v, c = 1)
}


\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{c}{
Weights: matrix of objective weights for pairs of datasets.}
}
\details{
The weight matrix \code{c} can be specified as a single number or as a square matrix with as many rows and columns as there are datasets, i.e. \code{length(x)}. This matrix must have non-negative entries. It is rescaled so that its values add up to one. 
}

\value{The sum of correlations between canonical scores.}

\seealso{
\code{objective.cov}
}


\examples{
x <- list(array(runif(60), c(2,3,10)), array(runif(160), c(4,2,2,10)))
v <- init.v.random(x)
objective.cor(x, v)
}
