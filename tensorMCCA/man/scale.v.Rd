\name{scale.v}
\alias{scale.v}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale Canonical Tensors
}
\description{
Scale canonical tensors according to (unit) norm or variance constraints. 
}

\usage{
scale.v(v, x = NULL, type = c("norm", "var"), scope = c("block", "global"), balance = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{Canonical tensors: }
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
\item{type}{Type of scaling constraints: unit norm or unit variance.}
\item{scope}{Scope of scaling constraints: block-level (i.e. for each dataset) or global.}
\item{balance}{Should the vectors making up the rank-1 tensors have equal norm?}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
A list of canonical tensors of the same length as \code{v}.
}

\examples{
}
