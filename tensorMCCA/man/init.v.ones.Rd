\name{init.v.ones}
\alias{init.v.ones}

\title{
Initialize Canonical Tensors With Ones
}

\description{
Given a list of datasets, create rank-1 tensors with compatible dimensions filled with ones and rescale the tensor according to norm or variance constraints. 
}

\usage{
init.v.ones(x, type = c("norm", "var"), scope = c("block", "global"), balance = TRUE)
}

\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
\item{type}{Type of scaling constraints: unit norm or unit variance.}
\item{scope}{Scope of scaling constraints: block-level (i.e. for each dataset) or global.}
\item{balance}{Should the vectors making up the rank-1 tensors have equal norm?}
}

\details{
The argument \code{type} specifies whether (unit) norm constraints or (unit) variance constraints should be applied to the created tensors. The argument \code{scope} specifies whether the unit norm constraint (if any) should be applied to each tensor separately or to the average of tensor norms. 
}

\value{
A list of the same length as \code{x} containing rank-1 tensors. Each list component is itself a list of vectors whose outer product form a rank-1 tensor. 
}

\seealso{
\code{\link{init.v.random}}, \code{\link{init.v.svd}}
}

% \examples{
% }

