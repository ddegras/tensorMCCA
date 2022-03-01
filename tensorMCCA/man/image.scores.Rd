\name{image.scores}
\alias{image.scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Individual Scores Along Canonical Tensors
}
\description{
Calculate individual canonical scores for given datasets and canonical tensors. 
}
\usage{
image.scores(x, v)
}

\arguments{
  \item{x}{Data: list of arrays.}
  \item{v}{Canonical tensors: list of lists of vectors.}
}

\details{
Each component in the list \code{x} should be a data array with the last mode representing individuals or objects. Each component in the list \code{v} should be a sublist of vectors specifying a rank-1 canonical tensor. The lists \code{x} and \code{v} should have the same length. Also, the lengths of the vectors in \code{v} must match the corresponding dimensions in the components of \code{x}. Finally, the arrays in \code{x} should all have the same number of dimensions in their last mode. 
}

\value{
Matrix of canonical scores with individuals in rows and datasets in columns.
}

\examples{
}

