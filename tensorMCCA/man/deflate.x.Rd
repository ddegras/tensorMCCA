\name{deflate.x}
\alias{deflate.x}

\title{
Deflate Data Arrays
}

\description{
Deflate data arrays with respect to canonical tensors or canonical scores.
}

\usage{
deflate.x(x, v = NULL, y = NULL, ortho = c("block.score", "global.score", "canon.t.1", "canon.t.all"), check.args = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
}
\item{y}{
Matrix of canonical scores (individuals/objects x datasets) 
}

\item{ortho}{
Type of orthogonality constraints used for deflation. 
}

\item{check.args}{
Should the arguments \code{x} and \code{v} be checked? 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
A list of deflated data arrays of the same size as \code{x}
}

\seealso{
\code{\link{deflate.v}}
}

\examples{
}
