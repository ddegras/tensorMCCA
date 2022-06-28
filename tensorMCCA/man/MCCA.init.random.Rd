\name{mcca.init.random}
\alias{mcca.init.random}

\title{
Initialize Canonical Tensors Randomly
}

\description{
Generate random rank-1 tensors of given size under scaling constraints and optional orthogonality constraints.}

\usage{
mcca.init.random(x, r = 1L, scale = c("norm", "var"), ortho.mode = NULL)}

\arguments{
\item{x}{list of data arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
\item{r}{number of sets of canonical tensors to generate.}
\item{scale}{character describing how to scale the canonical tensors: with unit norm (\code{norm}) or unit variance (\code{var})}
\item{ortho.mode}{integer or list indicating in which mode(s) -if any- to orthogonalize the canonical tensors of each dataset. This argument is recycled as needed to a list with same length as \code{x}.}

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
A list representing canonical tensors. Each sublist contains as many canonical vectors as there are modes in the corresponding dataset (not counting the last mode which corresponds to individuals or objects). The vectors of a sublist have length equal to the numbers of dimensions in the modes of the corresponding dataset.}

\seealso{
\code{\link{mcca.init.cca}}, \code{\link{mcca.init.svd}}
}

\examples{
}
