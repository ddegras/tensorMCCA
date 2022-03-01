\name{init.v.svd}
\alias{init.v.svd}

\title{Initialize Canonical Tensors by SVD (Method 1)}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
init.v.svd(x, type = c("norm", "var"), scope = c("block", "global"), balance = TRUE)
}

\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}
\item{r}{Number of canonical components to calculate.}
\item{type}{Type of scaling constraints: unit norm or unit variance.}
\item{scope}{Scope of scaling constraints: block-level (i.e. for each dataset) or global.}
\item{balance}{Should the vectors making up the rank-1 tensors have equal norm?}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
A list of \code{r} sublists representing canonical tensors. Each sublist contains as many canonical vectors as there are modes in the corresponding dataset (not counting the last mode which corresponds to individuals or objects). The vectors of a sublist have length equal to the numbers of dimensions in the modes of the corresponding dataset.}

\seealso{
\code{\link{init.v.ones}}, \code{\link{init.v.random}}
}

\examples{
}

