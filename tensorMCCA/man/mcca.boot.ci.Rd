\name{mcca.boot.ci}
\alias{mcca.boot.ci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap Confidence Intervals for Multiple Canonical Correlation Analysis
}
\description{
Calculate standard bootstrap CIs (basic, percentile, normal) 
for various quantities involved in tensor mCCA: 
canonical weights, covariance/correlation between canonical scores, 
covariance/correlation structure of noise. 
}

\usage{
mcca.boot.ci(object, level = 0.95, type = c("basic", "normal", "percentile"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
The result of a call to the function \code{\link{mcca.boot}}}. 
  \item{level}{
confidence level (default = 0.95)
}
  \item{type}{
type of CIs: any combination of the strings \code{basic}, \code{normal}, and \code{percentile}. 
Partial matching is enabled. 
}
}

\details{
For each of the inferential targets stored in \code{object}, the function calculates bootstrap estimates of the estimation bias and standard error, as well as pointwise confidence intervals based on the basic, percentile, or normal method. The normal CIs are bias-corrected. 
}


\value{
A list whose possible fields \code{v}, \code{score.cov.block}, \code{score.cov.global}, code{score.cor.block}, \code{score.cor.global}, \code{noise.cov}, and \code{noise.cor}. Each of these fields represents a target for bootstrap inference and contains some or all of these subfields according to the contents of \code{object}: \code{bias}, \code{se}, \code{basic}, \code{percentile}, \code{normal}. 
}

\references{
\url{https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Methods_for_bootstrap_confidence_intervals}
}

\seealso{
\code{\link{mcca.boot}}
}

\examples{
## Simulate data
n <- 100 # number of samples
m <- 3 # number of datasets
p <- list(10, c(5,10), c(6,8,10)) # data dimensions
d <- sapply(p, length) # tensor order
dimx <- lapply(p, append, value = n) # full data dimensions
r <- 4 # number of factors
score.cov <- seq(4, 1, len = r) # score variance 
noise.cov <- 0.5 # noise variance (will be recycled)
x <- simulate.factor.model(dimx, r, score.cov, noise.cov, ortho = TRUE)$x

## Run mCCA
mcca.fit <- mcca.cor(x, r)

## Perform bootstrap based on raw data
library(foreach)
library(doParallel)
registerDoParallel()
boot.object <- mcca.boot(x, mcca.fit)

## Calculate bootstrap CIs
boot.ci <- mcca.boot.ci(boot.object) 

## Examine some results
boot.ci$score.cov.block
# Because the simulation model utilizes independent global scores, 
# it is liekly that the above CIs for off-diagonal covariance terms
# should contain zero and those for diagonal terms (score variances) 
# should not contain zero 	
}

