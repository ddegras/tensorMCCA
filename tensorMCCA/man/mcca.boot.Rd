\name{mcca.boot}
\alias{mcca.boot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap Multiple Canonical Correlation Analysis
}
\description{
Calculate the bootstrap distribution of estimators of the canonical weights, 
the covariance/correlation matrices of the canonical scores, 
and covariance/correlation matrices of the noise.
}

\usage{
mcca.boot(x, object, 
	target = c("v", "score.cov", "noise.cov", "score.cor", "noise.cor"),
	resample = c("data", "residuals"), nboot = 100, parallel = TRUE, 
	control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
List of data arrays. All arrays must have the same size in their last dimension (= common instances, cases, or subjects).
}
  \item{object}{
The result of a call to \code{\link{mcca.cor}} or \code{\link{mcca.cov}} on \code{x}.
}
  \item{resample}{
Type of boostrap: based on raw \code{data} or on \code{residuals}.
}
  \item{nboot}{
Number of bootstrap replications. Default = 100.}
 
 \item{parallel}{
Logical value specifying whether calculations should be run in parallel (\code{TRUE}, default) or not. A back-end for parallel computations must be registered before calling this function. 
}
  \item{target}{
  	Targets for bootstrap inference. Any subset of "v", "score.cov", "noise.cov", "score.cor", and "noise.cor". Partial matching is enabled. 
}
\item{.errorhandling}{
	Error handling in the bootstrap loop.
}
}
\details{
The result of a call to \code{mcca.boot} is not meant to be looked at directly. 
Rather, it served as a intermediate step before calling the function \code{\link{mcca.boot.ci}}. 

If \code{method} = \code{raw}, the entire data are resampled (with replacement). If \code{method} = \code{residual}, all fitted values are retained while the residuals are resampled (with replacement) and added to the fitted values. 

If \code{parallel} = \code{TRUE} and execution in parallel is possible, the argument \code{.errorhandling} is passed as-is to the function \code{foreach} from package \code{foreach}. Otherwise the same functionalities (stop on errors, pass errors and keep error messages, or remove results containing errors) are implemented in the serial loop. 


}
\value{
A LIST with components 
\describe{
\item{bootstrap}{list of length \code{nboot} containing the required bootstrap estimates}
\item{original}{corresponding list of estimates based on the original data}
}
The sublist \code{original} and each field of the sublist \code{bootstrap} contain the same field names which may be: \code{v}, \code{score.cov.block}, \code{score.cov.global}, \code{score.cor.block}, \code{score.cor.global}, \code{noise.cov}, and \code{noise.cor}.
}


\seealso{
\code{\link{mcca.boot.ci}}
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

## Example 1: simple bootstrap based on raw dat
boot.object <- mcca.boot(x, mcca.fit, nboot = 20)

## Example 2: residual-based bootstrap with parallel computation
library(foreach)
library(doParallel)
registerDoParallel()
boot.object <- mcca.boot(x, mcca.fit, method = "resid")

## Examine the output
names(boot.object)
names(boot.object$original)
stopifnot(length(boot.object$bootstrap) == 100)
}
