\name{mcca.boot.ci}
\alias{mcca.boot.ci}
\title{
Bootstrap Confidence Intervals for Multiple Canonical Correlation Analysis
}
\description{
Calculate pointwise bootstrap confidence intervals (basic, percentile, normal,
and bias-corrected normal)
for various quantities involved in tensor mCCA: 
canonical weights, covariance/correlation between canonical scores, 
covariance/correlation structure of noise. 
}

\usage{
mcca.boot.ci(object, level = 0.95,
	type = c("basic", "normal", "normal-debiased", "percentile"))
}
\arguments{
  \item{object}{
The result of a call to \code{\link{mcca.boot}}, containing successful
bootstrap replicates with the same targets and dimensions as the original estimates.
}
  \item{level}{
A finite confidence level strictly between 0 and 1 (default = 0.95).
}
  \item{type}{
type of CIs: any combination of the strings \code{basic}, \code{normal}, \code{normal-debiased} and \code{percentile}. 
Partial matching is enabled. 
}
}

\details{
For each inferential target, the function calculates the bootstrap mean minus
the original estimate as its estimated bias. The standard error is the square
root of the mean squared deviation from the bootstrap mean, with the number of
replicates as denominator.

Basic intervals reflect the bootstrap quantiles about the original estimate.
Percentile intervals use the bootstrap quantiles directly. Normal intervals
are centered at the original estimate. The \code{normal-debiased} intervals
subtract the estimated bias from that center. Both normal methods use the
bootstrap standard error and standard normal quantiles.

Correlation matrices have fixed unit diagonals. Their bias and standard error
matrices therefore have zero diagonals, and their interval bounds have unit
diagonals. Failed, incomplete, or non-finite bootstrap replicates produce an
error identifying the affected replicate; they are not silently omitted.
}


\value{
A list with possible fields \code{v}, \code{block.score.cov},
\code{global.score.cov}, \code{block.score.cor}, \code{global.score.cor},
\code{noise.cov}, and \code{noise.cor}, matching the target names in
\code{object$original}. Each target contains \code{bias}, \code{se}, and the
requested interval fields among \code{basic}, \code{percentile}, \code{normal},
and \code{normal-debiased}. Interval arrays have a final dimension of length
two for the lower and upper bounds.
}

\references{
\url{https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Methods_for_bootstrap_confidence_intervals}
}

\seealso{
\code{\link{mcca.boot}}
}

\examples{
set.seed(42)
x <- list(matrix(rnorm(3 * 30), 3, 30),
          matrix(rnorm(4 * 30), 4, 30))
mcca.fit <- mcca.cov(x, r = 1, init = "svd")
boot.object <- mcca.boot(x, mcca.fit, target = "v", nboot = 19,
                         parallel = FALSE)
boot.ci <- mcca.boot.ci(boot.object)
boot.ci$v[[1,1]]$basic
boot.ci$v[[1,1]][["normal-debiased"]]
}

