\name{mcca.boot}
\alias{mcca.boot}
\title{Bootstrap Multiple Canonical Correlation Analysis}
\description{
Calculate bootstrap distributions of canonical weights, canonical score
covariance or correlation matrices, and noise covariance or correlation matrices.
}
\usage{
mcca.boot(x, object,
    target = c("v", "score.cov", "noise.cov", "score.cor", "noise.cor"),
    resample = c("data", "residuals"), nboot = 100, parallel = TRUE,
    control = list())
}
\arguments{
  \item{x}{List of numeric data arrays with the same size in the last
    dimension, which indexes observations.}
  \item{object}{The result of fitting \code{\link{mcca.cor}} or
    \code{\link{mcca.cov}} to \code{x}.}
  \item{target}{Any nonempty subset of \code{"v"}, \code{"score.cov"},
    \code{"noise.cov"}, \code{"score.cor"}, and \code{"noise.cor"}.
    Partial matching is enabled.}
  \item{resample}{Use \code{"data"} to resample observations or
    \code{"residuals"} to resample residuals while retaining fitted values.}
  \item{nboot}{Positive integer number of bootstrap replications.}
  \item{parallel}{Whether to use a registered \pkg{foreach} parallel backend.
    If no backend is registered, computations run sequentially with a warning.}
  \item{control}{List of options. \code{.errorhandling} is \code{"stop"}
    (default), \code{"remove"}, or \code{"pass"}. The logical options
    \code{unbiased.score.cov} and \code{unbiased.noise.cov} select the existing
    covariance adjustments and both default to \code{TRUE}.}
}
\details{
Each replication uses one common observation index vector across all datasets.
For residual resampling, the same index vector is applied to each residual
matrix and the resampled residuals are added to the retained fitted values.
Refits retain the original fit's options, including explicit initial weights
when these were supplied.

One common sign is chosen for each component by the inner product between
its stacked bootstrap and original weight tensors. That sign is applied across
all datasets. Cross-component global score covariances and correlations are
transformed by the same signs. Components retain their fitted output order;
this function does not perform an additional component assignment.

The error-handling option is supplied inside \code{control}. With
\code{"stop"}, the first failed replication raises an error. With
\code{"remove"}, failed replications are omitted, so the returned list can
be shorter than \code{nboot}, including empty when all replications fail.
With \code{"pass"}, error objects remain in the bootstrap list and are
not sign-aligned. Inspect or remove failed replications before constructing
confidence intervals.
}
\value{
A list with components:
\describe{
  \item{bootstrap}{The bootstrap estimates, subject to the selected error policy.}
  \item{original}{Corresponding estimates from the original fitted object.}
}
Estimate fields are \code{v}, \code{block.score.cov},
\code{global.score.cov}, \code{block.score.cor}, \code{global.score.cor},
\code{noise.cov}, and \code{noise.cor}, according to \code{target}.
Covariance fields store the lower triangle including its diagonal;
correlation fields store the strict lower triangle. Block score fields have
one column per fitted component; noise fields are lists indexed by dataset.
}
\seealso{\code{\link{mcca.boot.ci}}}
\examples{
set.seed(42)
n <- 40
signal <- rnorm(n)
x <- lapply(1:2, function(i)
    matrix(rnorm(3 * n), 3, n) + rep(signal, each = 3))
fit <- mcca.cov(x, init = "svd", maxit = 40)
boot <- mcca.boot(x, fit, target = c("v", "score.cov"),
                  nboot = 4, parallel = FALSE)
residual.boot <- mcca.boot(x, fit, target = "v", resample = "residuals",
                           nboot = 4, parallel = FALSE)
stopifnot(length(boot$bootstrap) == 4,
          length(residual.boot$bootstrap) == 4)
}
