\name{kkt}
\alias{kkt}

\title{
Calculate Residuals in the Karush-Kuhn-Tucker (KKT) Equations
}

\description{
Calculate the residuals in the KKT equations for tensor MCCA
}

\usage{
kkt(x, object)
}

\arguments{
  \item{x}{
List of tensor/array datasets
}
  \item{object}{
Result of a call to \code{\link{mcca.cov}} or \code{\link{mcca.cor}} 
}
}

\details{
At present, the function only calculates the KKT residuals for the **first** canonical components. 
}

\value{
A LIST with the same structure as the first set of canonical weights in \code{object}, that is, \code{object$v[,1]}.
}

\references{
%% ~put references to the literature/web site here ~
}

\seealso{
\code{\link{mcca.cov}}, \code{\link{mcca.cor}} 
}

\examples{
}

