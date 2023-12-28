\name{kkt}
\alias{kkt}

\title{
Solve Karush-Kuhn-Tucker (KKT) Equations for Tensor MCCA
}

\description{
Returns the gradients of the objective and constraint functions (scaling + orthogonality constraints) and the associated Lagrange multipliers and residuals of the 
the KKT equations for tensor MCCA.
}

\usage{
kkt(x, fit)
}

\arguments{
  \item{x}{
List of tensor/array datasets.
}
  \item{fit}{
Result of a call to \code{\link{mcca.cov}} or \code{\link{mcca.cor}} on \code{x}. 
}
}

\details{
The KKT equations are first-order necessary conditions for a solution to a constrained optimization problem to be locally optimal. 
	
The function attempts to solve the KKT equations of the MCCA problem whose objective, constraints, and solutions are in \code{fit}. For each set of canonical components, after the gradients of the objective and constraint functions have been calculated, the former are regressed against the latter by ordinary least squares (OLS). The regression coefficients are taken as Lagrange multipliers and the regression residuals are termed KKT residuals. 

The most important output for the user is probably the KKT residuals. Indeed, a first-order necessary condition for an MCCA solution to be locally optimal is that these residuals are close to zero (up to optimization error). The other way around, if the KKT residuals are not close to zero, then the solution found is certainly not a local optimizer. 

Other outputs of the function are the gradient of the objective function, the gradients of the constraint functions, and the associated Lagrange multipliers for each stage of optimization. The two latter quantities are necessary to calculate the Hessian of the Lagrange function and its projection (see function \code{\link{hessian}}).  

Let \eqn{m} be the number of tensor datasets and \eqn{r} be the number of canonical components calculated in the MCCA. Let's also denote by \eqn{(p_{i1},\ldots,p_{id_i},n)} the dimensions of the \eqn{i}th dataset (\eqn{i=1,\ldots,m}). The MCCA objective is parameterized as a function of \eqn{v = (v_{11},\ldots, v_{1d_1},\ldots, v_{md_m})} of length \eqn{P=\sum_{i=1}^m \sum_{k=1}^{d_i} p_{ik}}. 


}

\value{
A list with components
\item{grad.objective}{Gradient of the objective function for each set of canonical components (\eqn{P\times r} matrix).}
\item{grad.constraints}{Gradients of the constraint functions for each set of canonical components. List of length \eqn{r} whose \eqn{l}th component is a matrix with \eqn{P} rows and \eqn{ml} columns, resp. \eqn{l} columns, if the constraints are block-specific, resp. global.}
\item{multipliers}{List of Lagrange multipliers (length\eqn{r}).}
\item{residuals}{KKT residuals (\eqn{P\times r} matrix).}

In the matrices \code{grad.objective} and \code{residuals}, the \eqn{l}th column corresponds to the \eqn{l}th stage of optimization/set of canonical components (\eqn{l=1,\ldots, r})

}

\references{
Degras, Chu, and Shen (2023). Multiple Canonical Correlation Analysis for Tensor Data. \emph{In progress.}

Nocedal and Wright (2006). Numerical Optimization [Chapter 12]. 
}

\seealso{
\code{\link{mcca.cov}}, \code{\link{mcca.cor}}, \code{\link{hessian}}
}

\examples{
}

