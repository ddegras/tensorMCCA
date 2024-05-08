\name{hessian}
\alias{hessian}

\title{
Hessian Matrix of Lagrange Function And Its Projection
}
\description{
Calculate the Hessian matrix of the Lagrange function for MCCA solutions. 
One Hessian matrix is calculated for each set of canonical weights. 
The projection of the Hessian matrix on the orthogonal of the gradients of the constraint functions is also provided to help determine whether a solution is locally optimal. 
}

\usage{
hessian(x, fit, grad = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
List of data arrays. Each data array should have the same size in its last dimension.
}

\item{fit}{
Fitted MCCA object: the result of a call to \code{\link{mcca.cor}} or \code{\link{mcca.cov}}.
}

\item{grad}{
Object containing gradients and Lagrange multipliers associated with MCCA: 
the result of a call to \code{\link{kkt}}. }
}

\details{

If not provided, \code{grad} will be calculated inside the function via the call \code{kkt(x, fit)}.  

A necessary (second-order) condition for a MCCA solution \eqn{\hat{v}} to be a local maximizer is that for all vector \eqn{v} orthogonal to the gradients of the constraints at \eqn{\hat{v}}, \eqn{v^{T}Hv \le 0} where \eqn{H} is the Hessian of the objective function at \eqn{\hat{v}}. 

If \eqn{v^{T}Hv < 0} for all such (nonzero) vector \eqn{v}, then \eqn{\hat{v}} is a local maximizer (sufficient condition). A convenient way to assess this condition is to 1) find an orthogonal basis of the orthogonal space of the gradient constraints (e.g. by QR decomposition) and collect it in a matrix \eqn{Q}, 2) multiply \eqn{H} left and right by \eqn{Q^T} and \eqn{Q}, and 3) check that \eqn{Q^T H Q} is negative definite (i.e. has negative eigenvalues).  


}
\value{
List with components
\item{H}{Hessian matrices (\eqn{P\times P \times r} array)}
\item{projH}{List of projected Hessian matrices (length \eqn{r})}
}

\references{
%% ~put references to the literature/web site here ~
}


\seealso{
\code{\link{kkt}}
}

\examples{}

