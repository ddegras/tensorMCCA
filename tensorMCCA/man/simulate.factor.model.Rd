\name{simulate.factor.model}
\alias{simulate.factor.model}

\title{
Simulate Tensor Datasets based on Factor Model 
}

\description{
This function generates tensor datasets from the block factor model or the global factor model of Degras et al. (2023)
}

\usage{
simulate.factor.model(dimx, r, scale.v = 1, score.cov = NULL, noise.cov = NULL, 
	ortho.v = c("cyclical", "random", "alldim", "maxdim", "none"))
}

\arguments{
  \item{dimx}{
List of tensor dimensions as integer vectors. The last integer of each list element should be the same for all elements.}
  \item{r}{
Number of canonical components in model.}
\item{scale.v}{
	Norms of canonical weight tensors. Can be specified as a single number (to be recycled) or as a numerical vector of same length as \code{dimx}. 
}
  \item{score.cov}{
Covariance structure of random scores. Can be specified as a numeric vector of length \code{r}, an array of dimensions \eqn{m}-by-\eqn{m}-by-\code{r} with \eqn{m} the length of \code{dimx} (number of datasets), or a list of \code{r} covariance matrices of size \eqn{m}-by-\eqn{m}. 
}
  \item{noise.cov}{
Covariance structure of noise terms. Can be specified as a numeric vector of length \eqn{m} or a list of \eqn{m} covariance matrices of size \eqn{P_i\times P_i}{(Pi,Pi)} 
for \eqn{i=1,...,m,} where \eqn{P_i}{Pi} is the product of all elements in \code{dimx[[i]]} except the last.}
    \item{ortho.v}{modes in which the canonical weight tensors are orthogonal to each other. Can be specified as a string or as an integer vector of same length as \code{dimx}. 
}
}

\details{
The function simulates the global factor model 
\deqn{X_{it} = s_{t}^{(1)} v_{i}^{(1)} + \cdots + s_{t}^{(r)} v_{i}^{(r)} + E_{it}}{X(i,t) = s(t,1) v(i,1) + ... + s(t,r) v(i,r) + E(i,t)}
or the block factor model 
\deqn{X_{it} = s_{it}^{(1)} v_{i}^{(1)} + \cdots + s_{it}^{(r)} v_{i}^{(r)} + E_{it}}{{X(i,t) = s(i,t,1) v(i,1) + ... + s(i,t,r) v(i,r) + E(i,t)}}
where \eqn{i=1,...,m} denotes the dataset and \eqn{t=1,...,n} denotes the subject/case. 

The \eqn{X_{it}}{X(i,t)}, \eqn{t=1,...,n,} are independent random tensors with dimensions given by \code{dimx}, 
say \eqn{p_{i1}\times\cdots\times p_{id_i}}{pi1 x ... x pidi}. \cr
The \eqn{s_{t}^{(l)}}{s(t,l)} (global factor model) or \eqn{s_{it}^{(l)}}{s(i,t,l)}
(block factor model) are random factors or scores. In both models, these random variables have mean zero and are independent across \eqn{t} and \eqn{l}. \cr
The \eqn{v_i^{(l)}}{v(i,l)} are (deterministic) rank-1 tensors representing canonical weights (directions of shared variation across datasets). \cr 
Finally, the \eqn{E_{it}}{E(i,t)} are random noise tensors of general order. 


IF \code{score.cov} is given as a single number (to be recycled) or a vector, then the \eqn{l}th entry in this vector (\eqn{l=1,...,r}) represents the variance \eqn{\sigma_l^2}{sigma(l)^2} of the factor \eqn{s_t^{(l)}}{s(t,l)} and \code{simulate.factor.model} simulates the global factor model. 

IF \code{score.cov} is given as a list of matrices, then the \eqn{l}th matrix (\eqn{l=1,...,r}) is taken as the covariance matrix between the block factors \eqn{s(i,t,l)} and \eqn{s(j,t,l)} where \eqn{i,j=1,...,m}, and \code{simulate.factor.model} simulates the block factor model. 

Note: the global factor model can also be obtained as a special case of the block factor model by making the \eqn{m\times m}{(m,m)} covariance matrix of the \eqn{s_{it}^{(l)}}{s(i,t,l)}, \eqn{i=1,...,m} constant and equal to \eqn{\sigma_l^2}{sigma(l)^2} (degenerate covariance of rank 1). 

The usage of arguments \code{scale.v} and \code{ortho.v} is described in 
\code{\link{simulate.v}} (arguments \code{scale} and \code{ortho} there). 

}

\value{
A LIST with components 
\describe{
\item{x}{List of tensor datasets with dimensions given by \code{dimx}.}
\item{v}{List of dimensions \eqn{n\times r}{(m,r)} containing canonical weights. Each element is a list of vectors whose outer product is a rank-1 tensor of canonical weights.}
\item{score}{Array of scores/factors of size \eqn{n\times m\times r}{(n,m,r)} where \eqn{n} is the last integer of each element in \code{dimx} (number of samples).}
}

Note: for the global factor model, the scores are identical along the second dimension (\eqn{s_{1t}^{(l)}=\ldots=s_{mt}^{(l)}:=s_{t}^{(l)}}{s(1,t,l)=...=s(m,t,l):=s(t,l)}). 
}

\references{
%% ~put references to the literature/web site here ~
}

\examples{
n <- 100 # number of samples
m <- 3 # number of datasets
p <- list(10, c(5,10), c(6,8,10)) # data dimensions
d <- sapply(p, length) # tensor order
dimx <- lapply(p, append, value = n) # full data dimensions
r <- 4 # number of factors

## Example 1: global factor model
## X(i,t) = s(t,1) v(i,1) + ... + s(t,r) v(i,r) + E(i,t) 
score.cov <- seq(4, 1, len = r) # factor/score variances 
noise.cov <- 0.5 
# Each tensor E(i,t) is a white noise with variance 0.5 (recycled)
scale.v <- 1 # set norm of each canonical weight tensor v(i,l) to 1 (default) 
ortho.v <- 1 
# Impose orthogonality in first mode between v(i,l) and v(i,l') for l != l'
sim1 <- simulate.factor.model(dimx, r, scale.v, score.cov, noise.cov, ortho.v)
x <- sim1$x
v <- sim1$v
score <- sim1$score

## Example 2: block factor model 
## X(i,t) = s(i,t,1) v(i,1) + ... + s(i,t,r) v(i,r) + E(i,t) 
scale.v <- c(1.2, 1, 0.8) 
# || v(1,l) || = 1.2, || v(2,l) || = 1, || v(3,l) || = 0.8
score.cov <- lapply(c(0.1, 0.2, 0.3) / m, 
	function(sig2) sig2 * crossprod(matrix(rnorm(m^2), m, m))) 
noise.cov <- list(diag(0.5, prod(p[[1]])), 
	diag(0.5, prod(p[[2]])), diag(1, prod(p[[3]]))) # white noise covariances
ortho.v <- "alldim" 
# For each dataset, canonical weights are orthogonal to each other in all modes 
sim2 <- simulate.factor.model(dimx, r, scale.v, score.cov, noise.cov, ortho.v)
x <- sim2$x
v <- sim2$v
score <- sim2$score
}

