\name{simulate.factor.model}
\alias{simulate.factor.model}

\title{
Simulate Factor Model for Multiple Tensor Data
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
simulate.factor.model(dimx, r, score.cov, noise.cov, 
	scale = c("block", "global"), ortho = FALSE)
}

\arguments{
  \item{dimx}{
List of tensor dimensions as integer vectors. The last integer of each list element should be the same for all elements.}
  \item{r}{
Number of canonical components in model.}
  \item{score.cov}{
Covariance structure of random scores. Can be specified as a numeric vector of length \code{r}, an array of dimensions \eqn{m}-by-\eqn{m}-by-\code{r} with \eqn{m} the length of \code{dimx} (number of datasets), or a list of \code{r} covariance matrices of size \eqn{m}-by-\eqn{m}. 
}
  \item{noise.cov}{
Covariance structure of noise terms. Can be specified as a numeric vector of length \eqn{m} or a list of \eqn{m} covariance matrices of size \eqn{P[i]}-by-\eqn{P[i]} 
for \eqn{i=1,...,m,} where \eqn{P[i]} is the product of all elements in \code{dimx[[i]]} except the last.}
  \item{scale}{
Type of scaling for canonical weights: 
scale to unit norm for each dataset (\code{block})
or scale the weights so their squared norm is 1 on average across datasets (\code{global}). 
}
  \item{ortho}{
Logical: should the canonical weights be orthogonal? 
}
}

\details{
The model to simulate is 
\deqn{ X[i,t] = s(i,t,1) v(i,1) + ... + s(i,t,r) v(i,r) + E(i,t) }
Through proper specification of the score covariances, the scores \eqn{s(1,t,l), ... , s(m,t,l)} can be made equal (unique score for each sample and canonical component).



}
\value{
A LIST with components 
\describe{
\item{x}{Data tensors with dimensions given by \code{dimx}.}
\item{v}{Canonical weights as a list of size \eqn{n}-by-\code{r}. Each list element is itself a list specifying the vector components of the canonical weights which are rank-1 tensors.}
\item{score}{Scores in an array of size \eqn{n}-by-\eqn{m}-by-\code{r} where \eqn{n} is the last integer of each element in \code{dimx} (number of samples).}
}
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

## Model 1: X(i,t) = s(i,t,1) v(i,1) + ... + s(i,t,r) v(i,r) + E(i,t) 
score.cov <- seq(4, 1, len = r) # score variance 
noise.cov <- 0.5 # noise variance (will be recycled)
sim1 <- simulate.factor.model(dimx, r, score.cov, noise.cov, "block", TRUE)
x <- sim1$x
v <- sim1$v
score <- sim1$score

## Model 2: X(i,t) = s(i,t,1) v(i,1) + ... + s(i,t,r) v(i,r) + E(i,t) 
score.cov <- lapply(sqrt(seq(4, 1, len = r)), 
	function(sig) crossprod(matrix(rnorm(m^2, sd = sig), m, m))) 
noise.cov <- list(diag(0.5, prod(p[1])), 
	diag(0.5, prod(p[1])), diag(1, prod(p[1]))) # white noise 
sim2 <- simulate.factor.model(dimx, r, score.cov, noise.cov, "block", TRUE)
x <- sim2$x
v <- sim2$v
score <- sim2$score


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
