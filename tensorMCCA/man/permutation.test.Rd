\name{permutation.test}
\alias{permutation.test}

\title{
Permutation Test for Nonzero Canonical Components 
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
permutation.test(x, obj, nperm = 100, parallel = TRUE, ncores = NULL)
}

\arguments{
\item{x}{Data: list of arrays. All arrays must have the same number of dimensions in their last mode (= individual or object).}

\item{obj}{The result of a call to \code{\link{mcca.cor}} or \code{\link{mcca.cov}}.}

\item{nperm}{Number of random permutations to use for the test.}

\item{parallel}{Should calculations be conducted in parallel?}

\item{ncores}{Number of cores for parallel computation.}
}

\details{
The permutation tests are conducted sequentially over canonical components. 
To test the first component, the tensor data are each permuted separately in their common mode (the last). After that, mCCA is fitted to the permuted data using the same parameters as those in \code{obj}, yielding a single set of canonical components. The corresponding objective value (sum of canonical covariances or correlations) is recorded. The process is repeated \code{nperm} times. 

To test the \eqn{l}th component (\eqn{l>1}), the tensor datasets are first deflated to remove the influence of previous canonical components in the fit. After that, each dataset is permuted independently in its last mode and one round of mCCA is conducted. The objective value is recorded and the process is repeated \code{nperm} times. 
 
For each component (\eqn{l=1,...,r}), the uncorrected p-value is calculated as the proportion of permutations for which the (absolute) objective value exceeds the corresponding objective value for the original data. This leads to \eqn{r} p-values where \eqn{r} is the number of canonical components obtained in the original fit. 
The corrected p-values are obtained by taking the cumulative maximum of the uncorrected values. These values provide an admissible testing procedure (later components cannot be retained as significant if earlier ones have been rejected) and yield -in principle- accurate significance level. 

For large datasets, it may be useful to first carry out the permutation test with a small number of permutations. This may help calibrate the number of canonical components to test with a larger number of permutations - and reduce computation time. 
}

\value{
A LIST with components: 
\describe{
\item{pval.corrected}{vector of corrected p-values for each fitted canonical component}
\item{pval.uncorrected}{vector of uncorrected p-values for each fitted canonical component} 
\item{objective}{vector of objective values for original mCCA fit (one per canonical component)}
\item{objective.perm}{matrix of objective values for all permutations and canonical  components}
}
}

\references{
Winkler et al (2020). Permutation inference for canonical correlation analysis. \doi{10.1016/j.neuroimage.2020.117065}	
}


\seealso{
\code{\link{mcca.cov}}, \code{\link{mcca.cor}}
}

% \examples{
% }
