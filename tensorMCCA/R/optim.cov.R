###########################################
# FUNCTIONS TO MAXIMIZE SUM OF COVARIANCES 
# IN MCCA UNDER GLOBAL CONSTRAINTS
###########################################


###############################
# Solve max { v'AA'v + 2 b'v } 
# under constraint v'v = c
###############################


optim1D.cov <- function(A, b, c)
{
eps <- 1e-14
n <- ncol(A)
p <- nrow(A)

## SVD of A
svdA <- svd(A, nu = p, nv = 0)
P <- svdA$u
delta <- svdA$d^2
delta[delta < eps] <- 0
udelta <- unique(delta)
b <- crossprod(P, b)
idx.bz <- which(abs(b) < eps)

## Trivial case 
if (length(idx.bz) == p) {
	return(sqrt(c) * P[,1])
}

## Inspect cases where Lagrange multiplier is among the 
## (squared) singular values of A
if (length(idx.bz) > 0 && length(udelta) < p) {
	b[idx.bz] <- 0
	dd <- outer(delta, delta[idx.bz], "-")
	test <- abs(dd[-idx.bz,]) < eps
	cand <- which(colSums(test) == 0)
	if (length(cand) > 0) {
		z0 <- b / dd[,cand]
		nrm <- sqrt(colSums(z0^2))
		z0 <- t(t(z0) * (c/nrm))
		objective <- 
	}

}

	
	
}