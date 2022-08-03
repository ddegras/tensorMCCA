#############################
# HOSVD-based initialization
# of canonical vectors 
# (Quick Rank 1)
#############################

# Development version with multiple starting points

mcca.init.svd <- function(x, r = 1L, objective = c("cov", "cor"), 
	ortho = c("score", "canon.tnsr"), center = TRUE, check.args = TRUE,
	reduce.storage = FALSE)
{
if (check.args) 
	test <- check.arguments(x)
eps <- 1e-14
r <- as.integer(r)

## Data dimensions
m <- length(x) 
dimx <- lapply(x, dim) 
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
d <- sapply(p, length)
n <- tail(dimx[[1]], 1)

objective <- match.arg(objective) # norm or variance constraints
norm <- match.arg(norm) # block or global norm constraints

## Initialize canonical vectors 
r <- min(n, r)
v <- vector("list", m * r)
dim(v) <- c(m, r)

## Initial SVD in reduced space (incremental)
if (low.storage) {
	xx <- matrix(0, n, n)
	for (i in 1:m) {
		xi <- x[[i]]
		dim(xi) <- c(prod(p[[i]]), n)
		if (center) xi <- xi - rowMeans(xi)
		xx <- xx + crossprod(xi)
	}
	
}
## Unfold data along last mode (individuals/objects) and concatenate
xmat <- x
for (i in 1:m) 
	dim(xmat[[i]]) <- c(prod(p[[i]]), n)
xmat <- do.call(rbind, xmat)

## Data centering
if (center) 
	xmat <- xmat - rowMeans(xmat)

## Initial SVD
xmat <- if (all(dim(xmat) > 2)) {
	svds(xmat, k = 1, nv = 0)$u 
} else svd(xmat, nu = 1, nv = 0)$u


## Iteratively unfold singular vectors and recalculate SVD
lenx <- sapply(p, prod)
start <- c(0, cumsum(lenx[-m])) + 1
end <- cumsum(lenx)
for (i in 1:m) {
	pi <- p[[i]]
	xi <- xmat[start[i]:end[i]]
	for (k in 1:d[i]) {	
		if (k < d[i]) {
			dim(xi) <- c(pi[k], length(xi) / pi[k])
			svdx <- if (all(dim(xi) > 2)) {
				svds(xi, k = 1)
			} else svd(xi, nu = 1, nv = 1)
			v[[i]][[k]] <- svdx$u
			xi <- svdx$v
		} else {
			v[[i]][[k]] <- xi 
		}	
	}
}

## Scale initial canonical vectors as required	
if (objective == "cor") {
	y <- canon.scores(x, v)
	nrm <- sqrt(colMeans(y^2))
	nz <- which(nrm > eps)
	for (i in 1:m) {
		v[[i]] <- if (nz[i]) {
			lapply(v[[i]], "/", y = nrm[i]^(1/d[i]))
		} else { lapply(p[[i]], numeric) }
	}
}

v
}

