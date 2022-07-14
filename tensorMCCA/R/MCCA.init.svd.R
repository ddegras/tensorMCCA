#############################
# HOSVD-based initialization
# of canonical vectors 
# (Quick Rank 1)
#############################


mcca.init.svd <- function(x, objective = c("cov", "cor"), 
	norm = c("block", "global"), center = TRUE)
{
## Check argument x if required
test <- check.arguments(x)
eps <- 1e-14

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
d <- sapply(p, length)
n <- tail(dimx[[1]], 1)

objective <- match.arg(objective) # norm or variance constraints
norm <- match.arg(norm) # block or global norm constraints

## Initialize canonical vectors 
v <- vector("list", m)

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

