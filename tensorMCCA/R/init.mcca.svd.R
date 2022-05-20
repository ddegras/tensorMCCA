




#############################
# HOSVD-based initialization
# of canonical vectors 
# (Quick Rank 1)
#############################


init.mcca.svd <- function(x, objective = c("covariance", "correlation"), 
	cnstr = c("block", "global"), center = TRUE)
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
cnstr <- match.arg(cnstr) # block or global constraints

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
svdfun <- if (max(dim(xmat)) > 2) {
	function(x) svds(x, k = 1) } else {
	function(x) svd(x, nu = 1, nv = 1) } 
xmat <- svdfun(xmat)$u

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
			svdx <- svdfun(xi)
			v[[i]][[k]] <- svdx$u
			xi <- svdx$v
		} else {
			v[[i]][[k]] <- xi 
		}	
	}
}

## Scale initial canonical vectors as required	
if (objective == "correlation") {
	scores <- image.scores(x, v)
	s <- sqrt(colMeans(scores^2))
	for (i in 1:m) 
	for (k in 1:d[i]) 
		v[[i]][[k]] <- if (s[i] <= eps) {
			numeric(p[[i]][k]) } else {
				v[[i]][[k]] / s[i]^(1/d[i]) }
}

v
}

