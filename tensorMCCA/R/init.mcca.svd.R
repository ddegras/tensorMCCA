




###########################
# SVD-based initialization
# of canonical vectors 
# (Quick Rank 1)
###########################


init.mcca.svd <- function(x, objective = c("cov", "cor"), 
	cnstr = c("block", "global"), center = TRUE)
{
## Check argument x if required
test <- check.arguments(x)
tol <- 1e-14

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
for (i in 1:m) 
	dim(x[[i]]) <- c(prod(p[[i]]), n)
x <- do.call(rbind, x)

## Data centering
if (center) 
	x <- x - rowMeans(x)

## Initial SVD
svdfun <- if (max(dim(x)) > 2) {
	function(x) svds(x, k = 1) } else {
	function(x) svd(x, nu = 1, nv = 1) } 
x <- svdfun(x)$u

## Iteratively unfold singular vectors and recalculate SVD
lenx <- sapply(p, prod)
start <- c(0, cumsum(lenx[-m])) + 1
end <- cumsum(lenx)
for (i in 1:m) {
	pi <- p[[i]]
	xi <- x[start[i]:end[i]]
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
if (objective == "cor")
	v <- scale.v(v, x,"var", cnstr)

v
}

