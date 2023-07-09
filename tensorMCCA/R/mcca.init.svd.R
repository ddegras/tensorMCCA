mcca.init.svd <- function(x, w = NULL, objective = c("cov", "cor"), 
	scale = c("block", "global"), center = TRUE)
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
scale <- match.arg(scale)
if (is.null(w)) w <- matrix(1, m, m)

## Initialize canonical weights 
v <- vector("list", m)

## Unfold data along last mode (individuals/objects) and concatenate
xmat <- x
for (i in 1:m) 
	dim(xmat[[i]]) <- c(prod(p[[i]]), n)
xmat <- do.call(rbind, xmat)

## Data centering
if (center) 
	xmat <- xmat - rowMeans(xmat)

## Catch trivial case
if (all(xmat == 0)) {
	for (i in 1:m) 
		v[[i]] <- lapply(p[[i]], 
			function(len) rep(1/sqrt(len), len))
	return(v)
}

## Initial SVD
dim.flag <- all(dim(xmat) > 2)
error.flag <- FALSE
if (dim.flag) {
	svdx <- tryCatch(svds(xmat, k = 1, nv = 0)$u, 
		error = function(e) NULL)
	error.flag <- is.null(svdx)
	if (!error.flag) xmat <- svdx 
} 
if ((!dim.flag) || error.flag)
	xmat <- svd(xmat, nu = 1, nv = 0)$u
dim(xmat) <- NULL


## Iteratively unfold singular vectors and recalculate SVD
lenx <- sapply(p, prod)
end <- cumsum(lenx)
start <- c(0, end[-m]) + 1
for (i in 1:m) {
	xi <- xmat[start[i]:end[i]]
	dim(xi) <- if (d[i] == 1) NULL else p[[i]]
	v[[i]] <- if (all(xi == 0)) {
		lapply(p[[i]], function(len) rep(1/sqrt(len), len))
	} else { lapply(hosvd(xi, 1)$factors, as.vector) }
}


## Scale initial canonical weights as needed	
if (objective == "cov" && scale == "global") {
	y <- canon.scores(x, v)
	y <- scale(y, scale = FALSE)
	M <- crossprod(y) * (w/n)
	a <- eigen(M, TRUE)$vectors[,1] * sqrt(m) 
	s <- a^(1/d) # scaling factors
	for (i in 1:m)
		v[[i]] <- lapply(v[[i]], "*", y = s[i])
} else if (objective == "cor") {
	y <- canon.scores(x, v)
	y <- scale(y, scale = FALSE)
	nrm <- sqrt(colMeans(y^2))
	for (i in 1:m) {
		v[[i]] <- if (nrm[i] > eps) {
			lapply(v[[i]], "/", y = nrm[i]^(1/d[i]))
		} else { lapply(p[[i]], numeric) }
	}
}

v
}

