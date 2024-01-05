mcca.init.svd <- function(x, w = NULL, objective = c("cov", "cor"), 
	scale = c("block", "global"), center = TRUE)
{
## Check argument x if required
test <- check.arguments(x)
eps <- 1e-14

## Data dimensions
m <- length(x) # number of datasets 
dimfun <- function(x) if (is.vector(x)) c(1,length(x)) else dim(x)
dimx <- lapply(x, dimfun) # full data dimensions
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
d <- sapply(p, length)
n <- tail(dimx[[1]], 1)

objective <- match.arg(objective) # norm or variance constraints
scale <- match.arg(scale)
if (is.null(w)) w <- matrix(1, m, m)

## Initialize canonical weights 
v <- vector("list", m)

## Unfold data along last mode (individuals/objects) 
xmat <- x
for (i in 1:m) {
	dim(xmat[[i]]) <- c(prod(p[[i]]), n)
	if (center) 
		xmat[[i]] <- xmat[[i]] - rowMeans(xmat[[i]])	
}

## Calculate pseudo-covariance covariance matrix
pp <- sapply(p, prod)
cpp <- cumsum(c(0, pp))
covx <- matrix(0, cpp[m+1], cpp[m+1])
for (i in 1:m) {
	idxi <- (cpp[i]+1):cpp[i+1]
	for (j in 1:i) {
		if (w[i,j] == 0) next
		idxj <- (cpp[j]+1):cpp[j+1]
		covx[idxi,idxj] <- w[i,j] * tcrossprod(xmat[[i]], xmat[[j]])
		if (j < i) covx[idxj,idxi] <- t(covx[idxi,idxj])
	}
}
rm(xmat)

## Initial EVD
phi <- tryCatch(eigs_sym(covx, k = 1, which = "LA")$vectors, 
		error = function(e) NULL)
if (is.null(phi))
	phi <- eigen(covx, TRUE)$vectors[,1]
rm(covx)

## Split eigenvector according to dataset dimensions, tensorize, 
## and perform HOSVD
covblock <- (objective == "cov" && scale == "block")
for (i in 1:m) {
	idxi <- (cpp[i]+1):cpp[i+1]
	phii <- array(phi[idxi], p[[i]])
	svdi <- hosvd(phii, 1)
	v[[i]] <- lapply(svdi$factors, as.numeric) 
	if (!covblock)
		v[[i]][[1]] <- v[[i]][[1]] * as.numeric(svdi$core)
}


## Scale canonical weights and reorient as needed
if (!covblock) {
	v <- if (objective == "cov") {
		scale.v(v, scale = scale, check.args = FALSE)
	} else scale.v(v, "var", scale, x, FALSE)
}
score <- canon.scores(x, v, FALSE)
flip <- reorient(score, w)$flip
for (i in which(flip))
	v[[i]][[1]] <- (-v[[i]][[1]])



v
}

