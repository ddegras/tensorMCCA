mcca.init.svd <- function(x, w = NULL, objective = c("cov", "cor"), 
	scale = c("block", "global"), center = TRUE)
{
test <- check.arguments(x, w)
objective <- match.arg(objective) 
scale <- match.arg(scale)
eps <- 1e-14
w <- if (is.null(w)) {
	matrix(1/m^2, m, m) 
} else {
	(w + t(w)) / (2 * sum(w))
}

## Data dimensions
m <- length(x) # number of datasets 
dimfun <- function(x) if (is.vector(x)) c(1,length(x)) else dim(x)
dimx <- lapply(x, dimfun) # full data dimensions
n <- tail(dimx[[1]], 1)
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
pp <- sapply(p, prod)
cpp <- cumsum(c(0,pp))


## Canonical weights 
v <- vector("list", m)

## Unfold data along last mode (individuals/objects)
## and transform matricized data by SVD if needed 
test <- (objective == "cor" || all(pp > n))
if (test) {
	svdx <- vector("list", m)
	rankx <- numeric(m)
}
xmat <- x
for (i in 1:m) {
	dim(xmat[[i]]) <- c(pp[i], n)
	if (center)
		xmat[[i]] <- xmat[[i]] - rowMeans(xmat[[i]])	
	if (test) {
		svdx[[i]] <- svd(xmat[[i]])
		pos <- (svdx[[i]]$d > max(svdx[[i]]$d[1] * 1e-8, eps))
		svdx[[i]]$u <- svdx[[i]]$u[,pos,drop=FALSE]
		svdx[[i]]$d <- svdx[[i]]$d[pos]
		svdx[[i]]$v <- svdx[[i]]$v[,pos,drop=FALSE]
		xmat[[i]] <- switch(objective, 
			cor = t(svdx[[i]]$v)	, 
			cov = t(svdx[[i]]$v) * svdx[[i]]$d)
	}	
}
nrx <- sapply(xmat, nrow) # numbers of rows in transformed datasets
cnrx <- cumsum(c(0,nrx))


## Check separability of objective weights
if (qr(w)$rank == 1) {
	# Separable case: concatenate data + SVD
	a <- colSums(w) # w = aa'
	if (!all(a == a[1])) {
		for (i in 1:m) 
			xmat[[i]] <- a[i] * xmat[[i]]
	}
	xmat <- do.call(rbind, xmat)
	phi <- tryCatch(svds(xmat, k = 1, nv = 0)$u, 
		error = function(e) svd(xmat, 1, 0)$u)
	rm(xmat)
} else {
	# Nonseparable case: form pseudo-covariance matrix + EVD
	covx <- matrix(0, cnrx[m+1], cnrx[m+1])
	for (i in 1:m) {
		idxi <- (cnrx[i]+1): cnrx[i+1]
		for (j in 1:i) {
			if (w[i,j] == 0) next
			idxj <- (cnrx[j]+1):cnrx[j+1]
			covx[idxi,idxj] <- w[i,j] * tcrossprod(xmat[[i]], xmat[[j]])
			if (j < i) covx[idxj,idxi] <- t(covx[idxi,idxj])
		}
	}
	rm(xmat)
	phi <- tryCatch(eigs_sym(covx, k = 1, which = "LA")$vectors, 
		error = function(e) eigen(covx, TRUE)$vectors[,1])
	rm(covx)
}

## Transform back singular/eigen-vector in original space if needed 
if (test) {
	phitmp <- phi
	phi <- numeric(cpp[m+1])
	for (i in 1:m) {
		idxfull <- (cpp[i]+1):cpp[i+1]
		idxsmall <- (cnrx[i]+1):cnrx[i+1]
		phi[idxfull] <- switch(objective,
			cov = svdx[[i]]$u %*% phitmp[idxsmall],
			cor = svdx[[i]]$u %*% (phitmp[idxsmall] / svdx[[i]]$d)
		)
	}
	rm(svdx)	
}

## Split eigenvector according to dataset dimensions, tensorize, 
## and perform rank-1 approximation by HOSVD
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
if (center) score <- score - matrix(colMeans(score), n, m, TRUE)
flip <- reorient(score, w)$flip
for (i in which(flip))
	v[[i]][[1]] <- (-v[[i]][[1]])



v
}

