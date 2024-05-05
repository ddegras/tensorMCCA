mcca.init.svd <- function(x, w = NULL, objective = c("cov", "cor"), 
	scope = c("block", "global"))
{
	
################
# Preprocessing
################

test <- check.arguments(x, w = w)
objective <- match.arg(objective) 
scope <- match.arg(scope)
eps <- 1e-14

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dimfun) # full data dimensions
n <- tail(dimx[[1]], 1)
p <- lapply(dimx, function(idx) idx[-length(idx)])
d <- sapply(d, length)
pp <- sapply(p, prod)
cumpp <- cumsum(c(0,pp))

## Objective weights
if (is.null(w)) {
	w <- 1 - diag(m)
} else if (length(w) == 1) {
	w <- matrix(1, m, m)
}
if (objective == "cor") diag(w) <- 0
w <- (w + t(w)) / (2 * sum(w)) 

## Check for constant datasets
cnst.set <- sapply(x, function(a) all(a == a[1]))
w[,cnst.set] <- w[cnst.set,] <- 0
wzero <- (colSums(w = 0) == m)
wnzero <- !wzero
cnstfun <- switch(objective, 
	cov = function(len) numeric(len),
	cor = function(len) rep(1,len))

## Trivial case: all datasets are constant or associated weights are zero 
if (all(wzero)) {
	v <- vector("list", m)
	for (i in 1:m)
		v[[i]] <- lapply(p[[i]], cnstfun)
	if (objective == "cor")
		v <- scale.v(v, type = "var", x = x, check.args = FALSE)
	return(v)
}

## Data means for centering
xbar <- vector("list", m)
uncentered <- logical(m)
for(i in 1:m) {
    xbar[[i]] <- as.vector(rowMeans(x[[i]], dims = d[i]))
	uncentered[i] <- any(abs(xbar[[i]]) > eps)
}

## Trivial case: all but one dataset are constant 
## or have associated weights equal to zero 
if (sum(wnzero) == 1) {
	v <- vector("list", m)
	for (i in which(wzero)) 
		v[[i]] <- lapply(p[[i]], numeric)
	i <- which(wnzero)
	v[[i]] <- tnsr.rk1(scale = TRUE,
		x = if (uncentered[i]) {x[[i]] - xbar[[i]]} else x[[i]])
		v[[i]] <- scale.v(v[[i]])
	return(v)
}

## Remove any constant dataset
if (any(wzero)) {
	vfull <- vector("list", m)
	for (i in which(wzero)) 
		vfull[[i]] <- lapply(p[[i]], cnstfun)
	if (objective == "cor")
		vfull[wzero] <- scale.v(vfull[wzero], "var", x = x[wzero],
			check.args = FALSE)
	x <- x[wnzero]
	d <- d[wnzero]
	m <- sum(wnzero)
	p <- p[wnzero]
	w <- w[wnzero, wnzero]
}
pp <- sapply(p, prod)
k <- if (is.null(k)) 
	pmin(pp, n) else pmin(rep_len(k, m), pp, n)



#######################################
# Perform SVD/EVD of concatenated data 
#######################################

## Canonical weights 
v <- vector("list", m)

## Unfold data along last mode (individuals/objects)
## and transform matricized data by SVD if needed 
test <- (objective == "cor" || all(pp > n))
if (test) svdx <- vector("list", m)
xmat <- x
for (i in 1:m) {
	dim(xmat[[i]]) <- c(pp[i], n)
	xbar <- rowMeans(xmat[[i]])
	if (any(abs(xbar) > eps))
		xmat[[i]] <- xmat[[i]] - xbar	
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
nrowx <- sapply(xmat, nrow) # numbers of rows in transformed datasets
cumnrowx <- cumsum(c(0,nrowx))

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
	covx <- matrix(0, cumnrowx[m+1], cumnrowx[m+1])
	for (i in 1:m) {
		idxi <- (cumnrowx[i]+1):cumnrowx[i+1]
		for (j in 1:i) {
			if (w[i,j] == 0) next
			idxj <- (cumnrowx[j]+1):cumnrowx[j+1]
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
	phi <- numeric(cumpp[m+1])
	for (i in 1:m) {
		idxfull <- (cumpp[i]+1):cumpp[i+1]
		idxsmall <- (cumnrowx[i]+1):cumnrowx[i+1]
		phi[idxfull] <- switch(objective,
			cov = svdx[[i]]$u %*% phitmp[idxsmall],
			cor = svdx[[i]]$u %*% (phitmp[idxsmall] / svdx[[i]]$d)
		)
	}
	rm(svdx)	
}



#######################################################
# Split eigenvector according to dataset dimensions, 
# tensorize, and perform rank-1 approximation by HOSVD
#######################################################

covblock <- (objective == "cov" && scope == "block")
for (i in 1:m) {
	idxi <- (cumpp[i]+1):cumpp[i+1]
	phii <- array(phi[idxi], p[[i]])
	svdi <- hosvd(phii, 1)
	v[[i]] <- lapply(svdi$factors, as.numeric) 
	if (!covblock)
		v[[i]][[1]] <- v[[i]][[1]] * as.numeric(svdi$core)
}

## Scale canonical weights and reorient as needed
if (!covblock) {
	v <- if (objective == "cov") {
		scale.v(v, scope = scope, check.args = FALSE)
	} else scale.v(v, "var", scope, x, FALSE)
}
score <- canon.scores(x, v, FALSE)
if (center) score <- score - matrix(colMeans(score), n, m, TRUE)
flip <- reorient(score, w)$flip
for (i in which(flip))
	v[[i]][[1]] <- (-v[[i]][[1]])

## Put back any canonical weights for constant datasets
if (any(wzero)) {
	vfull[wnzero] <- v
	v <- vfull
	if (objective == "cov" && scope == "global")
		v <- scale.v(v, "norm", "global")
}

v
}

