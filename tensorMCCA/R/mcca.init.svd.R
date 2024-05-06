mcca.init.svd <- function(x, w = NULL, objective = c("cov", "cor"), 
	scope = c("block", "global"), k = NULL)
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
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
p[d == 0] <- 1
n <- tail(dimx[[1]], 1)

## Truncation order in SVD
pp <- sapply(p, prod)
if (!is.null(k)) {
	k <- as.integer(k)
	stopifnot(all(k > 0))
	k <- rep_len(k, m)
	k <- pmin(k, pp, n)
} else {
	k <- pmin(pp, n)
}

## Objective weights
if (is.null(w)) {
	w <- 1 - diag(m)
} else if (length(w) == 1) {
	w <- matrix(1, m, m)
}
w <- (w + t(w)) / (2 * sum(w)) 

## Check for constant datasets
constant <- sapply(x, function(a) all(a == a[1]))
w[, constant] <- w[constant,] <- 0
wzero <- apply(w == 0, 2, all)
wnzero <- which(!wzero)

## Trivial case: all datasets have associated weights zero 
if (all(wzero)) {
	v <- relist(lapply(unlist(p), numeric), p)
	return(v)
}

## Data means for centering
xbar <- vector("list", m)
uncentered <- logical(m)
for (i in 1:m) {
	xbar[[i]] <- if (d[i] == 0L) { 
		mean(x[[ii]])
    } else {
		as.vector(rowMeans(x[[ii]], dims = d[i]))
	}
	uncentered[i] <- any(abs(xbar[[i]]) > eps)
}

## Remove any constant dataset
if (any(wzero)) {
	vfull <- vector("list", m)
	for (i in which(wzero)) 
		vfull[[i]] <- lapply(p[[i]], numeric)
	d <- d[wnzero]
	k <- k[wnzero]
	m <- length(wnzero)
	p <- p[wnzero]
	pp <- pp[wnzero]
	w <- w[wnzero, wnzero]
}
cumpp <- c(0, cumsum(pp))



#######################################
# Perform SVD/EVD of concatenated data 
#######################################

## Canonical weights 
v <- vector("list", m)

## Unfold data along last mode (individuals/objects)
## and transform matricized data by SVD if needed 
reduce <- (k < pmin(pp,n)) | (k <= pp/2) | (objective == "cor")
u <- xmat <- vector("list", m)
for (i in 1:m) {
	ii <- wnzero[i]
	xmat[[i]] <- x[[ii]]
	dim(xmat[[i]]) <- c(pp[i], n)
	if (uncentered[ii])
		xmat[[i]] <- xmat[[i]] - xbar[[ii]]	
	if (reduce[i]) {
		svdx <- svd(xmat[[i]])
		pos <- (svdx$d > max(1e-8 * svdx$d[1], 1e-14))
		u[[i]] <- svdx$u[,pos,drop=FALSE]
		svdx$d <- svdx$d[pos]
		svdx$v <- svdx$v[,pos,drop=FALSE]
		if (objective == "cov") {
			u[[i]] <- svdx$u
			xmat[[i]] <- t(svdx$v) * svdx$d			
		} else if (objective == "cor") {
			u[[i]] <- sweep(u[[i]], 2, svdx$d, "/")
			xmat[[i]] <- t(svdx$v)
		} 
	}	
}
suppressWarnings(rm(svdx))
nrowx <- sapply(xmat, NROW) # numbers of rows in transformed datasets
cumnrowx <- cumsum(c(0,nrowx))

## Check separability of objective weights
sepw <- separable(w, objective == "cor")

# Separable case: concatenate data + SVD
if (sepw$separable) {
	a <- sepw$a # w = aa'
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

## Map back singular/eigen-vector to original space if needed 
if (any(reduce)) {
	phi <- split(phi, rep(1:m, nrowx))
	for (i in which(reduce)) 
		phi[[i]] <- u[[i]] %*% phi[[i]]
	phi <- unlist(phi)
}



#######################################################
# Split eigenvector according to dataset dimensions, 
# tensorize, and perform rank-1 approximation by HOSVD
#######################################################

covblock <- (objective == "cov" && scope == "block")
for (i in 1:m) {
	idxi <- (cumpp[i]+1):cumpp[i+1]
	phii <- array(phi[idxi], p[[i]])
	v[[i]] <- tnsr.rk1(phii, scale = covblock, 
		maxit = maxit, tol = tol)
}

## Scale canonical weights 
if (objective == "cor") {
	v <- scale.v(v, type = "var", x = x[wnzero], check.args = FALSE)
} else if (objective == "cov" && scope == "global") {
	v <-  scale.v(v, scope = "global", check.args = FALSE)
}

## Calculate scores
score <- canon.scores(x[wnzero], v)
if (any(uncentered[wnzero])) 
	score <- sweep(score, 2L, colMeans(score), check.margin = FALSE)

## Find best rescaling or orientation for weights
if (objective == "cov" && scope == "global") {
	s <- eigen(w * cov(score) * w, TRUE)$vectors[,1]
	for (i in 1:m) v[[i]][[1]] <- s[i] * v[[i]][[1]]
} else {
	flip <- reorient(score, w)$flip
	for (i in which(flip))
		v[[i]][[1]] <- (-v[[i]][[1]])
}

## Put back any canonical weights for constant datasets
if (any(wzero)) {
	vfull[wnzero] <- v
	v <- vfull
	if (objective == "cov" && scope == "global")
		v <- scale.v(v, "norm", "global")
}

v
}

