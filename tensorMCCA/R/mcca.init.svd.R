mcca.init.svd <- function(x, w = NULL, objective = c("cov", "cor"), 
	scale = c("block", "global"), center = TRUE)
{
## Check argument x if required
test <- check.arguments(x, w)
eps <- 1e-14

## Data dimensions
m <- length(x) # number of datasets 
dimfun <- function(x) if (is.vector(x)) c(1,length(x)) else dim(x)
dimx <- lapply(x, dimfun) # full data dimensions
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
pp <- sapply(p, prod)
d <- sapply(p, length)
n <- tail(dimx[[1]], 1)

objective <- match.arg(objective) # norm or variance constraints
scale <- match.arg(scale)
w <- if (is.null(w)) {
	matrix(1/m^2, m, m) 
} else {
	(w + t(w)) / (2 * sum(w))
}

## Initialize canonical weights 
v <- vector("list", m)

## Unfold data along last mode (individuals/objects) 
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
		pos <- (svdx[[i]]$d > max(svdx[[i]]$d[1] * 1e-8, 1e-14))
		rankx[i] <- sum(pos)
		svdx[[i]]$u <- svdx[[i]]$u[,pos]
		svdx[[i]]$d <- svdx[[i]]$d[pos]
		svdx[[i]]$v <- svdx[[i]]$v[,pos]
		xmat[[i]] <- switch(objective, 
			cor = t(svdx[[i]]$v)	, 
			cov = t(svdx[[i]]$v) * svdx[[i]]$d)
	}	
}

#  max a'x y'b    a'xx'a = 1.  
# if p > n, x = QR  , set a = Q (R')^- aa.  aa' R^- Q' QR 
# if n > p, x' = QR, a = R^- aa

# x = UDV' a = U D^- aa
# x' = QR. x = R^- z

## Check separability of objective weights
if (qr(w)$rank == 1) {
	a <- colSums(w) # w = aa'
	if (!all(a == a[1])) {
		for (i in 1:m) 
			xmat[[i]] <- xmat[[i]] * a[i]
	}
	xmat <- do.call(rbind, xmat)
	phi <- tryCatch(svds(xmat, k = 1, nv = 0)$u, 
		error = function(e) svd(xmat, 1, 0)$u)
	rm(xmat)
} else {
	## Calculate pseudo-covariance matrix
	nr <- sapply(xmat, nrow)
	cumnr <- cumsum(c(0, nr))
	covx <- matrix(0, cumnr[m+1], cumnr[m+1])
	for (i in 1:m) {
		idxi <- (cumnr[i]+1): cumnr[i+1]
		for (j in 1:i) {
			if (w[i,j] == 0) next
			idxj <- (cumnr[j]+1):cumnr[j+1]
			covx[idxi,idxj] <- w[i,j] * tcrossprod(xmat[[i]], xmat[[j]])
			if (j < i) covx[idxj,idxi] <- t(covx[idxi,idxj])
		}
	}
	rm(xmat)
	## EVD
	phi <- tryCatch(eigs_sym(covx, k = 1, which = "LA")$vectors, 
			error = function(e) eigen(covx, TRUE)$vectors[,1])
	rm(covx)
}

if (objective == "cor") {
	cpp <- cumsum(c(0, pp))
	phil <- vector("list", m)
	for (i in 1:m) {
		idxi <- (cpp[i]+1):cpp[i+1]
		philid
		phi[idxi] <- 
	}		
}	




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

