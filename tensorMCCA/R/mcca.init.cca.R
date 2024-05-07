mcca.init.cca <- function(x, w = NULL, objective = c("cov", "cor"), 
	scope = c("block", "global"), k = NULL, optim = NULL, 
	maxit = 100, tol = 1e-4)
{
	
################
# Preprocessing
################


## Check arguments x and w 
test <- check.arguments(x, w = w)
eps <- 1e-14

## Scaling constraints
objective <- match.arg(objective)   
scope <- match.arg(scope) # block or global constraints
if (objective == "cor" && scope == "global")
	stop(paste("Argument values 'objective = cor'", 
		"and 'scope = global' are incompatible."))

## Search method in combinatorial optimization
optim <- if (is.null(optim)) {
	ifelse(m <= 5, "exact", "greedy")
} else {
	match.arg(optim, c("exact", "greedy")) 
}

## Data dimensions
m <- length(x)
dimx <- lapply(x, dimfun)
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

## Identify constant datasets and set their weights to zero
constant <- sapply(x, function(a) all(a == a[1]))
w[, constant] <- w[constant,] <- 0
wzero <- apply(w == 0, 2, all)
wnzero <- which(!wzero)

## Trivial case: all weights are zero 
if (all(wzero)) {
	v <- relist(lapply(unlist(p), numeric), p)
	return(v)
}

## Remove any constant dataset
if (any(wzero)) {
	x <- x[wnzero] 
	w <- w[wnzero, wnzero]
	vfull <- vector("list", m)
	for (i in which(wzero)) 
		vfull[[i]] <- lapply(p[[i]], numeric)
	d <- d[wnzero]
	k <- k[wnzero]
	m <- length(wnzero)
	p <- p[wnzero]
	pp <- pp[wnzero]
}

## Calculate data means
xbar <- vector("list", m)
uncentered <- logical(m)
for(i in 1:m) {
    xbar[[i]] <- if (d[i] == 0L) { 
    		mean(x[[i]])
    	} else {
	    	as.vector(rowMeans(x[[i]], dims = d[i]))
	}
	uncentered[i] <- any(abs(xbar[[i]]) > 1e-16)
}





###########################
# SVD of unfolded datasets
###########################


## Calculate truncated SVD of each unfolded dataset
## if necessary or computationally efficient
ux <- vx <- vector("list", m)
reduce <- (k < pmin(pp,n)) | (k <= pp/2) | 
	(objective == "cor") | (diag(w) > 0)
for (i in 1:m) {
	if (!reduce[i]) next	
	xmat <- x[[i]] 
	dim(xmat) <- c(pp[i], n)
	if (uncentered[i]) xmat <- xmat - xbar[[i]] 
	svdx <- tryCatch(suppressWarnings(svds(xmat, k[i])), 
		error = function(e) svd(xmat, k[i], k[i]))
	pos <- (svdx$d > eps)
	if (!all(pos)) {
		k[i] <- sum(pos)
		svdx$u <- svdx$u[,pos,drop=FALSE]
		svdx$d <- svdx$d[pos]
		svdx$v <- svdx$v[,pos,drop=FALSE]
	}
	if (objective == "cov") {
		ux[[i]] <- svdx$u
		vx[[i]] <- sweep(svdx$v, 2, svdx$d, "*")
	} else if (objective == "cor") {
		ux[[i]] <- sweep(svdx$u, 2, svdx$d, "/")
		vx[[i]] <- svdx$v		
	}
}
rm(xmat, svdx)



##########################################
# Calculate CCA for each pair of datasets 	
##########################################


## Canonical vectors
a <- lapply(pp, function(nr) matrix(0, nr, m)) 
# a[[i]][,j] is the first left canonical vector associated to Xi Xj'
# For j > i, a[[j]][,i] is calculated as the first right canonical vector 
# associated to Xi Xj' 
 
for (i in 1:m) {	
	if (reduce[i]) {
		xi <- vx[[i]]
	} else {
		xi <- x[[i]]
		dim(xi) <- c(pp[i],n)
		if (uncentered[i]) 
			xi <- xi - xbar[[i]]
	}	
	for (j in 1:i) {
		if (w[i,j] == 0) next
		if (i == j) {
			if (objective == "cov")
				a[[i]][,i] <- ux[[i]][,1] 
			next
		}
		if (reduce[j]) {
			xj <- vx[[j]]
		} else {
			xj <- x[[j]]
			dim(xj) <- c(pp[j],n)
			if (uncentered[j]) 
				xj <- xj - xbar[[j]]
		}
		svdij <- tryCatch(
			supressWarnings(svds(crossprod(xi, xj), k = 1)), 
			error = function(e) svd(crossprod(xi, xj), 1, 1))
		a[[i]][,j] <- if (reduce[i]) {
			ux[[i]] %*% svdij$u } else svdij$u	
		a[[j]][,i] <- if (reduce[j]) {
			ux[[j]] %*% svdij$v } else svdij$v
	}
}

rm(ux, vx, xi, xj, svdij)



#####################################
# Approximate long canonical vectors 
# by rank-1 tensors 
#####################################


v <- vector("list", m^2)
dim(v) <- c(m, m)
covblock <- (objective == "cov" && scope == "block")
for (i in 1:m) {
	for (j in 1:m) {
		aij <- a[[i]][,j]
		if (d[i] > 1) dim(aij) <- p[[i]]
		v[[i,j]] <- tnsr.rk1(aij, scale = covblock, 
			maxit = maxit, tol = tol)
		# v[[i,j]] <- switch(objective,
			# cov = tnsr.rk1(aij, scale = TRUE, maxit = maxit, tol = tol),
			# cor = tnsr.rk1.score(aij, cnstr = x[[i]] - xbar[[i]], 
				# maxit = maxit, tol = tol))
	}
}

## Scale tensors according to norm (objective = covariance)
## or variance (objective = correlation)
if (objective == "cov" && scope == "global") {
	v <- scale.v(v, check.args = FALSE)
} else if (objective == "cor") {
	v <- scale.v(v, type = "var", x = x, check.args = FALSE)
}




#############################
# Calculate canonical scores 
#############################


score <- canon.scores(x, v)
score <- sweep(score, 2:3, colMeans(score))



#####################################
# Find best combination of canonical 
# vectors across datasets 
#####################################


opt <- switch(optim, 
	exact = optim.combn.exact(score, w),
	greedy = optim.combn.greedy(score, w))
v <- v[cbind(1:m, opt$idx)]
flip <- which(opt$sign == -1)
for (i in flip) 
	v[[i]][[1]] <- -v[[i]][[1]] 



####################################
# Case: maximize sum of covariances
# with global norm constraint
####################################


if (objective == "cov" && scope == "global") {
	score <- canon.scores(x, v)
	s <- sqrt(m) * eigen(w * cov(score), TRUE)$vectors[,1]
	sgns <- sign(s)
	for (i in 1:m) {
		if (d[i] <= 1) {
			si <- s[i]
		} else {
			si <- rep(abs(s[i])^(1/d[i]), d[i])
			si[1] <- si[1] * sgns[i]
		}
		v[[i]] <- mapply("*", x = v[[i]], y = si, SIMPLIFY = FALSE)
	}
}



#################
# Postprocessing
#################


## Drop singleton dimensions
for (i in 1:m)
	v[[i]] <- lapply(v[[i]], drop)

## Put back any canonical weights for constant datasets
if (any(wzero)) {
	vfull[wnzero] <- v
	v <- vfull
}


v

}



