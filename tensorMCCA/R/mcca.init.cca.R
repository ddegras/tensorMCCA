mcca.init.cca <- function(x, w = NULL, objective = c("cov", "cor"), 
	scope = c("block", "global"), k = NULL, optim = NULL, 
	maxit = 100, tol = 1e-4)
{
	
################
# Preprocessing
################


m <- length(x) # number of datasets 
for (i in 1:m) { # Matricize any vector dataset
	if (is.vector(x[[i]])) 
		dim(x[[i]]) <- c(1, length(x[[i]]))	
}

## Check arguments x and w 
test <- check.arguments(x, w = w)
eps <- 1e-14

## Scaling constraints
objective <- match.arg(objective)   
scope <- match.arg(scope) # block or global constraints
if (objective == "cor" && scope == "global")
	stop(paste("Argument values 'objective = cor'", 
		"and 'scope = global' are incompatible."))

## Data dimensions
m <- length(x)
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
pp <- sapply(p, prod)
n <- tail(dimx[[1]], 1)

## Reshape and center data
for (i in 1:m) {
	dim(x[[i]]) <- c(pp[i], n)
	xbar <- rowMeans(x[[i]])
	if (any(abs(xbar) > 1e-16)) 
		x[[i]] <- sweep(x[[i]], 2, xbar, "-")
}

## Search method in combinatorial optimization
optim <- if (is.null(optim)) {
	ifelse(m <= 5, "exact", "greedy")
} else {
	match.arg(optim, c("exact", "greedy")) 
}

## Truncation order in SVD
if (is.null(k)) {
	k <- pmin(pp, n)
} else {	
	k <- as.integer(k)
	stopifnot(all(k > 0))
	k <- rep_len(k, m)
	k <- pmin(k, pp, n)
}
		
## Objective weights
if (is.null(w)) {
	w <- 1 - diag(m)
} else if (length(w) == 1) {
	w <- matrix(1, m, m)
} else {
	w <- (w + t(w)) / 2  
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
	svdx <- tryCatch(suppressWarnings(svds(x[[i]], k[i])), 
		error = function(e) svd(x[[i]], k[i], k[i]))
	pos <- (svdx$d > eps)
	if (!all(pos)) {
		k[i] <- sum(pos)
		svdx$u <- svdx$u[,pos,drop=FALSE]
		svdx$d <- svdx$d[pos]
		svdx$v <- svdx$v[,pos,drop=FALSE]
	}
	if (objective == "cov") {
		ux[[i]] <- svdx$u
		vx[[i]] <- t(svdx$v) * svdx$d
	} else if (objective == "cor") {
		ux[[i]] <- sweep(svdx$u, 2, svdx$d, "/")
		vx[[i]] <- t(svdx$v)		
	}
}
if (any(reduce))	
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
	xi <- if (reduce[i]) vx[[i]] else x[[i]]
	for (j in 1:i) {
		if (w[i,j] == 0) next
		if (i == j) {
			a[[i]][,i] <- ux[[i]][,1]
			next
		}
		xj <- if (reduce[j]) vx[[j]] else x[[j]]
		svdij <- tryCatch(
			suppressWarnings(svds(tcrossprod(xi, xj), k = 1)), 
			error = function(e) svd(tcrossprod(xi, xj), 1, 1))
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
	}
}

## Reshape original data
for (i in 1:m) 
	dim(x[[i]]) <- dimx[[i]]

## Scale canonical tensors if needed
if (objective == "cor") {
	v <- scale.v(v, type = "var", x = x, check.args = FALSE)
}



#####################################
# Find best combination of canonical 
# vectors across datasets 
#####################################


score <- canon.scores(x, v)
best <- switch(optim, 
	exact = optim.combn.exact(score, w),
	greedy = optim.combn.greedy(score, w))
v <- v[cbind(1:m, best$idx)]
flip <- which(best$sign == -1)
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




v

}



