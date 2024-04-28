mcca.init.cca <- function(x, k = NULL, w = NULL, 
	objective = c("cov", "cor"), scale = c("block", "global"), 
	optim = NULL, maxit = 100, tol = .0001)
{
	
################
# Preprocessing
################


## Check arguments x and w 
test <- check.arguments(x, w = w)
eps <- 1e-14

## Scaling constraints
objective <- match.arg(objective)   
scale <- match.arg(scale) # block or global constraints
if (objective == "cor" && scale == "global")
	stop(paste("Argument values 'objective = cor'", 
		"and 'scale = global' are incompatible."))

## Search method in combinatorial optimization
optim <- if (is.null(optim)) {
	ifelse(m <= 5, "exact", "greedy")
} else {
	match.arg(optim, c("exact", "greedy")) 
}

## Data dimensions
m <- length(x)
dimfun <- function(x) if (is.vector(x)) c(1,length(x)) else dim(x)
dimx <- lapply(x, dimfun) 
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
n <- tail(dimx[[1]], 1)

## Objective weights
if (is.null(w)) {
	w <- 1 - diag(m)
} else if (length(w) == 1) {
	w <- matrix(1, m, m)
}
if (objective == "cor") diag(w) <- 0
w <- (w + t(w)) / (2 * sum(w)) 

## Check for constant datasets
cnst.set <- sapply(x, function(xx) all(abs(xx - xx[1]) < eps))
if (any(cnst.set)) ones <- function(len) rep(1/sqrt(len), len)
if ((objective == "cov" && all(cnst.set)) || 
	(objective == "cor" && sum(cnst.set) >= (m-1))) {
	v <- vector("list", m)
	for (i in 1:m) 
		v[[i]] <- lapply(p[[i]], ones)
	return(v)
} else if (any(cnst.set)) {
	m0 <- m
	dimx0 <- dimx
	p0 <- p
	x <- x[!cnst.set]
	m <- sum(!cnst.set)
	p <- p[!cnst.set]
	w <- w[!cnst.set, !cnst.set]
	w <- w / sum(w)	
}
pp <- sapply(p, prod)
d <- sapply(p, length)
k <- if (is.null(k)) 
	pmin(pp, n) else pmin(rep_len(k, m), pp, n)

## Data means for centering
xbar <- vector("list", m)
for(i in 1:m) 
    xbar[[i]] <- as.vector(rowMeans(x[[i]], dims = d[i]))



###########################
# SVD of unfolded datasets
###########################


## Calculate compact/truncated SVD of each unfolded dataset
ux <- vx <- vector("list", m)
for (i in 1:m) {
	xmat <- 	matrix(x[[i]] - xbar[[i]], pp[i], n)
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
	for (j in 1:i) {
		if (w[i,j] == 0) next
		if (j == i) {
			a[[i]][,i] <- ux[[i]][,1] 
			next
		}
		vxivxj <- crossprod(vx[[i]], vx[[j]])
		svdij <- tryCatch(supressWarnings(svds(vxivxj, k = 1)), 
			error = function(e) svd(vxivxj, 1, 1))
		a[[i]][,j] <- ux[[i]] %*% svdij$u
		a[[j]][,i] <- ux[[j]] %*% svdij$v
	}
}

rm(ux, vx, svdij)



#####################################
# Approximate long canonical vectors 
# by rank-1 tensors 
#####################################


v <- vector("list", m^2)
dim(v) <- c(m, m)
for (i in 1:m) {
	for (j in 1:m) {
		aij <- a[[i]][,j]
		if (d[i] > 1) dim(aij) <- p[[i]]
		v[[i,j]] <- switch(objective,
			cov = tnsr.rk1(aij, scale = TRUE, maxit = maxit, tol = tol),
			cor = tnsr.rk1.score(aij, cnstr = x[[i]] - xbar[[i]], 
				maxit = maxit, tol = tol))
	}
}


# After this stage, all canonical tensor weights 
# are scaled to have unit norm or unit variance 
# as required



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


if (objective == "cov" && scale == "global") {
	score <- canon.scores(x, v, FALSE)
	s <- eigen(w * cov(score) * w, TRUE)$vectors[,1]
	for (i in 1:m) v[[i]][[1]] <- s[i] * v[[i]][[1]]
}



#################
# Postprocessing
#################


## Drop singleton dimensions
for (i in 1:m)
	v[[i]] <- lapply(v[[i]], drop)

## Put back canonical tensors for constant datasets
if (any(cnst.set)) {
	vfull <- vector("list", m0)
	vfull[!cnst.set] <- v
	for (i in which(cnst.set))
		vfull[[i]] <- lapply(p0[[i]], ones)
	v <- vfull
}

v

}



