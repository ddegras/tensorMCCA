mcca.init.cca <- function(x, w = NULL, objective = c("cov", "cor"), 
	scope = c("block", "global"), k = NULL, optim = NULL, 
	maxit = 100, tol = .0001)
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
	m <- sum(wnzero)
	p <- p[wnzero]
	w <- w[wnzero, wnzero]
}
pp <- sapply(p, prod)
d <- sapply(p, length)
k <- if (is.null(k)) 
	pmin(pp, n) else pmin(rep_len(k, m), pp, n)





###########################
# SVD of unfolded datasets
###########################


## Calculate compact/truncated SVD of each unfolded dataset
ux <- vx <- vector("list", m)
for (i in 1:m) {
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


if (objective == "cov" && scope == "global") {
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

## Put back any canonical weights for constant datasets
if (any(wzero)) {
	vfull[wnzero] <- v
	v <- vfull
}

## Scale solution
if (objective == "cor") {
	v[wnzero] <- scale.v(v, x = x, check.args = FALSE,
		type = switch(objective, cov = "norm", cor = "var"))
} else {
	v <- scale.v(v, "norm", scope, check.args = FALSE)
}

v

}



