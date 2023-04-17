
###################################################
# FUNCTIONS FOR MAXIMIZING THE SUM OF CORRELATIONS 
# IN MCCA ONE DATA BLOCK AT A TIME
###################################################



####################################
# Wrapper function for optimization
####################################


optim.block.cor <- function(v, obj, scale, ortho, maxit, tol)
{
if (length(v) == 1L) {
	optim1D.cor(obj, scale, ortho)
} else if (length(v) == 2L) {
	optim2D.cor(v, obj, scale, ortho, maxit, tol)
} else if (length(v) == 3L) {
	optim3D.cor(v, obj, scale, ortho, maxit, tol)
} else {
	optim.gen.cor(v, obj, scale, ortho, maxit, tol)
}
}





##########################################
# Maximize linear form a'v subject to
# (1/n) v'BB'v = 1 where B has n columns 
# and to C'v = 0 
##########################################


## Inputs: 
## a:	objective (vector)
## b:	scaling constraints (matrix)
## cc:	orthogonality constraints (list of vectors or matrix)

optim1D.cor <- function(a, b, cc)
{
p <- length(a)
if (!is.null(cc)) {
	if (is.list(cc)) 
		cc <- matrix(unlist(cc), p, length(cc))
	qrc <- qr(cc)
	if (qrc$rank == p) return(list(numeric(p)))
	if (qrc$rank > 0) {
		Q <- qr.Q(qrc, complete = TRUE)[, -(1:qrc$rank), drop = FALSE]
		a <- crossprod(Q, a)
		b <- crossprod(Q, b)
	} else cc <- NULL
}
if (length(a) == 1) {
	v <- as.numeric(sign(a)) / sqrt(mean(b^2))
	if (is.nan(v)) v <- 0
	if (!is.null(cc)) v <- v * Q
	return(list(as.vector(v)))
}
eps <- 1e-14
n <- ncol(b)
sigma <- tcrossprod(b) / n
v <- tryCatch(solve(sigma, a), error = function(e) numeric(0))
if (length(v) == 0) v <- ginv(sigma) %*% a
s <- sum(a * v)
v <- if (s > eps) v / sqrt(s) else numeric(length(v))
if (!is.null(cc)) v <- Q %*% v
dim(v) <- NULL
list(v)
}



##########################################
# Maximize bilinear form v1' A v2 under
# constraint (1/n) sum_t (v1' Bt v2)^2 = 1
# and v1' C(l) v2 = 0 for l = 1,2,...
##########################################


## Inputs: 
## v:	initial solution (list of 2 vectors)
## a:	objective (matrix)
## b:	scaling constraints (3D array)
## cc:	orthogonality constraints (list of matrices)


optim2D.cor <- function(v, a, b, cc, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
b <- aperm(b, c(1, 3, 2))
dimb <- dim(b)
p <- dimb[c(1,3)]
n <- dimb[2]
if (!is.matrix(a)) dim(a) <- p
if (!is.null(cc)) {
	if (is.list(cc)) 
		cc <- array(unlist(cc), c(p, length(cc)))
	if (all(abs(cc) <= 1e-14)) {
		cc <- NULL
	} else if (any(p == 1)) {
		 return(lapply(p, numeric))
	} else { cc <- aperm(cc, c(1, 3, 2)) }
	northo <- dim(cc)[2]
}


## MAIN LOOP
objective <- -Inf
ccmat <- NULL
for (it in 1:maxit) {
	objective.old <- objective
	
	## Update canonical vector in dimension 1 
	if (p[2] == 1) {
		aa <- as.vector(a) * v[[2]]		
		bb <- matrix(b * v[[2]], p[1], n)
	} else {	
		aa <- a %*% v[[2]]
		dim(b) <- c(p[1] * n, p[2])
		bb <- b %*% v[[2]]
		dim(bb) <- c(p[1], n)
	}
	if (!is.null(cc)) {
		dim(cc) <- c(p[1] * northo, p[2])
		ccmat <- cc %*% v[[2]]
		dim(ccmat) <- c(p[1], northo)
	}
	v[1] <- optim1D.cor(aa, bb, ccmat)
		
	## Update canonical vector in dimension 2
	if (p[1] == 1) {
		aa <- as.vector(a) * v[[1]]
		bb <- matrix(b * v[[1]], p[2], n, byrow = TRUE)
	} else {
		aa <- crossprod(a, v[[1]])
		dim(b) <- c(p[1], n * p[2])
		bb <- crossprod(v[[1]], b)
		dim(bb) <- c(n, p[2])
		bb <- t(bb)
	}
	if (!is.null(cc)) {
		dim(cc) <- c(p[1], northo * p[2])
		ccmat <- crossprod(v[[1]], cc)
		dim(ccmat) <- c(northo, p[2])
		ccmat <- t(ccmat)
	}
	v[2] <- optim1D.cor(aa, bb, ccmat)
	
	## Calculate objective
	objective <- sum(aa * v[[2]])
	
	## Check convergence 
	if (it > 1 && abs(objective - objective.old) <= 
	    	tol * max(1, objective.old)) break
}

v
}


	

###################################################
# Maximize inner product of rank-1 3D tensor 
# with fixed 3D tensor under quadratic scaling 
# constraints and orthogonality constraints
# max_v < v,a > with v rank 1 subject to 
# (1/n) sum_t (< v,b(t) >)^2 = 1
# and < C(l), v > = 0 for l = 1, 2, ...
###################################################

## Inputs: 
## v:	initial solution (list of 3 vectors)
## a:	objective (3D array) 
## b:	scaling constraints (4D array) 
## cc:	orthogonality constraints (list of 3D arrays)

optim3D.cor <- function(v, a, b, cc, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
dimb <- dim(b)
p <- dimb[1:3]
n <- dimb[4]

## Reshape orthogonality constraints
if (!is.null(cc)) {
	if (is.list(cc)) 
		cc <- array(unlist(cc), c(p, length(cc)))
	northo <- dim(cc)[4]
	if (all(abs(cc) <= 1e-14)) cc <- NULL
}

## MAIN LOOP
objective <- -Inf
ccmat <- NULL
for (it in 1:maxit) {
	objective.old <- objective
	
	## Update canonical vector in dimension k = 1,2,3
	for (k in 1:3) {
		i1 <- if (k == 1L) 2L else 1L
		i2 <- if (k == 3L) 2L else 3L
		aa <- aperm(a, c(i1, k, i2))
		bb <- aperm(b, c(i1, k, 4, i2))
		if (p[i1] == 1L && p[i2] == 1L) {
			aa <- (v[[i1]] * v[[i2]]) * drop(aa)
			bb <- (v[[i1]] * v[[i2]]) * drop(b)
		} else if (p[i1] == 1L) {
			aa <- if (p[k] == 1L) {
				v[[i1]] * sum(aa * v[[i2]])
			} else {
				drop(aa) %*% (v[[i1]] * v[[i2]])
			}
			dim(bb) <- c(p[k] * n, p[i2])
			bb <- bb %*% (v[[i1]] * v[[i2]])
		} else if (p[i2] == 1L) {
			aa <- if (p[k] == 1) {
				v[[i2]] * sum(aa * v[[i1]])
			} else {
				crossprod(v[[i1]] * v[[i2]], drop(aa)) 
			}
			dim(bb) <- c(p[i1], p[k] * n)
			bb <- crossprod(v[[i1]] * v[[i2]], bb)
		} else {		
			dim(aa) <- c(p[i1], p[k] * p[i2])
			aa <- crossprod(v[[i1]], aa) 
			dim(aa) <- c(p[k], p[i2])
			aa <- aa %*% v[[i2]]
			dim(bb) <- c(p[i1], p[k] * n * p[i2])
			bb <- crossprod(v[[i1]], bb)
			dim(bb) <- c(p[k] * n, p[i2])
			bb <- bb %*% v[[i2]]
		}
		dim(bb) <- c(p[k], n)	
		if (!is.null(cc)) {
			ccmat <- aperm(cc, c(i1, k, 4, i2))
			dim(ccmat) <- c(p[i1], p[k] * northo * p[i2])
			ccmat <- crossprod(v[[i1]], ccmat)
			dim(ccmat) <- c(p[k] * northo, p[i2])
			ccmat <- ccmat %*% v[[i2]]
			dim(ccmat) <- c(p[k], northo)
		}
		v[k] <- optim1D.cor(aa, bb, ccmat)
	}
		
	## Calculate objective
	objective <- sum(aa * v[[3]]) 
	
	## Check convergence 
	if (abs(objective - objective.old) <= 
	    	tol * max(1, objective.old)) break
}

v
}



########################################################
# Maximize inner product of rank-1 tensor 
# with fixed tensor under quadratic scaling constraints
# and orthogonality constraints
# max_v < v,a > subject to (1/n) sum_t (< v,b(t) >)^2 = 1
# and < c(l),v > = 0 for l = 1, 2, ...
########################################################


## Inputs: 
# v:	list of d vectors
# a:	array with d dimensions 
# b:	array with (d+1) dimensions
# cc:	list of d-dimensional arrays
# The lengths of the vectors in v, the dimensions of a, 
# and the first d dimensions of b must all match

optim.gen.cor <- function(v, a, b, cc, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
dima <- dim(a)
d <- length(v)

## MAIN LOOP
objective <- -Inf
ccmat <- NULL
for (it in 1:maxit) {

	objective.old <- objective

	## Update canonical vector in dimension k = 1,...,d
	for (k in 1:d) { 
		aa <- tnsr.vec.prod(a, v[-k], (1:d)[-k])
		bb <- tnsr.vec.prod(b, v[-k], (1:d)[-k])
		if (!is.null(cc)) 
			ccmat <- sapply(cc, tnsr.vec.prod, 
				v = v[-k], modes = (1:d)[-k])
		v[k] <- optim1D.cor(aa, bb, ccmat)
	}

	## Calculate objective
	objective[it] <- if (dima[d] == 1) {
		v[[d]]^2 * mean(aa^2) + 2 * v[[d]] * bb
	} else { mean(crossprod(v[[d]], aa)^2) + 
		2 * sum(bb * v[[d]]) }
	
	## Check convergence 
	if (it > 1 && abs(objective - objective.old) <= 
	    	tol * max(1, objective.old)) break
}

return(v)
}


