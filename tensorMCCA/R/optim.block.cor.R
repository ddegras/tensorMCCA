
###################################################
# FUNCTIONS FOR MAXIMIZING THE SUM OF CORRELATIONS 
# IN MCCA ONE DATA BLOCK AT A TIME
###################################################



####################################
# Wrapper function for optimization
####################################


optim.block.cor <- function(v, a, b, maxit, tol)
{
if (length(v) == 1L) {
	optim1D.cor(a, b)
} else if (length(v) == 2L) {
	optim2D.cor(v, a, b, maxit, tol)
} else if (length(v) == 3L) {
	optim3D.cor(v, a, b, maxit, tol)
} else {
	optim.gen.cor(v, a, b, maxit, tol)
}
}





##########################################
# Maximize linear form v'a subject to
# (1/n) v'BB'v = 1 where B has n columns  
##########################################


## Inputs: 
## a:	objective (vector)
## b:	constraints (matrix)

optim1D.cor <- function(a, b)
{
if (length(a) == 1) 
	return(sign(a) / sqrt(mean(b^2)))
eps <- 1e-14
n <- ncol(b)
sigma <- tcrossprod(b) / n
v <- tryCatch(solve(sigma, a), error = function(e) numeric(0))
if (length(v) == 0) v <- ginv(sigma) %*% a
dim(v) <- NULL
s <- sum(a * v)
v <- if (s > eps) v / sqrt(s) else numeric(length(v))
list(v)
}



##########################################
# Maximize bilinear form v1' A v2 under
# constraint (1/n) sum_t (v1' Bt v2)^2 = 1
##########################################


## Inputs: 
## v:	initial solution (list of 2 vectors)
## a:	objective (2D matrix)
## b:	constraints (3D array)

optim2D.cor <- function(v, a, b, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
b <- aperm(b, c(1, 3, 2))
dimb <- dim(b)
p <- dimb[c(1,3)]
n <- dimb[2]
if (!is.matrix(a)) dim(a) <- p

## MAIN LOOP
objective <- -Inf
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
	v[1] <- optim1D.cor(aa, bb)
		
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
	v[2] <- optim1D.cor(aa, bb)
	
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
# with fixed 3D tensor under quadratic constraints
# max < v,a > subject to (1/n) sum_t (< v,b(t) >)^2 = 1
###################################################


## Inputs: 
## v:	initial solution (list of 3 vectors)
## a:	objective (3D array) 
## b:	constraints (4D array) 

optim3D.cor <- function(v, a, b, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
dimb <- dim(b)
p <- dimb[1:3]
n <- dimb[4]

## MAIN LOOP
objective <- -Inf

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
		v[k] <- optim1D.cor(aa, bb)
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
# with fixed tensor under quadratic constraints
# max < v,a > subject to (1/n) sum_t (< v,b(t) >)^2 = 1
########################################################


## Inputs: 
# v:	 list of d vectors
# a:	 array with d dimensions 
# b:	 array with (d+1) dimensions

# The lengths of the vectors in v, the dimensions of a, 
# and the first d dimensions of b must all match

optim.gen.cor <- function(v, a, b, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
dima <- dim(a)
d <- length(v)


## MAIN LOOP
objective <- -Inf

for (it in 1:maxit) {

	objective.old <- objective

	## Update canonical vector in dimension k = 1,...,d
	for (k in 1:d) { 
		aa <- tnsr.vec.prod(a, v[-k], (1:d)[-k])
		bb <- tnsr.vec.prod(b, v[-k], (1:d)[-k])
		v[k] <- optim1D.cor(aa, bb)
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


