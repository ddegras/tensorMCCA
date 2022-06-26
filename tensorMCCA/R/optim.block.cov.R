###########################################
# FUNCTIONS TO MAXIMIZE SUM OF COVARIANCES 
# IN MCCA UNDER BLOCK CONSTRAINTS
###########################################



####################################
# Wrapper function for optimization
####################################



## Inputs
# v:	initial solution for canonical tensor
# a:	entire data tensor for quadratic term
# b:	reduced data tensor (weighted average over individuals) for linear term

optim.block.cov <- function(v, a, b, maxit = 1000, tol = 1e-6)
{
if (length(v) == 1) { return(optim1D.cov(a, b)) }
if (length(v) == 2) { return(optim2D.cov(v, a, b, maxit, tol)) }
if (length(v) == 3) { return(optim3D.cov(v, a, b, maxit, tol)) }
return(optim.gen.cov(v, a, b, maxit, tol))
}




#####################################
# Maximize { v'AA'v + 2 b'v } 
# subject to v'v = 1 
#####################################


optim1D.cov <- function(a, b)
{
eps <- 1e-14
n <- NCOL(a)
p <- length(b)
dim(b) <- NULL

## Trivial case: A = 0
if (all(a == 0)) {
	nrm <- sqrt(sum(b^2))
	if (nrm > eps) {
		return(list(b/nrm))
	} 
	return(list(numeric(p)))
}

## SVD of A
svdA <- svd(a, nu = p, nv = 0)
P <- svdA$u
delta <- (svdA$d)^2 
if (n < p) delta[(n+1):p] <- 0

## Trivial case: no linear term in objective
if (all(abs(b) <= eps)) {
	return(list(P[,1]))
}
## Trivial case: AA' proportional to identity matrix
if (abs(delta[p]-delta[1]) <= eps) {
	return(list(b/sqrt(sum(b^2))))
}
## Change of coordinates
b <- as.vector(crossprod(P, b))
bzero <- (abs(b) <= eps)
b[bzero] <- 0

objective.best <- -Inf

## Find stationary points of the Lagrange function in the case
## where the Lagrange multiplier is an eigenvalue of AA'
if (any(bzero)) {	
	diff.delta <- outer(-delta, delta[bzero], "+")
	equal.delta <- abs(diff.delta) <= eps
	test <- (colSums((!bzero) | equal.delta) == p)
	if (any(test)) {
		diff.delta <- diff.delta[, test, drop = FALSE]
		equal.delta <- equal.delta[, test, drop = FALSE]
		z <- b / diff.delta
		z[equal.delta] <- 0
		nrmz2 <- colSums(z^2)
		test <- (nrmz2 <= 1)
		if (any(test)) {
			if (!all(test)) {
				z <- z[,test, drop = FALSE]
				equal.delta <- equal.delta[, test, drop = FALSE]
				nrmz2 <- nrmz2[test]
			}
			nreplace <- colSums(equal.delta)
			vals <- sqrt((1 - nrmz2) / nreplace)
			z[equal.delta] <- rep.int(vals, nreplace)
			objective <- colSums(delta * z^2) + 2 * colSums(b * z)
			idx <- which.max(objective)
			v <- P %*% z[,idx]
			objective.best <- objective[idx]
		}	
	}			
}

## Find stationary points when the Lagrange multipliers
## are strictly between eigenvalues of AA'

## Formulate stationarity equation with unique eigenvalues 
i <- 1
count <- 1
bplus <- dplus <- numeric(0)
while(i <= p) {
	idx <- which(abs(delta - delta[i]) <= eps)
	dplus[count] <- delta[i]
	bplus[count] <- sqrt(sum(b[idx]^2))
	i <- max(idx) + 1
	count <- count + 1	
}
bplus <- rev(bplus)
dplus <- rev(dplus)
pplus <- length(bplus)

## Solve secular equation g(lambda) = 0 
## where lambda = Lagrange multiplier
g <- function(lambda) 
{
if (length(lambda) == 1)
	return(sum((bplus/(lambda-dplus))^2) - 1)
colSums((bplus / outer(-dplus, lambda, "+"))^2) - 1
}
lambda <- rep(NA, pplus + 1)
h <- 1.01 * sqrt(sum(bplus^2))
lb <- c(dplus[1] - h, dplus + 1e-3)
ub <- c(dplus - 1e-3, dplus[pplus] + h)
grid.size <- 21
for (i in 1:(pplus + 1)) {
	## Grid search 
	lam <- seq(lb[i], ub[i], length.out = grid.size)
	glam <- g(lam)
	sign.change <- (glam[-grid.size] * glam[-1] <= 0)
	## Refined search
	if (any(sign.change)) {
		idx <- which(sign.change)[1]
		lambda[i] <- uniroot(g, interval = c(lam[idx], lam[idx+1]),
			tol = 1e-6)$root
	} 
}
lambda <- lambda[!is.na(lambda)]

## Find corresponding solutions z and objective values
z <- b / outer(-delta, lambda, "+")
objective <- colSums(delta * z^2) + 2 * colSums(b * z)
idx <- which.max(objective)
if (objective[idx] > objective.best) 
	v <- P %*% z[,idx]
dim(v) <- NULL
v <- v / sqrt(sum(v^2))

list(v)
}


	

######################################################
# Maximize (1/n) sum_{t=1:n} < v,a(t) >^2 + 2 < v,b > 
# subject to < v,v > = 1  with a(t), b, v 2nd order 
# tensors (v of rank 1)
######################################################

## Inputs: 
# a:	3D array
# b:	matrix
# v:	list of 2 vectors

optim2D.cov <- function(v, a, b, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
dima <- dim(a)
p <- dima[1:2]
n <- dima[3]

## Trivial case
if (all(a == 0)) {
	svdb <- if (max(p) > 2) svds(b, 1) else svd(b, 1, 1)
	return(svdb[c("u", "v")])
} 

## MAIN LOOP
objective <- numeric(maxit)
for (it in 1:maxit) {

	## Update canonical vector in dimension 1 
	aa <- matrix(0, p[1], n)
	for (t in 1:n) aa[,t] <- a[,,t] %*% v[[2]]
	bb <- b %*% v[[2]]
	v[1] <- optim1D.cov(aa, bb)
		
	## Update canonical vector in dimension 2
	aa <- matrix(0, p[2], n)
	for (t in 1:n) aa[,t] <- crossprod(a[,,t], v[[1]])
	bb <- crossprod(b, v[[1]])
	v[2] <- optim1D.cov(aa, bb)
	
	## Calculate objective
	objective[it] <- mean(crossprod(v[[2]], aa)^2) + 
		2 * sum(bb * v[[2]])
	
	## Check convergence 
	if (it > 1 && abs(objective[it]-objective[it-1]) <= 
	    	tol * max(1,objective[it-1])) break
}

return(v)
}


	

######################################################
# Maximize (1/n) sum_{t=1:n} < v,a(t) >^2 + 2 < v,b > 
# subject to < v,v > = 1  with a(t), b, v 3rd order 
# tensors (v of rank 1)
######################################################


## Inputs: 
# v:	list of 3 vectors
# a:	4D array 
# b:	3D array

optim3D.cov <- function(v, a, b, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
dima <- dim(a)
p <- dima[1:3]
n <- dima[4]

## Trivial case
if (all(a == 0))
	return(tnsr3d.rk1(b, maxit, tol))

## MAIN LOOP
objective <- numeric(maxit)

for (it in 1:maxit) {

	## Update canonical vector in dimension 1 
	aa <- matrix(0, p[1], n)
	bb <- numeric(p[1])
	for (i in 1:p[1]) { 
		bb[i] <- crossprod(v[[2]], b[i,,] %*% v[[3]])
		for (t in 1:n) 
			aa[i,t] <- crossprod(v[[2]], a[i,,,t] %*% v[[3]])
	}
	v[[1]] <- unlist(optim1D.cov(aa, bb))
	
	## Update canonical vector in dimension 2
	aa <- matrix(0, p[2], n)
	bb <- numeric(p[2])
	for (j in 1:p[2]) { 
		bb[j] <- crossprod(v[[1]], b[,j,] %*% v[[3]])
		for (t in 1:n) 
			aa[j,t] <- crossprod(v[[1]], a[,j,,t] %*% v[[3]])
	}
	v[[2]] <- unlist(optim1D.cov(aa, bb))

	## Update canonical vector in dimension 3
	aa <- matrix(0, p[3], n)
	bb <- numeric(p[3])
	for (k in 1:p[3]) { 
		bb[k] <- crossprod(v[[1]], b[,,k] %*% v[[2]])
		for (t in 1:n) 
			aa[k,t] <- crossprod(v[[1]], a[,,k,t] %*% v[[2]])
	}
	v[[3]] <- unlist(optim1D.cov(aa, bb))

	## Calculate objective
	objective[it] <- mean(crossprod(v[[3]], aa)^2) + 
		2 * sum(bb * v[[3]])
	
	## Check convergence 
	if (it > 1 && abs(objective[it]-objective[it-1]) <= 
	    	tol * max(1,objective[it-1])) break
}

return(v)
}


######################################################
# Maximize (1/n) sum_{t=1:n} < v,a(t) >^2 + 2 < v,b > 
# subject to < v,v > = 1  with a(t), b, v general order 
# tensors (v of rank 1)
######################################################


## Inputs: 
# v:	list of d vectors
# a:	array with (d+1) dimensions 
# b:	array with d dimensions

optim.gen.cov <- function(v, a, b, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
dima <- dim(a)
d <- length(dima) - 1L
p <- dima[1:d]
n <- dima[d + 1L]

## Trivial case
if (all(a == 0))
	return(tnsr.rk1(b, maxit, tol))

## MAIN LOOP
objective <- numeric(maxit)

for (it in 1:maxit) {

	## Update canonical vector in dimension 1 
	aa <- matrix(0, p[1], n)
	bb <- numeric(p[1])
	for (i in 1:p[1]) { 
		bb[i] <- crossprod(v[[2]], b[i,,] %*% v[[3]])
		for (t in 1:n) 
			aa[i,t] <- crossprod(v[[2]], a[i,,,t] %*% v[[3]])
	}
	v[[1]] <- unlist(optim1D.cov(aa, bb))
	
	## Update canonical vector in dimension 2
	aa <- matrix(0, p[2], n)
	bb <- numeric(p[2])
	for (j in 1:p[2]) { 
		bb[j] <- crossprod(v[[1]], b[,j,] %*% v[[3]])
		for (t in 1:n) 
			aa[j,t] <- crossprod(v[[1]], a[,j,,t] %*% v[[3]])
	}
	v[[2]] <- unlist(optim1D.cov(aa, bb))

	## Update canonical vector in dimension 3
	aa <- matrix(0, p[3], n)
	bb <- numeric(p[3])
	for (k in 1:p[3]) { 
		bb[k] <- crossprod(v[[1]], b[,,k] %*% v[[2]])
		for (t in 1:n) 
			aa[k,t] <- crossprod(v[[1]], a[,,k,t] %*% v[[2]])
	}
	v[[3]] <- unlist(optim1D.cov(aa, bb))

	## Calculate objective
	objective[it] <- mean(crossprod(v[[3]], aa)^2) + 
		2 * sum(bb * v[[3]])
	
	## Check convergence 
	if (it > 1 && abs(objective[it]-objective[it-1]) <= 
	    	tol * max(1,objective[it-1])) break
}

return(v)
}



