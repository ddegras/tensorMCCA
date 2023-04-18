1###########################################
# FUNCTIONS TO MAXIMIZE SUM OF COVARIANCES 
# IN MCCA UNDER BLOCK CONSTRAINTS
###########################################



####################################
# Wrapper function for optimization
####################################



## Inputs
# v: initial solution for canonical tensor
# a:	 entire data tensor for quadratic term
# b:	 reduced data tensor (weighted average over individuals) for linear term

optim.block.cov <- function(v, a, b, ortho, maxit = 1000, tol = 1e-6)
{
if (length(v) == 1L) { return(optim1D.cov(a, b, ortho)) }
if (length(v) == 2L) { return(optim2D.cov(v, a, b, ortho, maxit, tol)) }
if (length(v) == 3L) { return(optim3D.cov(v, a, b, ortho, maxit, tol)) }
return(optim.gen.cov(v, a, b, ortho, maxit, tol))
}




#####################################
# Maximize { v'AA'v + 2 v'b } 
# subject to v'v = 1 and C'v = 0
#####################################

## Inputs
# a: matrix of dimensions p-by-n
# b: vector of length p
# cc: matrix with p rows

optim1D.cov <- function(a, b, cc)
{
eps <- 1e-14
n <- NCOL(a)
p <- length(b)
dim(b) <- NULL

## Trivial case: v scalar
if (p == 1L) {
	test <- (is.null(cc) || all(abs(cc) <= eps))
	v <- if (test && b >= 0) { 1
	} else if (test && b < 0) { -1
	} else { 0 }
	return(list(v))
}

## Handle orthogonality constraints by changing variables
if (!is.null(cc)) {
	if (is.list(cc)) 
		cc <- matrix(unlist(cc), p, length(cc))	
	qrc <- qr(cc)
	if (qrc$rank == p) return(list(numeric(p)))
	qq <- qr.Q(qrc, complete = TRUE)[, -(1:qrc$rank), drop = FALSE]
	if (is.matrix(a)) a <- crossprod(qq, a)
	b <- as.vector(crossprod(qq, b))
}

## Trivial case: A = 0
if (all(abs(a) <= eps)) {
	nrm <- sqrt(sum(b^2))
	v <- if (nrm > eps && is.null(cc)) { 
		b / nrm 
	} else if (nrm > eps && (!is.null(cc))) {
		as.vector(qq %*% b) / nrm	
	} else if (nrm <= eps && is.null(cc)) {
		rep(1/sqrt(p), p)
	} else if (nrm <= eps && (!is.null(cc))) {
		qq[,1]
	}
	return(list(v))
}

## SVD of A
svda <- svd(a, nv = 0)
pos <- (svda$d >= 1e-8 * svda$d[1])
P <- svda$u[, pos, drop = FALSE]
delta <- (svda$d[pos])^2 
pp <- sum(pos) 
# pp = number of optimization variables in transformed problem 

## Change variables
b <- as.vector(crossprod(P, b))

## Trivial case: no linear term in objective
bzero <- (abs(b) <= eps)
if (all(bzero)) {
	v <- if (is.null(cc)) {
		P[,1] } else { as.vector(qq %*% P[,1]) }
	return(list(v))
}
b[bzero] <- 0

objective.best <- -Inf

## Find stationary points of the Lagrange function in the case
## where the Lagrange multiplier is an eigenvalue of AA'
if (any(bzero)) {	
	diff.delta <- outer(-delta, delta[bzero], "+")
	equal.delta <- abs(diff.delta) <= eps
	test <- (colSums((!bzero) | equal.delta) == pp)
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
while(i <= pp) {
	idx <- which(abs(delta - delta[i]) <= eps)
	dplus[count] <- delta[i]
	bplus[count] <- sqrt(sum(b[idx]^2))
	i <- max(idx) + 1
	count <- count + 1	
}
bplus <- rev(bplus)
dplus <- rev(dplus)
nz <- (bplus > eps)
bplus <- bplus[nz]
dplus <- dplus[nz]
pplus <- length(bplus)

## Solve secular equation g(lambda) = 0 
## where lambda = Lagrange multiplier
g <- function(lambda) 
{
if (length(lambda) == 1)
	return(sum((bplus/(lambda-dplus))^2) - 1)
colSums((bplus / outer(-dplus, lambda, "+"))^2) - 1
}
lambda <- lambda.opt <- numeric()
h <- 1.01 * sqrt(sum(bplus^2))
lb <- c(dplus[1] - h, dplus + 1e-5)
ub <- c(dplus - 1e-5, dplus[pplus] + h)
grid.size <- 21
for (i in 1:(pplus + 1)) {
	## Grid search 
	lam <- seq(lb[i], ub[i], length.out = grid.size)
	glam <- g(lam)
	sgn.glam <- sign(glam)
	if (any(sgn.glam == 0)) {
		lambda <- c(lambda, lam[sgn.lam == 0])
		next }
	sign.change <- (sgn.glam[-grid.size] != sgn.glam[-1])
	## Refined search
	if (any(sign.change)) {
		idx <- which(sign.change)[1]
		root <- uniroot(g, interval = c(lam[idx], lam[idx+1]),
			tol = 1e-7)$root
		lambda <- c(lambda, root) 			
	} else {
		lambda.opt <- c(lambda.opt, lam[which.min(abs(glam))])
	}
}
if (length(lambda) == 0) { lambda <- lambda.opt }


## Find corresponding solutions z and objective values
z <- b / outer(-delta, lambda, "+")
nrmz <- sqrt(colSums(z^2))
z <- sweep(z, 2, nrmz, "/")
objective <- colSums(delta * z^2) + 2 * colSums(b * z)
idx <- which.max(objective)
if (objective[idx] > objective.best) { v <- P %*% z[,idx] }
v <- v / sqrt(sum(v^2))
if (!is.null(cc)) v <- qq %*% v
dim(v) <- NULL
list(v)
}


	

##########################################################
# Maximize (1/n) sum_{t=1:n} (v1' A(t) v2)^2 + 2 v1' B v2  
# subject to || v1 ||^2 = || v2 ||^2 = 1 
# and v1' C(l) v2 = 0 for l = 1,2,... 
# with A(t), B, C(l) matrices and v1, v2 vectors
##########################################################

## Inputs: 
# v:	list of 2 vectors (initial solution)
# a:	3D array (quadratic component of objective)
# b:	matrix (linear component)
# cc:	list of matrices or NULL (orthogonality constraints)

optim2D.cov <- function(v, a, b, cc, maxit = 1000, tol = 1e-6)
{
## Data dimensions
p <- sapply(v, length)
if (length(dim(a)) == 3L) {
	a <- aperm(a, c(1, 3, 2))
	dima <- dim(a)
	n <- dima[2]
}

## Reshape orthogonality constraints if any
if (!is.null(cc)) {
	northo <- length(cc)
	cc <- array(unlist(cc), c(p, northo))
}

## Trivial case: A(t) = 0 for all t and C(l) = 0 for all l
azero <- all(a == 0) 
if (azero && is.null(cc)) {
	if (all(b == 0)) 
		return(lapply(p, function(len) rep(1/sqrt(len), len)))
	svdb <- NULL
	if (min(p) > 2) 
		svdb <- tryCatch(svds(b, 1), error = function(e) NULL) 
	if (is.null(svdb) || any(is.nan(c(svdb$u, svdb$v))))
		svdb <- svd(b, 1, 1)
	return(list(as.numeric(svdb$u), as.numeric(svdb$v)))
}

## Trivial case: A(t) = 0 for all t and B = 0
if (azero && !is.null(cc) && all(b == 0) && northo < max(p)) {
	if (northo < p[1]) {
		qrc <- qr(cc[,1,])
		return(list(qr.Q(qrc, complete = TRUE)[, qrc$rank + 1],
			rep(c(1,0), c(1,p[2]-1))))
	} else {
		qrc <- qr(cc[1,,])
		return(list(rep(c(1,0), c(1,p[1]-1),
			qr.Q(qrc, complete = TRUE)[, qrc$rank + 1])))
	}
}

## MAIN LOOP
objective <- numeric(maxit)
aa <- 0
if (!is.null(cc)) cc <- aperm(cc, c(1,3,2)) # new dims p1 x northo x p2
ccmat <- NULL

for (it in 1:maxit) {

	## Update canonical vector in dimension 1 
	if (p[2] == 1L) {
		if (!azero) aa <- matrix(a * v[[2]], p[1], n)
		bb <- as.vector(b * v[[2]])
	} else {
		if (!azero) {
			dim(a) <- c(p[1] * n, p[2])
			aa <- a %*% v[[2]]
			dim(aa) <- c(p[1], n)
			dim(a) <- dima }
		bb <- b %*% v[[2]]
	}
	if (!is.null(cc)) {
		dim(cc) <- c(p[1] * northo, p[2])
		ccmat <- cc %*% v[[2]]
		dim(ccmat) <- c(p[1], northo)
	}
	v[1] <- optim1D.cov(aa, bb, ccmat)
		
	## Update canonical vector in dimension 2
	if (p[1] == 1L) {
		if (!azero) aa <- matrix(a * v[[1]], p[2], n)
		bb <- as.vector(b) * v[[1]]
	} else {
		if (!azero) {
			dim(a) <- c(p[1], n * p[2])
			aa <- crossprod(v[[1]], a)
			dim(aa) <- c(n, p[2])
			aa <- t(aa)
			dim(a) <- dima 
		}
		bb <- crossprod(b, v[[1]])
	}
	if (!is.null(cc)) {
		dim(cc) <- c(p[1], northo * p[2])
		ccmat <- crossprod(v[[1]], cc)
		dim(ccmat) <- c(northo, p[2])
		ccmat <- t(ccmat)
	}
	v[2] <- optim1D.cov(aa, bb, ccmat)
	
	## Calculate objective
	objective[it] <- if (p[2] == 1) {
		v[[2]]^2 * mean(aa^2) + 2 * v[[2]] * bb
	} else if (azero) {
		2 * sum(bb * v[[2]])
	} else { mean(crossprod(v[[2]], aa)^2) + 
		2 * sum(bb * v[[2]]) }
	
	## Check convergence 
	if (it > 1 && abs(objective[it] - objective[it-1]) <= 
	    	tol * max(1, objective[it-1])) break
}

return(v)
}


	

######################################################
# Maximize (1/n) sum_{t=1:n} < v,A(t) >^2 + 2 < v,B > 
# subject to < v,v > = 1 and < C(l),v > = 0 for l=1,2,...
# with A(t), B, v 3rd order tensor arrays, v rank-1
######################################################


## Inputs: 
# v:	list of 3 vectors (initial solution)
# a:	4D array (quadratic component of objective)
# b:	3D array (linear component of objective)
# cc: 	list of 3D arrays (orthogonality constraints)

optim3D.cov <- function(v, a, b, cc, maxit = 1000, tol = 1e-6)
{
	
## Trivial case
if (all(a == 0) && is.null(cc)) 
	return(tnsr3d.rk1(b, scale = TRUE, maxit, tol))

## Data dimensions
p <- dim(b)
anzero <- !all(a == 0)
if (!anzero) {	
	a <- 0
} else {
	n <- dim(a)[4]
}

## Reshape orthogonality constraints
if (!is.null(cc)) {
	if (is.list(cc)) 
		cc <- array(unlist(cc), c(p, length(cc)))
	northo <- dim(cc)[4]
	if (all(abs(cc) <= 1e-14)) cc <- NULL
}
 
## MAIN LOOP
objective <- numeric(maxit)
ccmat <- NULL
for (it in 1:maxit) {
	
	## Update canonical vector in dimension k = 1,2,3
	for (k in 1:3) {
		i1 <- if (k == 1L) 2L else 1L
		i2 <- if (k == 3L) 2L else 3L		 
		aa <- if (anzero) aperm(a, c(i1, k, 4, i2)) else 0
		bb <- aperm(b, c(i1, k, i2))
		
		if (p[i1] == 1 && p[i2] == 1) {
			if (anzero) 
				aa <- (v[[i1]] * v[[i2]]) * drop(aa)
			bb <- (v[[i1]] * v[[i2]]) * drop(bb)
		} else if (p[i1] == 1) {
			if (anzero) {
				dim(aa) <- c(p[k] * n, p[i2])
				aa <- aa %*% (v[[i1]] * v[[i2]])
			}
			bb <- if (p[k] == 1) {
				v[[i1]] * sum(bb * v[[i2]])
			} else {
				drop(bb) %*% (v[[i1]] * v[[i2]])
			}
		} else if (p[i2] == 1) {
			if (anzero) {
				dim(aa) <- c(p[i1], p[k] * n)
				aa <- crossprod(v[[i1]] * v[[i2]], aa)
			}
			bb <- if (p[k] == 1) {
				v[[i2]] * sum(bb * v[[i1]])
			} else {
				crossprod(v[[i1]] * v[[i2]], drop(bb)) 
			}
		} else {
			if (anzero) {	
				dim(aa) <- c(p[i1], p[k] * n * p[i2])
				aa <- crossprod(v[[i1]], aa)
				dim(aa) <- c(p[k] * n, p[i2])
				aa <- aa %*% v[[i2]]
			}
			dim(bb) <- c(p[i1], p[k] * p[i2])
			bb <- crossprod(v[[i1]], bb) 
			dim(bb) <- c(p[k], p[i2])
			bb <- bb %*% v[[i2]]
		}
		if (anzero) dim(aa) <- c(p[k], n)	
		if (!is.null(cc)) {
			ccmat <- aperm(cc, c(i1, k, 4, i2))
			dim(ccmat) <- c(p[i1], p[k] * northo * p[i2])
			ccmat <- crossprod(v[[i1]], ccmat)
			dim(ccmat) <- c(p[k] * northo, p[i2])
			ccmat <- ccmat %*% v[[i2]]
			dim(ccmat) <- c(p[k], northo)
		}
		v[k] <- optim1D.cov(aa, bb, ccmat)
	}
	
	## Calculate objective
	objective[it] <- if (p[[3]] == 1) {
		v[[3]]^2 * mean(aa^2) + 2  * v[[3]] * bb
	} else if (anzero) { 
		mean(crossprod(v[[3]], aa)^2) + 2 * sum(bb * v[[3]]) 
	} else { 
		2 * sum(bb * v[[3]]) 
	}
	
	## Check convergence 
	if (it > 1 && abs(objective[it] - objective[it-1]) <= 
	    	tol * max(1, objective[it-1])) break
}

return(v)
}


######################################################
# max_v (1/n) sum_{t=1:n} < v,a(t) >^2 + 2 < v,b > 
# subject to < v,v > = 1  and < c(l),v > = 0, l=1,2,...
# with a(t), b, c(l), and v tensors of same dimensions
# (v of rank 1)
######################################################


## Inputs: 
# v:	list of d vectors
# a:	(d+1)-dimensional array  
# b:	d-dimensional array 
# cc:	list of d-dimensional arrays

optim.gen.cov <- function(v, a, b, cc, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
dima <- dim(a)
d <- length(v)

## Trivial case
if (all(a == 0) && is.null(cc))
	return(tnsr.rk1(b, maxit, tol))

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
		v[k] <- optim1D.cov(aa, bb, ccmat)
	}

	## Calculate objective
	objective[it] <- if (dima[d] == 1) {
		v[[d]]^2 * mean(aa^2) + 2 * v[[d]] * bb
	} else { mean(crossprod(v[[d]], aa)^2) + 
		2 * sum(bb * v[[d]]) }
	
	## Check convergence 
	if (abs(objective - objective.old) <= 
	    	tol * max(1, objective.old)) break
}

v
}



