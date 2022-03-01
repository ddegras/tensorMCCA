
# FUNCTIONS FOR MAXIMIZING THE SUM OF CORRELATIONS 
# IN MCCA ONE DATA BLOCK AT A TIME




##########################################
# Maximize scalar product form v' a under
# constraint (1/n) sum_(t=1:n) (v' bt)^2 = 1
##########################################


## Inputs: 
## a:	objective (vector)
## b:	constraints (matrix)

optim1D.cor <- function(a, b)
{
n <- ncol(b)
sigma <- tcrossprod(b) / n
v <- tryCatch(solve(sigma, a), error = function(e) numeric(0))
if (length(v) == 0) v <- ginv(sigma) %*% a
s <- as.numeric(crossprod(a, v))
if (s > 0) v <- v / sqrt(s)
list(v)
}



##########################################
# Maximize bilinear form v1' A v2 under
# constraint sum_t (v1' Bt v2)^2 = 1
##########################################


## Inputs: 
## v:	initial solution (list of 2 vectors)
## a:	objective (2D matrix)
## b:	constraints (3D array)

optim2D.cor <- function(v, a, b, maxit = 1000, tol = 1e-6, debug = TRUE)
{
	
## Data dimensions
dimv <- sapply(v,length)
dimb <- dim(b)
p <- dimv[1]
q <- dimv[2]
n <- dimb[3]

## Check compatibility of dimensions
stopifnot(identical(dim(a),dimv))
stopifnot(identical(dim(b)[1:2],dimv))


## MAIN LOOP
objective <- numeric(maxit)
# browser()
for (it in 1:maxit) {

	## Update canonical vector in dimension 1 
	c <- a %*% v[[2]]
	d <- matrix(0, p, n)
	for (t in 1:n) d[,t] <- b[,,t] %*% v[[2]]
	# sigma <- tcrossprod(d) / n
	# v[[1]] <- tryCatch(solve(sigma, c), error = function(e) numeric(0))
	# if (length(v[[1]]) == 0) v[[1]] <- ginv(sigma) %*% c
	# cp <- as.numeric(crossprod(v[[1]], sigma %*% v[[1]]))
	# v[[1]] <- v[[1]] / sqrt(cp)
	v[[1]] <- unlist(optim1D.cor(c, d))
	
	if (debug) {
		check1 <- crossprod(c, v[[1]]) # debugging == sqrt(cp) ?
		cat("\n In optim2D.cor Iteration",it,"Objective",check1)
	}
	
	## Update canonical vector in dimension 2
	c <- crossprod(a, v[[1]])
	d <- matrix(0, q, n)
	for (t in 1:n) d[,t] <- crossprod(b[,,t], v[[1]])
	# sigma <- tcrossprod(d) / n
	# v[[2]] <- tryCatch(solve(sigma, c), error = function(e) numeric(0))
	# if (length(v[[2]]) == 0) v[[2]] <- ginv(sigma) %*% c
	# cp <- as.numeric(crossprod(v[[2]], c))
	# if (cp > 0) v[[2]] <- v[[2]] / sqrt(cp)
	v[[2]] <- unlist(optim1D.cor(c, d))
	
	## Calculate objective
	objective[it] <- crossprod(c, v[[2]]) # = sqrt(cp) ?
	if (debug) {
		cat(" ", objective[it])
	}
	
	## Check convergence 
	if (it > 1 && abs(objective[it]-objective[it-1]) <= 
	    	tol * max(1,objective[it-1])) break
}

return(v)
}


	

###################################################
# Maximize inner product of rank-1 3D tensor 
# with fixed 3D tensor under quadratic constraints
# max < v,a > subject to sum_t (< v,b(t) >)^2 = 1
###################################################


## Inputs: 
## v:	initial solution (list of 3 vectors)
## a:	objective (3D array) 
## b:	constraints (4D array) 

optim3D.cor <- function(v, a, b, maxit = 1000, tol = 1e-6, debug = TRUE)
{
	
## Data dimensions
dimv <- sapply(v, length)
dimb <- dim(b)
p <- dimv[1]
q <- dimv[2]
r <- dimv[3]
n <- dimb[4]

## Check compatibility of dimensions
stopifnot(identical(dim(a),dimv))
stopifnot(identical(dim(b)[1:3],dimv))


## MAIN LOOP
objective <- numeric(maxit)

for (it in 1:maxit) {

	## Update canonical vector in dimension 1 
	c <- numeric(p)
	d <- matrix(0,p,n)
	for (i in 1:p) { 
		c[i] <- crossprod(v[[2]], a[i,,] %*% v[[3]])
		for (t in 1:n) 
			d[i,t] <- crossprod(v[[2]], b[i,,,t] %*% v[[3]])
	}
		# sigma <- tcrossprod(d) / n
		# v[[1]] <- tryCatch(solve(sigma, c), error = function(e) numeric(0))
		# if (length(v[[1]]) == 0) v[[1]] <- ginv(sigma) %*% c
		# cp <- as.numeric(crossprod(v[[1]], sigma %*% v[[1]]))
		# if (cp > 0) v[[1]] <- v[[1]] / sqrt(cp)
	v[[1]] <- unlist(optim1D.cor(c, d))
	if (debug) {
		check1 <- crossprod(c, v[[1]]) # debugging == sqrt(cp) ?
		cat("\nIn optim3D.cor Iteration",it,"Objective",check1)
	}
	
	## Update canonical vector in dimension 2
	c <- numeric(q)
	d <- matrix(0,q,n)
	for (j in 1:q) { 
		c[j] <- crossprod(v[[1]], a[,j,] %*% v[[3]])
		for (t in 1:n) 
			d[j,t] <- crossprod(v[[1]], b[,j,,t] %*% v[[3]])
	}
	v[[2]] <- unlist(optim1D.cor(c, d))
		# sigma <- tcrossprod(d) / n
		# v[[2]] <- tryCatch(solve(sigma, c), error = function(e) numeric(0))
		# if (length(v[[2]]) == 0) v[[2]] <- ginv(sigma) %*% c
		# cp <- as.numeric(crossprod(v[[2]], sigma %*% v[[2]]))
		# if (cp > 0) v[[2]] <- v[[2]] / sqrt(cp)
	
	if (debug) {
		check2 <- crossprod(c, v[[2]]) # debugging == sqrt(cp)
		cat(" ",check2)
	}

	## Update canonical vector in dimension 3
	c <- numeric(r)
	d <- matrix(0,r,n)
	for (k in 1:r) { 
		c[k] <- crossprod(v[[1]], a[,,k] %*% v[[2]])
		for (t in 1:n) 
			d[k,t] <- crossprod(v[[1]], b[,,k,t] %*% v[[2]])
	}
	v[[3]] <- unlist(optim1D.cor(c, d))
		# sigma <- tcrossprod(d) / n
		# v[[3]] <- tryCatch(solve(sigma, c), error = function(e) numeric(0))
		# if (length(v[[3]]) == 0) v[[3]] <- ginv(sigma) %*% c
		# cp <- as.numeric(crossprod(v[[3]], sigma %*% v[[3]]))
		# if (cp > 0) v[[3]] <- v[[3]] / sqrt(cp)

	## Calculate objective
	objective[it] <- crossprod(c, v[[3]]) # == sqrt(cp)
	if (debug) {
		cat(" ", objective[it])
	}
	
	## Check convergence 
	if (it > 1 && abs(objective[it]-objective[it-1]) <= 
	    	tol * max(1,objective[it-1])) break
}

return(v)
}


####################################
# Wrapper function for optimization
####################################


optim.cor <- function(v, a, b, maxit = 1000, tol = 1e-6, debug = TRUE)
{
if (length(v) == 1) {
	optim1D.cor(a, b)
} else if (length(v) == 2) {
	optim2D.cor(v, a, b, maxit, tol, debug)
} else {
	optim3D.cor(v, a, b, maxit, tol, debug)
}
}



