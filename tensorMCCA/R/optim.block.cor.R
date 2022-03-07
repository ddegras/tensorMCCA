
# FUNCTIONS FOR MAXIMIZING THE SUM OF CORRELATIONS 
# IN MCCA ONE DATA BLOCK AT A TIME




##########################################
# Maximize linear form v'a subject to
# (1/n) v'BB'v = 1 where B has n columns  
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
# constraint (1/n) sum_t (v1' Bt v2)^2 = 1
##########################################


## Inputs: 
## v:	initial solution (list of 2 vectors)
## a:	objective (2D matrix)
## b:	constraints (3D array)

optim2D.cor <- function(v, a, b, maxit = 1000, tol = 1e-6)
{
	
## Data dimensions
dimb <- dim(b)
p <- dimb[1:2]
n <- dimb[3]

## MAIN LOOP
objective <- numeric(maxit)
for (it in 1:maxit) {

	## Update canonical vector in dimension 1 
	c <- a %*% v[[2]]
	d <- matrix(0, p[1], n)
	for (t in 1:n) d[,t] <- b[,,t] %*% v[[2]]
	v[[1]] <- unlist(optim1D.cor(c, d))
		
	## Update canonical vector in dimension 2
	c <- crossprod(a, v[[1]])
	d <- matrix(0, p[2], n)
	for (t in 1:n) d[,t] <- crossprod(b[,,t], v[[1]])
	v[[2]] <- unlist(optim1D.cor(c, d))
	
	## Calculate objective
	objective[it] <- crossprod(c, v[[2]])
	
	## Check convergence 
	if (it > 1 && abs(objective[it]-objective[it-1]) <= 
	    	tol * max(1,objective[it-1])) break
}

return(v)
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
objective <- numeric(maxit)

for (it in 1:maxit) {

	## Update canonical vector in dimension 1 
	c <- numeric(p[1])
	d <- matrix(0, p[1], n)
	for (i in 1:p[1]) { 
		c[i] <- crossprod(v[[2]], a[i,,] %*% v[[3]])
		for (t in 1:n) 
			d[i,t] <- crossprod(v[[2]], b[i,,,t] %*% v[[3]])
	}
	v[[1]] <- unlist(optim1D.cor(c, d))
	
	## Update canonical vector in dimension 2
	c <- numeric(p[2])
	d <- matrix(0, p[2], n)
	for (j in 1:p[2]) { 
		c[j] <- crossprod(v[[1]], a[,j,] %*% v[[3]])
		for (t in 1:n) 
			d[j,t] <- crossprod(v[[1]], b[,j,,t] %*% v[[3]])
	}
	v[[2]] <- unlist(optim1D.cor(c, d))

	## Update canonical vector in dimension 3
	c <- numeric(p[3])
	d <- matrix(0, p[3], n)
	for (k in 1:p[3]) { 
		c[k] <- crossprod(v[[1]], a[,,k] %*% v[[2]])
		for (t in 1:n) 
			d[k,t] <- crossprod(v[[1]], b[,,k,t] %*% v[[2]])
	}
	v[[3]] <- unlist(optim1D.cor(c, d))

	## Calculate objective
	objective[it] <- crossprod(c, v[[3]]) # == sqrt(cp)
	
	## Check convergence 
	if (it > 1 && abs(objective[it]-objective[it-1]) <= 
	    	tol * max(1,objective[it-1])) break
}

return(v)
}


####################################
# Wrapper function for optimization
####################################


optim.block.cor <- function(v, a, b, maxit, tol)
{
if (length(v) == 1) {
	optim1D.cor(a, b)
} else if (length(v) == 2) {
	optim2D.cor(v, a, b, maxit, tol)
} else {
	optim3D.cor(v, a, b, maxit, tol)
}
}



