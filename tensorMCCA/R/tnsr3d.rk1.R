##############################
# Rank-1 tensor approximation 
# to order 3 tensor
##############################


tnsr3d.rk1 <- function(x, maxit = 100, tol = 1e-6, verbose = FALSE)
{
stopifnot(is.array(x) && length(dim(x)) == 3)
p <- dim(x)
v <- vector("list",3)
test <- all(d >= 3)
maxit <- as.integer(maxit)
stopifnot(maxit >= 0)

## Quick Rank 1 initialization
## A.k.a. interlaced computation of rank-1 HOSVD
v <- hosvd(x, 1)$u
# dim(x) <- c(d[1], d[2] * d[3])
# svd1 <- if (test) svds(x, k = 1) else svd(x, nu = 1, nv = 1)
# mat <- matrix(svd1$v, d[2], d[3])
# svd2 <- if (test) svds(mat, k = 1) else svd(mat, nu = 1, nv = 1)
# s <- (svd1$d[1] * svd2$d[1])^(1/3)
# v[[1]] <- svd1$u * s
# v[[2]] <- svd2$u * s
# v[[3]] <- svd2$v * s

iters <- if (maxit > 0) {1:maxit} else NULL
for (it in iters) {
	v.old <- v
	
	## Update v1
	dim(x) <- c(d[1] * d[2], d[3])
	xv <- x %*% v[[3]]
	dim(xv) <- c(d[1], d[2])
	v[[1]] <- xv %*% v[[2]] 
	v[[1]] <- v[[1]] / sqrt(sum(v[[1]]^2))

	## Update v2
	v[[2]] <- crossprod(xv, v[[1]]) 
	v[[2]] <- v[[2]] / sqrt(sum(v[[2]]^2))
	 	
	## Update v3
	dim(x) <- c(d[1], d[2] * d[3])
	xv <- crossprod(x, v[[1]])
	dim(xv) <- d[2:3]
	v[[3]] <- crossprod(xv, v[[2]]) 
	nrm <- sqrt(sum(v[[3]]^2))
			
	## Check progress
	e <- kronecker(v[[3]], kronecker(v[[2]], v[[1]])) - 
		kronecker(v.old[[3]], kronecker(v.old[[2]], v.old[[1]]))
	if (sqrt(sum(e^2)) <= tol * nrm) break
}

## Balance norms 
s <- nrm^(c(1/3, 1/3, -2/3))
for (i in 1:3) v[[i]] <- v[[i]] * s[i] 

v	
	
}







########
# HOSVD
########


hosvd <- function(x, r = NULL)
{
stopifnot(is.numeric(x))
if (is.vector(x)) {
	d <- 1L
	p <- length(x) 
} else {
	d <- length(p) 
	p <- dim(x)
}
if (is.null(r)) {
	r <- if (is.vector(x)) 1L else p 
} else {
	stopifnot(all(r > 0))
	r <- rep_len(as.integer(r), d)
	r <- pmin(r, p)
}

## Trivial cases 
if (d == 1) {
	nrm <- sqrt(sum(x^2))
	if (nrm == 0) {
		return(list(u = numeric(p), core = 0))
	} else {
		return(list(u = x/nrm, core = nrm))
	}
}
if (d == 2) {
	RSpectra.flag <- all(p >= 3)
	svdx <- if (RSpectra.flag) {
		svds(x, max(r), r[1], r[2])	
	} else svd(x, r[1], r[2])
	return(list(factors = list(svdx$u, svdx$v), 
		core = diag(svdx$d, r[1], r[2])))
}

## Outputs
core <- x
u <- vector("list", d)

for (k in 1:d) {
	## Mode-k flattening
	if (k > 1) {
		perm <- 1:d
		perm[c(1,k)] <- perm[c(k,1)]
		core <- aperm(core, perm)
	}
	dimc <- dim(core)
	dim(core) <- c(p[k], prod(dimc[-1]))
	
	## SVD
	RSpectra.flag <- (min(dim(core)) > max(2, r[k]))
	svdk <- if (RSpectra.flag) {
		svds(core, r[k])
	} else { svd(core, nu = r[k], nv = r[k]) }
	u[[k]] <- svdk$u
	
	## Reshape and permute core
	core <- svdk$d[1:(r[k])] * t(svdk$v)
	dim(core) <- c(r[k], dimc[-1])
	if (k > 1) core <- aperm(core, perm)	
}	

list(factors = u, core = core)
	
}






##############################
# Rank-1 tensor approximation 
# to general tensor
##############################

# Not needed yet (will be needed to extend package to tensors of order 4 and higher)


# tnsr.rk1 <- function(x, maxit = 100, tol = 1e-6, verbose = FALSE)
# {
# ## Trivial case
# if (is.vector(x)) return(x)

# ## Dimensions of order-d tensor with d > 1
# stopifnot(is.array(x))
# p <- dim(x) # tensor dimensions
# pp <- prod(p) # total number of elements 
# d <- length(p) # tensor order
# v <- vector("list", d)

# ## Rank-1 HOSVD
# hosvdx <- hosvd(x, 1)


# for (it in 1:maxit) {
	# v.old <- v
	
	# ## Update v1
	# dim(x) <- c(d[1] * d[2], d[3])
	# xv <- x %*% v[[3]]
	# dim(xv) <- c(d[1], d[2])
	# v[[1]] <- xv %*% v[[2]] 
	# v[[1]] <- v[[1]] / sqrt(sum(v[[1]]^2))

	# ## Update v2
	# v[[2]] <- crossprod(xv, v[[1]]) 
	# v[[2]] <- v[[2]] / sqrt(sum(v[[2]]^2))
	 	
	# ## Update v3
	# dim(x) <- c(d[1], d[2] * d[3])
	# xv <- crossprod(x, v[[1]])
	# dim(xv) <- d[2:3]
	# v[[3]] <- crossprod(xv, v[[2]]) 
	# nrm <- sqrt(sum(v[[3]]^2))
			
	# ## Check progress
	# e <- kronecker(v[[3]], kronecker(v[[2]], v[[1]])) - 
		# kronecker(v.old[[3]], kronecker(v.old[[2]], v.old[[1]]))
	# if (sqrt(sum(e^2)) <= tol * nrm) break
# }

# ## Balance norms 
# s <- nrm^(c(1/3, 1/3, -2/3))
# for (i in 1:3) v[[i]] <- v[[i]] * s[i] 

# v	
	
# }


