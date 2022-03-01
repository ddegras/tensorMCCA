##############################
# Initialization of canonical
# vectors with ones
##############################


init.v.ones <- function(x, type = c("norm", "var"), 
	scope = c("block", "global"), balance = TRUE)
{
## Check argument x 
test <- check.arguments(x)

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
dim.img <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
ndim.img <- sapply(dim.img, length)
n <- tail(dimx[[1]], 1)

## Scaling constraints
type <- match.arg(type)   # norm or variance constraints
scope <- match.arg(scope) # block or global constraints

## Set canonical vectors to 1's
v <- vector("list",m)
for (i in 1:m)
	v[[i]] <- lapply(dim.img[[i]], function(len) rep_len(1,len))

## Scale canonical vectors
v <- scale.v(v, x, type, scope, balance)

return(v)
}



########################
# Random initialization 
# of canonical vectors 
########################


init.v.random <- function(x, type = c("norm", "var"), 
	scope = c("block", "global"), balance = TRUE)
{
## Check argument x if required
test <- check.arguments(x)

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
dim.img <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
ndim.img <- sapply(dim.img, length)

## Scaling constraints
type <- match.arg(type)   # norm or variance constraints
scope <- match.arg(scope) # block or global constraints

## Initialize canonical vectors randomly
v <- vector("list",m)
for (i in 1:m)
	v[[i]] <- lapply(dim.img[[i]], function(len) runif(len))

## Scale canonical vectors
v <- scale.v(v, x, type, scope, balance)
return(v)
}



###########################
# SVD-based initialization
# of canonical vectors
###########################


init.v.svd <- function(x, type = c("norm", "var"), 
	scope = c("block", "global"), balance = TRUE)
{
## Check argument x if required
test <- check.arguments(x)

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
dim.img <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
ndim.img <- sapply(dim.img, length)
n <- tail(dimx[[1]], 1)
	
## Scaling constraints
type <- match.arg(type)   # norm or variance constraints
scope <- match.arg(scope) # block or global constraints

## Initialize canonical vectors to 1's
v <- vector("list",m)
for (i in 1:m)
	v[[i]] <- lapply(dim.img[[i]], function(len) rep_len(1,len))
	
## Update canonical vector in mode 1 by SVD 	
dim.k <- sapply(dim.img, "[", 1)
csum <- cumsum(c(0,dim.k))
compk <- matrix(0, csum[m+1], n)			
idx.rows <- lapply(1:m, 
	function(i) seq.int(csum[i] + 1, csum[i+1]))
for (i in 1:m) {
	compk[idx.rows[[i]],] <- if (ndim.img[i] == 1) {
		x[[i]] } else {
		apply(x[[i]], c(1, ndim.img[i] + 1), mean) }
}
vk <- if (min(dim(compk)) >= 3) {
	svds(compk, 1, 1, 0)$u } else { svd(compk, 1, 0)$u }
for (i in 1:m) 
	v[[i]][[1]] <- vk[idx.rows[[i]]]	
v <- scale.v(v, x, type, scope, balance)

## Update canonical vector in mode 2 by SVD 	
idx.i <- which(ndim.img >= 2)
if (length(idx.i) > 0) {
	dim.k <- sapply(dim.img[idx.i], "[", 2)
	csum <- cumsum(c(0,dim.k))
	compk <- matrix(0, tail(csum,1), n)			
	idx.rows <- lapply(seq_along(idx.i), 
		function(i) seq.int(csum[i] + 1, csum[i+1]))
	for (i in seq_along(idx.i)) {
		ii <- idx.i[i]
		iprod <- x[[ii]]
		dims <- dimx[[ii]]
		if (ndim.img[ii] == 3) 
			iprod <- apply(iprod, c(1,2,4), mean)
		dim(iprod) <- c(dims[1], dims[2] * n)
		iprod <- crossprod(v[[ii]][[1]], iprod)
		dim(iprod) <- c(dims[2], n)
		compk[idx.rows[[i]],] <- iprod
	}
	vk <- if (min(dim(compk)) >= 3) {
		svds(compk, 1, 1, 0)$u } else { svd(compk, 1, 0)$u }
	for (i in seq_along(idx.i)) 
		v[[idx.i[i]]][[2]] <- vk[idx.rows[[i]]]
	v <- scale.v(v, x, type, scope, balance)
}
	
## Update canonical vector in mode 3 by SVD 	
idx.i <- which(ndim.img == 3)
if (length(idx.i) > 0) {
	dim.k <- sapply(dim.img[idx.i], "[", 3)
	csum <- cumsum(c(0,dim.k))
	compk <- matrix(0, tail(csum,1), n)			
	idx.rows <- lapply(seq_along(idx.i), 
		function(i) seq.int(csum[i] + 1, csum[i+1]))
	for (i in seq_along(idx.i)) {
		ii <- idx.i[i]
		iprod <- x[[ii]]
		dims <- dimx[[ii]]	
		dim(iprod) <- c(dims[1], prod(dims[-1]))
		iprod <- crossprod(v[[ii]][[1]], iprod)
		dim(iprod) <- c(dims[2], dims[3] * n)
		iprod <- crossprod(v[[ii]][[2]], iprod)
		dim(iprod) <- c(dims[3], n)
		compk[idx.rows[[i]],] <- iprod
	}
	vk <- if (min(dim(compk)) >= 3) {
		svds(compk, 1, 1, 0)$u } else { svd(compk, 1, 0)$u }
	for (i in seq_along(idx.i)) 
		v[[idx.i[i]]][[3]] <- vk[idx.rows[[i]]]
	v <- scale.v(v, x, type, scope, balance)
}		

return(v)	
}



###########################
# SVD-based initialization
# of canonical vectors 
# (Quick Rank 1)
###########################


init.v.svd2 <- function(x, r, type = c("norm", "var"), 
	scope = c("block", "global"), balance = TRUE)
{
## Check argument x if required
test <- check.arguments(x)

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
dim.img <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
ndim.img <- sapply(dim.img, length)
n <- tail(dimx[[1]], 1)
	
## Scaling constraints
type <- match.arg(type)   # norm or variance constraints
scope <- match.arg(scope) # block or global constraints

## Initialize canonical vectors to 1's
v <- vector("list",m)

## Unfold data along last mode (individuals/objects) and concatenate
for (i in 1:m) 
	dim(x[[i]]) <- c(prod(dim.img[[i]]), n)
x <- do.call(x, rbind)

## Initial SVD
u <- svds(x, r, nu = r, nv = 0)$u

## Iteratively unfold singular vectors and recalculate SVD
for (i in 1:m) {
	for (k in 1:ndim.img[i]) {
		if (i == m && k == ndim.img[m]) {
			break
		} else {
			d1 <- dim.img[[i]][k]
			d2 <- length(u) / d1 / r
			dim(u) <- c(d1, d2, r)
			unew <- matrix(nrow = d2, ncol = r)
			last <- (i == m && k == ndim.img[m]-1)
			for (l in 1:r) {
				svdu <- svds(u[,,r], 1)
				v[[i]][[k]][,l] <- svdu$u
				unew[,l] <- svdu$v
				if (last) 
					v[[i]][[k+1]][,l] <- svdu$v
			}
			u <- unew
		}		
	}	
}

# ISSUE: THE VECTORS FOUND ARE NOT ORTHOGONAL 
	
## Scale initial canonical vectors as required	
v <- scale.v(v, x, type, scope, balance)

return(v)	
}