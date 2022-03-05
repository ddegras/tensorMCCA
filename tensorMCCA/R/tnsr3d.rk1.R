tnsr3d.rk1 <- function(x, maxit = 100, tol = 1e-6, verbose = FALSE)
{
stopifnot(is.array(x) && length(dim(x)) == 3)
d <- dim(x)
v <- vector("list",3)
test <- all(d >= 3)

## Quick Rank 1 initialization
dim(x) <- c(d[1], d[2] * d[3])
svd1 <- if (test) svds(x, k = 1) else svd(x, nu = 1, nv = 1)
mat <- matrix(svd1$v, d[2], d[3])
svd2 <- if (test) svds(mat, k = 1) else svd(mat, nu = 1, nv = 1)
s <- (svd1$d[1] * svd2$d[1])^(1/3)
v[[1]] <- svd1$u * s
v[[2]] <- svd2$u * s
v[[3]] <- svd2$v * s

for (it in 1:maxit) {
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



# tensor.norm.3d.v2 <- function(x, maxit = 100, tol = 1e-6, verbose = FALSE)
# {
# stopifnot(is.array(x) && length(dim(x)) == 3)
# dimx <- dim(x)
# v <- vector("list",3)
# svdfun <- if (all(dimx >= 3)) {
	# function(x) RSpectra::svds(x, k = 1) 
	# } else { function(x) svd(x, nu = 1, nv = 1) }
# ## Quick Rank 1 initialization
# dim(x) <- c(dimx[1], prod(dimx[2:3])) # mode-1 unfolding
# svdx1 <- svdfun(x)
# v[[1]] <- svdx1$u
# xv <- matrix(svdx1$v, dimx[2], dimx[3])
# svdxv <- svdfun(xv)
# v[[2]] <- svdxv$u; v[[3]] <- svdxv$v
# xv <- matrix(crossprod(v[[1]], x), dimx[2], dimx[3])
# xnorm <- as.numeric(crossprod(v[[2]], xv %*% v[[3]]))
# if (verbose) print(xnorm)
# for (it in 1:maxit) {
	# xnorm.old <- xnorm
	# ## Update v1 and v2
	# dim(x) <- c(prod(dimx[1:2]), dimx[3])
	# xv <- matrix(x %*% v[[3]], dimx[1], dimx[2])
	# svdxv <- svdfun(xv)
	# v[[1]] <- svdxv$u; v[[2]] <- svdxv$v
	
	# ## Update v2 and v3
	# dim(x) <- c(dimx[1], prod(dimx[2:3]))
	# xv <- matrix(crossprod(x, v[[1]]), dimx[2], dimx[3])
	# svdxv <- svdfun(xv)
	# v[[2]] <- svdxv$u; v[[3]] <- svdxv$v
	
	# xnorm <- svdxv$d[1]
	# if (verbose) print(xnorm)
	
	# ## Check progress	
	# if (xnorm  < (1 + tol) * xnorm.old) break
# }

# xnorm	
	
# }
