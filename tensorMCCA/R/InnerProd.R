
######################################
# Inner product between rank 1 tensor 
# and general tensor (1D, 2D or 3D)
######################################

# x: general tensor of order d (d = 1, 2, or 3)
# v: rank-1 tensor specified via list of length d

InnerProd <- function(x, v)
{
	dimx <- dim(x)
	if (is.null(dimx)) {
		as.numeric(crossprod(x, v[[1]]))
	} else if (length(dimx) == 2) {
		as.numeric(crossprod(v[[1]], x %*% v[[2]]))
	} else if (length(dimx) == 3) {
		dim(x) <- c(dimx[1], dimx[2] * dimx[3])
		iprod <- crossprod(v[[1]], x)
		dim(iprod) <- dimx[2:3]
		as.numeric(crossprod(v[[2]], iprod %*% v[[3]]))
	} else NULL
}



##############################################
# Function to calculate tensor-vector product 
##############################################

tnsr.vec.prod <- function(tnsr, vec, mode)
{
	dims <- dim(tnsr) 
	d <- length(dims)
	if (mode == 1) {
		dim(tnsr) <- c(dims[1], prod(dims[-1]))
		tnsr <- crossprod(vec, tnsr)
		dim(tnsr) <- dims[-1]
	} else if (mode == d) {
		dim(tnsr) <- c(prod(dims[-d]), dims[d])
		tnsr <- tnsr %*% vec
		dim(tnsr) <- dims[-d]
	} else {
		perm <- c(1:(mode-1), (mode+1):d, mode)
		tnsr <- aperm(tnsr, perm)
		dim(tnsr) <- c(prod(dims[-mode]), mode)
		tnsr <- tnsr %*% vec
		dim(tnsr) <- dims[-mode]
	}
	tnsr
}

