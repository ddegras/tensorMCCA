
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



###############################################
# Function to calculate tensor-vector products 
###############################################

# tvec.prod generalizes tnsr.vec.prod and is as fast
# tvec.prod2 is never much faster than tvec.prod, sometimes markedly slower


## Product of tensor with vector in single mode
# tnsr.vec.prod <- function(tnsr, vec, mode)
# {
	# dims <- dim(tnsr) 
	# d <- length(dims)
	# if (mode == 1) {
		# out <- tnsr
		# dim(out) <- c(dims[1], prod(dims[-1]))
		# out <- crossprod(vec, out)
	# } else if (mode == d) {
		# out <- tnsr
		# dim(out) <- c(prod(dims[-d]), dims[d])
		# out <- out %*% vec
	# } else {
		# perm <- c((1:d)[-mode], mode)
		# out <- aperm(tnsr, perm)
		# dim(out) <- c(prod(dims[-mode]), dims[mode])
		# out <- out %*% vec
	# }
	# dim(out) <- dims[-mode]
	# out
# }


## Product of tensor with vectors in multiple modes
tvec.prod <- function(tnsr, vecs, modes = NULL)
{
	## Check argument types
	stopifnot(is.array(tnsr))
	stopifnot(is.numeric(vecs) || is.list(vecs))
	stopifnot(is.null(modes) || is.numeric(modes))
	modes <- as.integer(modes)
	if (is.numeric(vecs)) vecs <- list(vecs)
	
	## Argument dimensions
	dim.tnsr <- dim(tnsr) 
	tnsr.order <- length(dim.tnsr)
	nvecs <- length(vecs)
	nmodes <- length(modes)
	
	## Check compatibility of arguments
	if (nmodes == 0) {
		if (nvecs == tnsr.order) {
			modes <- 1:tnsr.order
		} else {
			stop(paste("Please specify the modes in which", 
			"to perform the multiplications"))
		}
	} else {
		if (nmodes != nvecs)
			stop("'modes' and 'vecs' must have equal lengths")
		if (any(modes <= 0 | modes > tnsr.order)) 
			stop(paste("'modes' must only contain integers",
			"between 1 and the number of dimensions of 'tnsr'"))
		if (nmodes > length(unique(modes)))
			stop("The values in 'modes' must be unique")
	}
			
	## Calculate tensor-vector products
	ord <- order(modes)
	if (any(diff(ord) < 0)) {
		vecs <- vecs[ord]
		modes <- modes[ord]
	}
	first.modes <- identical(modes, 1:nmodes)
	last.modes <- identical(modes, tail(1:tnsr.order, nmodes))
	if (first.modes) {
		out <- tnsr
		for (kk in 1:nmodes) {
			k <- modes[kk]
			nc <- if (k < tnsr.order) prod(dim.tnsr[-(1:k)]) else 1L
			dim(out) <- c(dim.tnsr[k], nc)
			out <- crossprod(vecs[[kk]], out)
		}		
	} else if (last.modes) {
		out <- tnsr
		for (kk in 1:nmodes) {
			k <- modes[nmodes - kk + 1L]
			dim(out) <- c(prod(dim.tnsr[1:(k - 1L)]), dim.tnsr[k])
			out <- out %*% vecs[[nmodes - kk + 1L]]
		}
	} else {
		perm <- c((1:tnsr.order)[-modes], modes)
		out <- aperm(tnsr, perm)
		dim.out <- dim.tnsr[perm]
		for (kk in 1:nmodes) {
			k <- modes[nmodes - kk + 1L]
			dim(out) <- c(prod(dim.out[1:(tnsr.order - kk)]), 
				dim.tnsr[k])
			out <- out %*% vecs[[nmodes - kk + 1L]]
		}
	}
	dim(out) <- if (nmodes == tnsr.order) {
		NULL } else if (nmodes == tnsr.order - 1L) {
		c(length(out), 1L) } else { 
		dim.tnsr[-modes] }
	out
}

# tvec.prod2 <- function(tnsr, vecs, modes = NULL)
# {
	# ## Check argument types
	# stopifnot(is.array(tnsr))
	# stopifnot(is.numeric(vecs) || is.list(vecs))
	# stopifnot(is.null(modes) || is.numeric(modes))
	# modes <- as.integer(modes)
	# if (is.numeric(vecs)) vecs <- list(vecs)
	
	# ## Argument dimensions
	# dim.tnsr <- dim(tnsr) 
	# tnsr.order <- length(dim.tnsr)
	# nvecs <- length(vecs)
	# nmodes <- length(modes)
	
	# ## Check compatibility of arguments
	# if (nmodes == 0) {
		# if (nvecs == tnsr.order) {
			# modes <- 1:tnsr.order
		# } else {
			# stop(paste("Please specify the modes in which", 
			# "to perform the multiplications"))
		# }
	# } else {
		# if (nmodes != nvecs)
			# stop("'modes' and 'vecs' must have equal lengths")
		# if (any(modes <= 0 | modes > tnsr.order)) 
			# stop(paste("'modes' must only contain integers",
			# "between 1 and the number of dimensions of 'tnsr'"))
		# if (nmodes > length(unique(modes)))
			# stop("The values in 'modes' must be unique")
	# }
			
	# ## Calculate tensor-vector products
	# ord <- order(modes)
	# if (any(diff(ord) < 0)) {
		# vecs <- vecs[ord]
		# modes <- modes[ord]
	# }
	# if (identical(modes, 1:tnsr.order)) {
		# vecs <- Reduce(kronecker, rev(vecs))
		# dim(vecs) <- NULL
		# out <- sum(tnsr * vecs)
	# } else if (identical(modes, 1:nmodes)) {
		# out <- tnsr
		# dim(out) <- c(prod(dim.tnsr[1:nmodes]),
			# prod(dim.tnsr[-(1:nmodes)]))
		# vecs <- Reduce(kronecker, rev(vecs[modes]))
		# dim(vecs) <- NULL
		# out <- crossprod(vecs, out)		
	# } else if (identical(modes, (tnsr.order-nmodes+1):tnsr.order)) {
		# out <- tnsr
		# dim(out) <- c(prod(dim.tnsr[-modes]), prod(dim.tnsr[modes]))
		# vecs <- Reduce(kronecker, rev(vecs[modes]))
		# dim(vecs) <- NULL
		# out <- out %*% vecs
	# } else {
		# perm <- c((1:tnsr.order)[-modes], modes)
		# out <- aperm(tnsr, perm)
		# dim(out) <- c(prod(dim.tnsr[-modes]), prod(dim.tnsr[modes]))
		# vecs <- Reduce(kronecker, rev(vecs[modes]))
		# dim(vecs) <- NULL
		# out <- out %*% vecs
	# }
	# dim(out) <- if (nmodes >= tnsr.order - 1) NULL else dim.tnsr[-modes]
	# out
# }



tmat.prod <- function(tnsr, mats, modes = NULL)
{
## Check argument types
stopifnot(is.array(tnsr))
stopifnot(is.numeric(mats) || is.list(mats))
stopifnot(is.null(modes) || is.numeric(modes))
modes <- as.integer(modes)
if (is.numeric(mats)) mats <- list(mats)
isvec <- sapply(mats, is.vector)
ismat <- sapply(mats, is.matrix)
if (!all(isvec | ismat))
	stop(paste("'mats' must be a list of numerical", 
	"vectors or matrices"))
if (any(isvec)) mats[isvec] <- lapply(mats[isvec], t)

## Argument dimensions
dim.tnsr <- dim(tnsr) 
tnsr.order <- length(dim.tnsr)
nmats <- length(mats)
nmodes <- length(modes)
nr <- sapply(mats, nrow)
nc <- sapply(mats, ncol)

## Check compatibility of arguments
if (nmodes == 0) {
	if (nmats == tnsr.order) {
		modes <- 1:tnsr.order
	} else {
		stop(paste("Please specify the modes in which", 
		"to perform the multiplications"))
	}
} else {
	if (nmodes != nmats)
		stop("'mats' and 'modes' must be of equal length")
	if (any(modes <= 0 | modes > tnsr.order)) 
		stop(paste("'modes' must only contain integers",
		"between 1 and the number of modes",
		"(i.e., dimensions) of 'tnsr'"))
	if (nmodes > length(unique(modes)))
		stop("The values in 'modes' must be unique")
	if (any(nc != dim.tnsr[modes]))
		stop(paste("The numbers of columns in 'mats'", 
		"do not match the corresponding dimensions in 'tnsr'"))
}
		
## Calculate tensor-matrix products
ord <- order(modes)
if (any(diff(ord) < 0)) {
	mats <- mats[ord]
	modes <- modes[ord]
}
out <- tnsr	
for (kk in 1:nmodes) {
	dim.out <- dim(out)
	k <- modes[kk]
	if (k > 1) {
		perm <- 1:tnsr.order
		perm[c(1,k)] <- c(k,1)
		out <- aperm(out, perm)
	}
	dim(out) <- c(dim.tnsr[k], prod(dim.out) / dim.tnsr[k])
	out <- mats[[kk]] %*% out
	dim(out) <- c(nr[kk], dim.out[perm[-1]])
	out <- aperm(out, perm)
}
out
}