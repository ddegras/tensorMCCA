tnsr.vec.prod <- function(x, v, modes)
{
if (length(v) == 0) return(x)
p <- dimfun(x) 
d <- length(p)
if (!is.list(v)) v <- list(v)
if (missing(modes) && length(v) == d) modes <- 1:d
stopifnot(length(v) == length(modes))
modes <- as.integer(modes)
nmodes <- length(modes)
if (nmodes > 1 && any(diff(modes) < 0)) {
	ord <- order(modes)
	v <- v[ord]
	modes <- modes[ord]
}
if (d == 1) { 
	out <- sum(x * v[[1]])
	return(out) 
}
if (d == 2) {
	out <- if (identical(modes, 1L)) {
		crossprod(x, v[[1]])
	} else if (identical(modes, 2L)) {
		x %*% v[[1]]
	} else {
		crossprod(v[[1]], x %*% v[[2]])
	}
	dim(out) <- NULL
	return(out)
}
first.modes <- identical(modes, 1:nmodes)
last.modes <- identical(modes, tail(1:d, nmodes))
if (first.modes) {
	out <- x
	for (kk in 1:nmodes) {
		k <- modes[kk]
		dim(out) <- if (k < d) {
			c(p[k], prod(p[-(1:k)]))
		} else NULL
		out <- crossprod(v[[kk]], out)
	}		
} else if (last.modes) {
	out <- x
	for (kk in 1:nmodes) {
		k <- modes[nmodes - kk + 1L]
		dim(out) <- c(prod(p[1:(k - 1L)]), p[k])
		out <- out %*% v[[nmodes - kk + 1L]]
	}
} else {
	perm <- c((1:d)[-modes], modes)
	out <- aperm(x, perm)
	dim.out <- p[perm]
	for (kk in 1:nmodes) {
		k <- modes[nmodes - kk + 1L]
		dim(out) <- c(prod(dim.out[1:(d - kk)]), p[k])
		out <- out %*% v[[nmodes - kk + 1L]]
	}
}
dim(out) <- if (nmodes == d) {
	NULL } else if (nmodes == (d - 1L)) {
	c(length(out), 1L) } else { 
	p[-modes] }
out
}
