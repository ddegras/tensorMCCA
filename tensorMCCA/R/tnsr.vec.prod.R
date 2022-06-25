tnsr.vec.prod <- function(x, v, modes = NULL)
{	
p <- dim(x) 
d <- length(dim.tnsr)
nmodes <- length(modes)
ord <- order(modes)
if (any(diff(ord) < 0)) {
	v <- v[ord]
	modes <- modes[ord]
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
