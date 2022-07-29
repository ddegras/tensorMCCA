match.list.dims <- function(x, dims)
{
stopifnot(is.list(x))
stopifnot((is.numeric(dims) || 
	(is.list(dims) && length(x) == length(dims)))
if (is.numeric(dims)) {
	dims <- as.integer(dims)
	dimx <- sapply(x, length, USE.NAMES = FALSE)
	if (identical(dimx, dims)) return(x)
	dims0 <- dims[dims > 1]
	if (!identical(dimx, dims0))
		stop("'x' and 'dims' have incompatible dimensions")
	out <- vector("list", length(dims))
	out[dims == 1] <- replicate(sum(dims == 1), 1, FALSE)
	out[dims != 1] <- x
} else {
	m <- length(x)
	out <- vector("list", m)
	for (i in 1:m) {
		dims <- as.integer(dims[[i]])
		dimx <- sapply(x[[i]], length, USE.NAMES = FALSE)
	if (identical(dimx, dims)) {
		out[[i]] <- x[[i]]
		next
	}
	dims0 <- dims[dims > 1]
	if (!identical(dimx, dims0)) 
		stop(paste("x[[", i, "]] and dims[[", i, 
			"]] have incompatible dimensions"))
	out[[i]] <- vector("list", length(dims))
	out[[i]][dims == 1] <- replicate(sum(dims == 1), 1, FALSE)
	out[[i]][dims != 1] <- x
	}
}
out	
}