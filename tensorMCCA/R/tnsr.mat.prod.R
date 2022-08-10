tnsr.mat.prod <- function(x, mat, modes = NULL)
{	
stopifnot(is.numeric(x))
if (!is.list(mat)) mat <- list(mat)
if (length(mat) == 0) return(x)
if (is.null(dim(x))) x <- as.matrix(x)
d <- length(dim(x))
if (is.null(modes) && length(mat) == d) modes <- 1:d
stopifnot(length(mat) == length(modes))
modes <- as.integer(modes)
nmodes <- length(modes)

if (d == 2) {
	if (identical(modes, 1L)) {
		return(mat[[1]] %*% x) 
	}
	if (identical(modes, 2L)) {
		return(tcrossprod(x, mat[[1]])) 
	}
	return(tcrossprod(mat[[1]] %*% x, mat[[2]]))
}

nrow.mat <- sapply(mat, nrow)
for (kk in 1:nmodes) {
	k <- modes[kk]
	perm <- 1:d
	perm[c(1,k)] <- c(k,1)
	dimx <- dim(x)
	if (k > 1) x <- aperm(x, perm)
	dim(x) <- c(dimx[k], prod(dimx[-k]))
	x <- mat[[kk]] %*% x
	dim(x) <- c(nrow.mat[kk], dimx[perm[-1]])
	if (k > 1) x <- aperm(x, perm)
}
x
}



tnsr.rk1.mat.prod <- function(v, mat, modes, transpose.mat = FALSE)
{
if (!is.list(v)) v <- list(as.numeric(v))
test.v.numeric <- all(sapply(v, is.numeric)) 
test.v.list <- all(sapply(v, is.list))
if (!(test.v.numeric || test.v.list))
	stop("'v' must be a list of vectors or a list of lists (of vectors)")

if (test.v.numeric) {
	modes <- as.integer(modes)
	if (!is.list(mat)) mat <- list(mat)
	stopifnot(length(modes) == length(mat))
	nmodes <- length(modes)
	for (kk in 1:nmodes) {
		k <- modes[kk]
		v[[k]] <- if (transpose.mat) {
			crossprod(mat[[kk]], v[[k]])
		} else { mat[[kk]] %*% v[[k]] }
	}
	v <- lapply(v, as.numeric)
	return(v)
}

m <- length(v)
test.mat.modes.list <- (is.list(mat) && is.list(modes))
test.mat.modes.length <- (length(mat) == m && length(modes) == m)
if (!(test.mat.modes.list && test.mat.modes.length))
	stop(paste("If 'v' is a list of lists, then 'mat' and 'modes'",
	"should be lists of same length as 'v'."))
nmodes <- sapply(modes, length)
for (i in 1:m) {
	for (kk in 1:nmodes[i]) {
		k <- modes[[i]][kk]
		# if (length(v[[i]][[k]]) == 1) browser() 
		v[[i]][[k]] <- if (transpose.mat) {
			crossprod(mat[[i]][[kk]], v[[i]][[k]])
		} else { mat[[i]][[kk]] %*% v[[i]][[k]] }
	}
	v[[i]] <- lapply(v[[i]], as.numeric)
}
v
}
