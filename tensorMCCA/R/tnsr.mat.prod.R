tnsr.mat.prod <- function(x, mat, modes = NULL)
{	
if (!is.list(mat)) mat <- list(mat)
if (length(mat) == 0) return(x)
if (is.null(dim(x))) x <- as.matrix(x)
d <- length(dim(x))
if (is.null(modes) && length(mat) == d) modes <- 1:d
stopifnot(length(mat) == length(modes))
modes <- as.integer(modes)
nmodes <- length(modes)
if (nmodes > 1 && any(diff(modes) < 0)) {
	ord <- order(modes)
	mat <- mat[ord]
	modes <- modes[ord]
}

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
	x <- mat[[k]] %*% x
	dim(x) <- c(nrow.mat[k], dimx[perm[-1]])
	if (k > 1) x <- aperm(x, perm)
}
x
}



tnsr.rk1.mat.prod <- function(v, mat, modes, transpose.mat = FALSE)
{
m <- length(v)
nmodes <- sapply(modes, length)
for (i in 1:m) {
	for (kk in nmodes[i]) {
		k <- modes[[i]][kk]
		v[[i]][[k]] <- if (transpose.mat) {
			crossprod(mat[[i]][[kk]], v[[i]][[k]])
		} else { mat[[i]][[kk]] %*% v[[i]][[k]] }
	}
}
v
}
