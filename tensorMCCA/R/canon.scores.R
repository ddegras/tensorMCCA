canon.scores <- function(x, v, check.args = TRUE)
{
if (!is.list(x)) x <- list(x)
if (!is.list(v)) v <- list(v)
if (check.args) {
	stopifnot(is.vector(v) || is.matrix(v))
	test1 <- (is.vector(v) && length(v) == length(x))
	test2 <- (is.matrix(v) && nrow(v) == length(x))
	if (!(test1 || test2)) {
		stop(paste("'v' must either be a vector of type 'list'",
		"having same length as 'x' or a matrix of type 'list'",
		"whose number of rows equals the length of 'x'")) }
}
m <- length(x)
dimx <- lapply(x, dim)
n <- tail(dimx[[1]], 1) 
d <- sapply(dimx, length) - 1  
r <- NCOL(v)
if (r == 1) {
	score <- matrix(nrow = n, ncol = m)
	for (i in 1:m) 
		score[, i] <- tnsr.vec.prod(x[[i]], v[[i]], 1:d[i])
} else {
	score <- array(dim = c(n, m, r))
	for (i in 1:m) {
		for (l in 1:r) {
			score[, i, l] <- tnsr.vec.prod(x[[i]], v[[i,l]], 1:d[i])
		}
	}
}
score	
}









# image.scores <- function(x, v)
# {
# m <- length(x) # number of datasets
# dimx <- lapply(x, dim)
# n <- tail(dimx[[1]], 1) # number of instances per dataset
# d <- sapply(dimx, length) - 1 # number of dimensions in images/features 
# y <- matrix(0, n, m) # image/individual scores on canonical components
# for (i in 1:m) { 		
	# iprod <- x[[i]]
	# pi <- dimx[[i]]
	# if (d[i] == 1) { # 1D case
		# y[,i] <- crossprod(iprod, v[[i]][[1]]) 
	# } else if (d[i] == 2) { # 2D case
		# dim(iprod) <- c(pi[1], pi[2] * n)
		# iprod <- crossprod(v[[i]][[1]], iprod)
		# dim(iprod) <- c(pi[2], n) 
		# y[,i] <- crossprod(iprod, v[[i]][[2]]) 
	# } else if (d[i] == 3) { # 3D case
		# dim(iprod) <- c(pi[1], prod(pi[2:4]))
		# iprod <- crossprod(v[[i]][[1]], iprod)
		# dim(iprod) <- c(pi[2], pi[3] * n)
		# iprod <- crossprod(v[[i]][[2]], iprod)
		# dim(iprod) <- c(pi[3], n) 
		# y[,i] <- crossprod(iprod, v[[i]][[3]])
	# } else warning("Function not yet supported for tensors of order 4+")
# }
# y
# }
