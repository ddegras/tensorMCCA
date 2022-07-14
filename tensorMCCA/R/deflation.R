

###########################################
# General function to deflate data tensors
###########################################


deflate.x <- function(x, v = NULL, score = NULL, 
	scope = c("block", "global"), check.args = TRUE)
{
scope <- match.arg(scope)
	
if (check.args) {
	test <- check.arguments(x, v) 
	m <- length(x)
	n <- tail(dim(x[[1]]), 1)
	if (is.null(v)  == is.null(score)) 
		stop(paste("Please specify exactly one of the arguments",
		"'v' and 'score' and leave the other set to NULL"))
	if (!is.null(score)) {
		stopifnot(is.vector(score) || is.matrix(score))
		test1 <- (is.vector(score) && length(score) == m)
		test2 <- (is.matrix(score) && all(dim(score) == c(n, m)))
		if (!(test1 || test2))
			stop(paste("If specified, 'score' should be either", 
				"a vector of length n or a matrix of dimensions",
				"(n,m)\nwhere m is the number of data tensors",
				"(length of 'x') and n is the number of dimensions",
				"in the last mode of the data tensors."))
	}
}	
	
## Data dimensions
m <- length(x)
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L
p <- lapply(1:m, function(i) dimx[[i]][1:d[i]])
n <- tail(dimx[[1]], 1)
eps <- 1e-14

if (!is.null(score)) {
	if (NCOL(score) != m) score <- matrix(score, n, m)
	score <- sweep(score, 2, colMeans(score), "-")
	nrm <- sqrt(colSums(score^2))	
	nzero <- which(nrm > eps)
	if (length(nzero) == 0) return(x)
	for (i in nzero) {
		dim(x[[i]]) <- c(prod(p[[i]]), n)
		x[[i]] <- x[[i]] - tcrossprod(x[[i]] %*% score[,i], 
			score[,i] / nrm[i]^2)
		dim(x[[i]]) <- dimx[[i]]
	}
}

if (!is.null(v) && scope == "block") {
	for (i in 1:m) {
		for (k in 1:d[i]) {
			if (is.null(v[[i]][[k]])) next
			## Basis of orthogonal space
			Q <- qr.Q(qr(v[[i]][[k]])) 
			## Permute data tensor (exchange modes 1 and k)
			if (k > 1) {
				perm <- 1:(d[i]+1)
				perm[c(1,k)] <- perm[c(k,1)]
				x[[i]] <- aperm(x[[i]], perm)
			}
			## Unfold permuted tensor along mode 1 
			## (original along mode k)
			dim(x[[i]]) <- c(dimx[[i]][k], prod(dimx[[i]][-k])) 
			## Orthogonal projection
			x[[i]] <- x[[i]] - Q %*% crossprod(Q, x[[i]]) 
			## Reshape x(i) to tensor
			dim(x[[i]]) <- if (k == 1) dimx[[i]] else dimx[[i]][perm]
			## Permute modes back
			if (k > 1) 
				x[[i]] <- aperm(x[[i]], perm)
		}
	}
}

if (!is.null(v) && scope == "global") {
	nmodes <- sapply(ortho.mode, length)
	nvecs <- prod(nmodes)
	fn <- function(v, i) 
	for (ii in 1:nvecs) {
		k <- arrayInd(ii, nmodes)	
		vv <- unlist(mapply("[[", v, k, SIMPLIFY = FALSE))
		for (i in 1:m) {
			
			
		}
		for (k in 1:d[i]) {
			if (is.null(v[[i]][[k]])) next
			## Basis of orthogonal space
			Q <- qr.Q(qr(v[[i]][[k]])) 
			## Permute data tensor (exchange modes 1 and k)
			if (k > 1) {
				perm <- 1:(d[i]+1)
				perm[c(1,k)] <- perm[c(k,1)]
				x[[i]] <- aperm(x[[i]], perm)
			}
			## Unfold permuted tensor along mode 1 
			## (original along mode k)
			dim(x[[i]]) <- c(dimx[[i]][k], prod(dimx[[i]][-k])) 
			## Orthogonal projection
			x[[i]] <- x[[i]] - Q %*% crossprod(Q, x[[i]]) 
			## Reshape x(i) to tensor
			dim(x[[i]]) <- if (k == 1) dimx[[i]] else dimx[[i]][perm]
			## Permute modes back
			if (k > 1) 
				x[[i]] <- aperm(x[[i]], perm)
		}
	}
}

x
}




################################################
# General function to deflate canonical vectors
# with respect to other canonical vectors
################################################


# Inputs
# v:	 list of lists of vectors or matrices (length m)
# vprev: list of lists of vectors or matrices (length m)

# deflate.v <- function(v, vprev)
# {
# stopifnot(is.list(v) && is.list(vprev))
# stopifnot(length(v) == length(vprev))
# m <- length(v)
# for (i in 1:m) {
	# stopifnot(length(v[[i]]) == length(vprev[[i]]))	
	# d <- length(v[[i]])
	# for (k in 1:d) {
		# if (length(vprev[[i]][[k]]) == 0) next
		# qrx <- qr(vprev[[i]][[k]])
		# if (qrx$rank == 0) next
		# v[[i]][[k]] <- lsfit(x = qr.Q(qrx)[,1:qrx$rank], 
			# y = v[[i]][[k]], intercept = FALSE)$residuals
	# }
# }	
# v
# }



#########################
# Deflation function for 
# 3D array of 2D data
#########################

# Inputs
# x:	3D array of data
# v:	list of two canonical vectors
# y:	vector of component scores (optional)  

# Calling (p,q,n) the dimensions of x (p,q = image dimensions, n = number of instances)
# the vectors in v should have respective length p and q, and the vector y (if provided)
# should have length n 

# deflate2D <- function(x, v, y = NULL)
# {
# n <- dim(x)[3]
# for (i in 1:2) {
	# s <- sqrt(sum(v[[i]]^2))
	# if (s > 0 && s != 1) v[[i]] <- v[[i]] / s
# }
# if (is.null(y)) {
	# y <- numeric(n)
	# for (t in 1:n) 
		# y[t] <- crossprod(v[[1]], x[,,t] %*% v[[2]])
# }
# vprod <- tcrossprod(v[[1]], v[[2]])
# for (t in 1:n) 
	# x[,,t] <- x[,,t] - y[t] * vprod
# return(x)
# }


#########################
# Deflation function for 
# 4D array of 3D data
#########################

# Inputs
# x:	data (array)
# v:	canonical vectors (list)
# y:	component scores (optional vector)  

# Writing (d1,...,dk,d(k+1)) for the dimensions of x 
# (first k dimensions = image dimensions, 
# (k+1)-th dimension = instance/individual)
# the list v must have length dk and contain 
# vectors of respective lengths d1,...,dk
# and the vector y (if provided) must have length d(k+1)

# deflate <- function(x, v, y = NULL)
# {
# ## Data dimensions
# dimx <- dim(x)
# k <- length(v)
# n <- dimx[k+1]
# ## Scale canonical vectors 
# s <- numeric(k)
# for (i in 1:k) {
	# s[i] <- sqrt(sum(v[[i]]^2))
	# if (s[i] > 0) v[[i]] <- v[[i]] / s[i]
# }
# s[s == 0] <- 1
# ## Kronecker product of canonical vectors
# kronv <- if (k == 2) {
	# kronecker(v[[2]], v[[1]])
# } else { kronecker(v[[3]], kronecker(v[[2]], v[[1]])) }
# ## Matricize data 
# dim(x) <- c(prod(dimx[1:k]), n)
# ## Calculate components if needed
# if (is.null(y)) {
	# y <- crossprod(kronv, x)
# } else {
	# y <- y / prod(s)
	# dim(y) <- c(1,n)
# }
# ## Return deflated data
# x <- x - kronv %*% y  
# dim(x) <- dimx
# x
# }

