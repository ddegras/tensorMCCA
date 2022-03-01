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
# # browser()
# ## Return deflated data
# x <- x - kronv %*% y  
# dim(x) <- dimx
# x
# }


################################################
# General function to deflate canonical vectors
# with respect to other canonical vectors
################################################


# Inputs
# v:	list of lists of single vectors 
# vprev: list of lists of matrices

deflate.v <- function(v, vprev, ortho = c("canon.t.1", "canon.t.all"))
{
stopifnot(is.list(v) && is.list(vprev))
stopifnot(length(v) == length(vprev))
ortho <- match.arg(ortho)
m <- length(v)
d <- if (ortho == "canon.t.1") {
	rep(1L, m) } else if (ortho == "canon.t.all") {
	sapply(v, length) }
for (i in 1:m) {
for (k in 1:d[i]) {
	qrx <- qr(vprev[[i]][[k]])
	if (qrx$rank == 0) next
	v[[i]][[k]] <- lsfit(x = qr.Q(qrx)[,1:qrx$rank], 
		y = v[[i]][[k]], intercept = FALSE)$residuals
	}
}
v
}



###########################################
# General function to deflate data tensors
###########################################

deflate.x <- function(x, v = NULL, y = NULL, ortho = c("block.score", 
	"global.score", "canon.t.1", "canon.t.all"), check.args = TRUE)
{
if (check.args) test <- check.arguments(x, v)
m <- length(x)
dimx <- lapply(x, dim)
n <- tail(dimx[[1]], 1)
ortho <- match.arg(ortho)

if (ortho == "block.score") {
	if(is.null(v) && is.null(y))
		stop(paste("At least one of 'v' or 'y' must be specified",
			"if 'ortho' is set to 'block.score'"))
	if (is.null(y)) y <- image.scores(x, v)
	nrmy <- sqrt(colSums(y^2))
	nrmy[nrmy == 0] <- 1
	for (i in 1:m) {
		y[,i] <- y[,i] / nrmy[i]
		dim(x[[i]]) <- c(length(x[[i]]) / n, n)
		x[[i]] <- x[[i]] - tcrossprod(x[[i]] %*% y[,i], y[,i])
		dim(x[[i]]) <- dimx[[i]]
	}
}

if (ortho == "global.score") {
	if(is.null(v) && is.null(y))
		stop(paste("At least one of 'v' or 'y' must be specified",
			"if 'ortho' is set to 'global.score'"))
	if (is.null(y)) y <- image.scores(x, v)
	y <- rowSums(y) # global score (unscaled)
	y <- y / sqrt(sum(y^2))
	for (i in 1:m) {
		dim(x[[i]]) <- c(length(x[[i]]) / n, n)
		x[[i]] <- x[[i]] - tcrossprod(x[[i]] %*% y, y)
		dim(x[[i]]) <- dimx[[i]]
	}	
}

if (ortho == "canon.t.1") {
	if (is.null(v)) 
		stop(paste("Argument 'v' must be specified",
			"if 'ortho' is set to 'canon.t.1'"))
	for (i in 1:m) {		
		vi1 <- v[[i]][[1]] # mode 1 canonical vector for i-th dataset
		vi1 <- vi1  / sqrt(sum(vi1^2)) # rescale to unit norm 
		## Matricize x(i) according to mode 1
		dim(x[[i]]) <- c(dimx[[i]][1], prod(dimx[[i]][-1])) 
		## Orthogonal projection
		x[[i]] <- x[[i]] - vi1 %*% crossprod(vi1, x[[i]]) 
		## Reshape x(i) to tensor
		dim(x[[i]]) <- dimx[[i]] 
	}
}

if (ortho == "canon.t.all") {
	if (is.null(v)) 
		stop(paste("Argument 'v' must be specified", 
			"if 'ortho' is set to 'canon.t.all'"))
	for (i in 1:m) {
		ndimxi <- length(dimx[[i]])
		nvi <- length(v[[i]]) 
		for (k in 1:nvi) {
			vik <- v[[i]][[k]]
			vik <- vik  / sqrt(sum(vik^2))
			if (k == 1) {
				# For orthogonal projection in mode 1, 
				# no need to permute array dimensions
				## Matricize x(i) along mode k
				dim(x[[i]]) <- c(dimx[[i]][1], prod(dimx[[i]][-1]))
				x[[i]] <- x[[i]] - vik %*% crossprod(vik, x[[i]]) 
				dim(x[[i]]) <- dimx[[i]]	 
			} else {
				## Permute modes 1 and k in x(i) to enable multiplication
				perm <- 1:(nvi+1)
				perm[c(1,k)] <- c(k,1) 
				x[[i]] <- aperm(x[[i]], perm)
				## Matricize x(i) along mode k
				dim(x[[i]]) <- c(dimx[[i]][k], prod(dimx[[i]][-k]))
				x[[i]] <- x[[i]] - vik %*% crossprod(vik, x[[i]]) 
				dim(x[[i]]) <- dimx[[i]]	
				## Permute dimensions back to original
				x[[i]] <- aperm(x[[i]], perm) 
			}				
		}	
	}
}

return(x)
}

