

###########################################
# General function to deflate data tensors
###########################################

# This function is used inside 'mcca.cov', 'mcca.cor', and 'permutation.test' 
# It is not meant to be directly called by the user

deflate.x <- function(x, ortho, scope = c("block", "global"))
{
scope <- match.arg(scope)
stopifnot(length(x) == NROW(ortho))
if (!is.matrix(ortho)) 
	dim(ortho) <- c(length(ortho), 1L)
toexpand <- which(sapply(ortho, is.list))
for (idx in toexpand)
	ortho[[idx]] <- outer.prod.nodim(ortho[[idx]])
tovectorize <- which(sapply(ortho, is.array))
for (idx in tovectorize)
	ortho[[idx]] <- as.vector(ortho[[idx]])
m <- length(x)
r <- ncol(ortho)
dimx <- lapply(x, dimfun)
n <- tail(dimx[[1]],1)
pp <- sapply(dimx, function(dims) prod(dims[-length(dims)]))

if (scope == "block") {	
	for (i in 1:m) {
		orthoi <- do.call(cbind, ortho[i,])
		qi <- qr(orthoi)
		if (qi$rank == 0) next
		qi <- qr.Q(qi)[, 1:qi$rank, drop = FALSE]
		dim(x[[i]]) <- c(pp[i],n)
		x[[i]] <- x[[i]] - qi %*% crossprod(qi, x[[i]])
		dim(x[[i]]) <- dimx[[i]]
	}
}

if (scope == "global") {
	orthoall <- matrix(, sum(pp), r)
	csumpp <- c(0, cumsum(pp))
	for (i in 1:m) {
		idx <- (csumpp[i]+1):csumpp[i+1]
		orthoall[idx,] <- do.call(cbind, ortho[i,])
	}
	qall <- qr(orthoall)
	if (qall$rank == 0) return(x)
	qall <- qr.Q(qall)[, 1:qall$rank, drop = FALSE]
	for (i in 1:m) {
		idx <- (csumpp[i]+1):csumpp[i+1]
		qi <- qall[idx,,drop=FALSE]
		dim(x[[i]]) <- c(pp[i], n)
		x[[i]] <- x[[i]] - qi %*% crossprod(qi, x[[i]])		
		dim(x[[i]]) <- dimx[[i]]
	}
}

x
}



# deflate.x <- function(x, v = NULL, score = NULL, 
	# scope = c("block", "global"), check.args = TRUE)
# {
# scope <- match.arg(scope)
# if (!is.null(v) && scope == "global")
	# stop(paste("Deflation with respect to a tensor 'v' is only",
	# "defined for 'scope' equal to 'block'."))
# if (is.numeric(x)) x <- list(x)
# if (check.args) test <- check.arguments(x, v) 
# if (is.null(v) && is.null(score)) return(x)
# if (!is.null(v) && !is.null(score)) 
	# stop("Please specify exactly one of the arguments ",
		# "'v' and 'score' and leave the other equal to NULL.")

# ## Data dimensions
# m <- length(x)
# dimfun <- function(x) if (is.vector(x)) c(1, length(x)) else dim(x)
# dimx <- lapply(x, dimfun)
# d <- sapply(dimx, length) - 1L
# n <- tail(dimx[[1]], 1) 
# p <- mapply(head, dimx, d, SIMPLIFY = FALSE) 
# eps <- 1e-14

# ## Deflate with respect to canonical weights if requested
# if (!is.null(v)) {
	# r <- NCOL(v)
	# if (!is.matrix(v)) dim(v) <- c(m,r)
	# for (i in 1:m) {
		# vi <- sapply(v[i,], outer.prod.nodim)
		# # vi <- sapply(v[i,], function(vi) Reduce(kronecker, rev(vi)))
		# qrvi <- qr(vi)
		# if (qrvi$rank == 0) next
		# vi <- qr.Q(qrvi)[,1:qrvi$rank]
		# dim(x[[i]]) <- c(prod(p[[i]]), n)
		# x[[i]] <- x[[i]] - vi %*% crossprod(vi, x[[i]])
		# dim(x[[i]]) <- dimx[[i]]
	# }
	# return(x)
# }

# ## Deflate with respect to canonical scores if requested
# stopifnot(is.numeric(score))
# if (scope == "block" && length(dim(score)) < 3)
	# dim(score) <- c(dim(score),1,1,1)[1:3]
# if (scope == "global" && length(dim(score)) < 2)
	# dim(score) <- c(dim(score),1,1)[1:2]
# dimscore <- dim(score)
# test1 <- (scope == "block" && length(dimscore) == 3 && 
	# dimscore[1] == n && dimscore[2] == m)
# test2 <- (scope == "global" && is.matrix(score) && 
	# dimscore[1] == n)
# if (!(test1 || test2))
	# stop("If specified, 'score' should be either", 
		# "a matrix with n rows (if 'scope' = 'global')\n", 
		# "or a 3D array with first two dimensions equal to (n,m)",
		# "(if 'scope' = 'block')\nwhere m is the number of data arrays",
		# "(length of 'x') and n is the last dimension",
		# "of the data arrays.")
# if (scope == "block") {
	# score <- aperm(score, c(1,3,2)) # dim = c(n,r,m)
	# ## Ensure scores are centered
	# score <- sweep(score, 2:3, colMeans(score))
	# for (i in 1:m) {
		# ## Ensure scores are orthogonal
		# qrscorei <- qr(score[,,i])
		# if (qrscorei$rank == 0) next
		# scorei <- qr.Q(qrscorei)[,1:qrscorei$rank] 
		# dim(x[[i]]) <- c(prod(p[[i]]), n)
		# x[[i]] <- x[[i]] - 
			# tcrossprod(x[[i]] %*% scorei, scorei)
		# dim(x[[i]]) <- dimx[[i]]
	# }
# }
# if (scope == "global") {
	# score <- sweep(score, 2, colMeans(score))
	# qrscore <- qr(score)
	# if (qrscore$rank == 0) return(x)			
	# score <- qr.Q(qrscore)[, 1:qrscore$rank]		
	# for (i in 1:m) {
		# dim(x[[i]]) <- c(prod(p[[i]]), n)
		# x[[i]] <- x[[i]] - 
			# tcrossprod(x[[i]] %*% score, score)
		# dim(x[[i]]) <- dimx[[i]]
	# }	
# }

# x

# }
