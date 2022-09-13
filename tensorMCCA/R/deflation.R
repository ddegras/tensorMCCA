

###########################################
# General function to deflate data tensors
###########################################


deflate.x <- function(x, v = NULL, score = NULL, # ortho.mode = NULL, 
	scope = c("block", "global"), check.args = TRUE)
{
scope <- match.arg(scope)
if (!is.null(v) && scope == "global")
	stop(paste("Deflation with respect to a tensor 'v' is only",
	"defined for 'scope' equal to 'block'."))
if (check.args) {
	test <- check.arguments(x, v) 
	m <- length(x)
	dimx <- lapply(x, dim)
	d <- sapply(dimx, length) - 1L
	n <- dimx[[1]][d[1] + 1L]
	if (is.null(v)  == is.null(score)) 
		stop(paste("Please specify exactly one of the arguments",
		"'v' and 'score' and leave the other set to NULL."))
	# if (!is.null(v)) {		
		# if (is.null(ortho.mode))
			# stop("If 'v' is specified, then so must 'ortho.mode'.")
		# r <- NCOL(v)
		# ortho.mode <- drop(ortho.mode)
		# if (!((r == 1L && length(ortho.mode) == m) || 
			# (r > 1L) && identical(dim(ortho.mode), c(m, r))))
			# stop(paste("'ortho.mode' must be a matrix of lists",
			# "with dimensions matching those of 'v'"))
		# for (i in 1:m) {
			# modes <- unlist(
				# if (r == 1L) { ortho.mode[[i]] 
				# } else ortho.mode[i,])
			# if (!all(modes %in% (1:d[i])))
				# stop(paste("Values in 'ortho.mode[[i]]' or",
				# "'ortho.mode[i,]' not compatible with the",
				# "dimensions of 'x'."))
		# }	
	# }
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
p <- mapply(head, dimx, d, SIMPLIFY = FALSE) 
n <- tail(dimx[[1]], 1)
eps <- 1e-14

## Deflate with respect to canonical weight tensors
if (!is.null(v)) {
	if (!is.matrix(v)) v <- as.matrix(v)
	r <- ncol(v)
	for (i in 1:m) {
		vi <- sapply(v[i,], function(vi) Reduce(kronecker, rev(vi)))
		qrvi <- qr(vi)
		if (qrvi$rank == 0) next
		vi <- qr.Q(qrvi)[,1:qrvi$rank]
		dim(x[[i]]) <- c(prod(p[[i]]), n)
		x[[i]] <- x[[i]] - vi %*% crossprod(vi, x[[i]])
		dim(x[[i]]) <- dimx[[i]]
	}
	return(x)
	# cnstr <- set.ortho.mat(v = v, modes = ortho.mode)
	# x <- mapply(tnsr.mat.prod, x = x, mat = cnstr$mat,
		# modes = cnstr$modes, SIMPLIFY = FALSE)
	# return(x)
}


## Deflate with respect to canonical scores
if (!is.null(score)) {
	score <- drop(score)
	if (scope == "block") {
		arr <- (length(dim(score)) > 2L)  
		if (arr) score <- aperm(score, c(1,3,2))
		dims <- 2:length(dim(score))
		score <- sweep(score, dims, colMeans(score))
		if (!arr) {
			nrm <- sqrt(colSums(score^2))
			zero <- (nrm < eps)
			if (all(zero)) return(x)
			score <- sweep(score, dims, nrm, "/")
		}
		for (i in 1:m) {
			if (arr) {
				qrscorei <- qr(score[,,i])
				if (qrscorei$rank == 0L) next
				scorei <- qr.Q(qrscorei)[,1:qrscorei$rank] 
			} else if (zero[i]) {
				next
			} else {
				scorei <- score[,i]
			}
			dim(x[[i]]) <- c(prod(p[[i]]), n)
			x[[i]] <- x[[i]] - 
				tcrossprod(x[[i]] %*% scorei, scorei)
			dim(x[[i]]) <- dimx[[i]]
		}
	}
	if (scope == "global") {
		if (is.vector(score)) {
			score <- score - mean(score)
			nrm <- sqrt(sum(score^2))
			score <- if (nrm < eps) {
				numeric(n) } else score / nrm 
		} else {
			score <- sweep(score, 2, colMeans(score), "-")
			qrscore <- qr(score)
			score <- qr.Q(qrscore)[, 1:qrscore$rank]
		}		
		for (i in 1:m) {
			dim(x[[i]]) <- c(prod(p[[i]]), n)
			x[[i]] <- x[[i]] - 
				tcrossprod(x[[i]] %*% score, score)
			dim(x[[i]]) <- dimx[[i]]
		}	
	} 
	return(x)
}

NULL

}
