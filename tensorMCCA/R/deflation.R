

###########################################
# General function to deflate data tensors
###########################################


deflate.x <- function(x, v = NULL, score = NULL, ortho.mode = NULL, 
	scope = c("block", "global"), check.args = TRUE)
{
scope <- match.arg(scope)
if (!is.null(v) && scope == "global")
	stop(paste("Deflation with respect to a tensor 'v' is only",
	"defined for 'scope' equal to 'block'."))
if (check.args) {
	test <- check.arguments(x, v) 
	m <- length(x)
	n <- tail(dim(x[[1]]), 1)
	if (is.null(v)  == is.null(score)) 
		stop(paste("Please specify exactly one of the arguments",
		"'v' and 'score' and leave the other set to NULL."))
	if (!is.null(v)) {		
		if (is.null(ortho.mode))
			stop("If 'v' is specified, then so must 'ortho.mode'.")
		r <- NCOL(v)
		ortho.mode <- drop(ortho.mode)
		if (!((r == 1L && length(ortho.mode) == m) || 
			(r > 1L) && identical(dim(ortho.mode), c(m, r))))
			stop(paste("'ortho.mode' must be a matrix of lists",
			"with dimensions matching those of 'v'"))			
	}
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

## Deflate with respect to canonical tensors
if (!is.null(v)) {
	r <- NCOL(v)
	if (r == 1L) {
		for (i in 1:m) {
			for (k in ortho.mode[[i]]) {
				vik <- v[[i]][[k]]
				nrm <- sqrt(sum(vik^2))
				if (nrm < eps) next
				vik <- vik / nrm
				perm <- 1:(d[i] + 1L)
				if (k > 1L) {
					perm[c(1,k)] <- c(k,1)
					x[[i]] <- aperm(x[[i]], perm)
				}
				dim(x[[i]]) <- c(dimx[[i]][k], 
					prod(dimx[[i]][-k]))
				x[[i]] <- x[[i]] - vik %*% crossprod(vik, x[[i]])
				dim(x[[i]]) <- dimx[[i]][perm]
				if (k > 1L) x[[i]] <- aperm(x[[i]], perm) 	
			}		
		}
	} else {
		for (i in 1:m) {
			for (k in 1:d[i]) {
				idxl <- which(sapply(ortho.mode[,i], 
					is.element, el = k))
				if (length(idxl) == 0) next
				vik <- sapply(v[i,idxl], "[[", k)
				qrvik <- qr(vik)
				if (qrvik$rank == 0L) next
				vik <- qr.Q(qrvik)[, 1:qrvik$rank]
				perm <- 1:(d[i] + 1L)
				if (k > 1L) {
					perm[c(1,k)] <- c(k,1)
					x[[i]] <- aperm(x[[i]], perm)
				}
				dim(x[[i]]) <- c(dimx[[i]][k], 
					prod(dimx[[i]][-k]))
				x[[i]] <- x[[i]] - (vik %*% crossprod(vik, x[[i]]))
				dim(x[[i]]) <- dimx[[i]][perm]
				if (k > 1L) x[[i]] <- aperm(x[[i]], perm) 		
			}
		}		
	}
	return(x)
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




# if (!is.null(v) && scope == "global") {
	# nmodes <- sapply(ortho.mode, length)
	# nvecs <- prod(nmodes)
	# fn <- function(v, i) 
	# for (ii in 1:nvecs) {
		# k <- arrayInd(ii, nmodes)	
		# vv <- unlist(mapply("[[", v, k, SIMPLIFY = FALSE))
		# for (i in 1:m) {
			
			
		# }
		# for (k in 1:d[i]) {
			# if (is.null(v[[i]][[k]])) next
			# ## Basis of orthogonal space
			# Q <- qr.Q(qr(v[[i]][[k]])) 
			# ## Permute data tensor (exchange modes 1 and k)
			# if (k > 1) {
				# perm <- 1:(d[i]+1)
				# perm[c(1,k)] <- perm[c(k,1)]
				# x[[i]] <- aperm(x[[i]], perm)
			# }
			# ## Unfold permuted tensor along mode 1 
			# ## (original along mode k)
			# dim(x[[i]]) <- c(dimx[[i]][k], prod(dimx[[i]][-k])) 
			# ## Orthogonal projection
			# x[[i]] <- x[[i]] - Q %*% crossprod(Q, x[[i]]) 
			# ## Reshape x(i) to tensor
			# dim(x[[i]]) <- if (k == 1) dimx[[i]] else dimx[[i]][perm]
			# ## Permute modes back
			# if (k > 1) 
				# x[[i]] <- aperm(x[[i]], perm)
		# }
	# }
# }






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



