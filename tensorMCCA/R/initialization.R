##############################
# Initialization of canonical
# vectors with ones
##############################


init.v.ones <- function(x, type = c("norm", "var"), 
	cnstr = c("block", "global"), balance = TRUE)
{
## Check argument x 
test <- check.arguments(x)

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
dim.img <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
ndim.img <- sapply(dim.img, length)
n <- tail(dimx[[1]], 1)

## Scaling constraints
type <- match.arg(type)   # norm or variance constraints
cnstr <- match.arg(cnstr) # block or global constraints

## Set canonical vectors to 1's
v <- vector("list",m)
for (i in 1:m)
	v[[i]] <- lapply(dim.img[[i]], function(len) rep_len(1,len))

## Scale canonical vectors
v <- scale.v(v, x, type, cnstr, balance)

return(v)
}



########################
# Random initialization 
# of canonical vectors 
########################


init.v.random <- function(x, type = c("norm", "var"), 
	cnstr = c("block", "global"), balance = TRUE)
{
## Check argument x if required
test <- check.arguments(x)

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
dim.img <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
ndim.img <- sapply(dim.img, length)

## Scaling constraints
type <- match.arg(type)   # norm or variance constraints
cnstr <- match.arg(cnstr) # block or global constraints

## Initialize canonical vectors randomly
v <- vector("list",m)
for (i in 1:m)
	v[[i]] <- lapply(dim.img[[i]], function(len) runif(len))

## Scale canonical vectors
v <- scale.v(v, x, type, cnstr, balance)
return(v)
}



###########################
# SVD-based initialization
# of canonical vectors
###########################


# init.v.svd <- function(x, type = c("norm", "var"), 
	# cnstr = c("block", "global"), balance = TRUE)
# {
# ## Check argument x if required
# test <- check.arguments(x)

# ## Data dimensions
# m <- length(x) # number of datasets 
# dimx <- lapply(x, dim) # full data dimensions
# dim.img <- lapply(dimx, function(idx) idx[-length(idx)]) 
# # image dimensions (ignore last dimension of datasets = replications)
# ndim.img <- sapply(dim.img, length)
# n <- tail(dimx[[1]], 1)
	
# ## Scaling constraints
# type <- match.arg(type)   # norm or variance constraints
# cnstr <- match.arg(cnstr) # block or global constraints

# ## Initialize canonical vectors to 1's
# v <- vector("list",m)
# for (i in 1:m)
	# v[[i]] <- lapply(dim.img[[i]], function(len) rep_len(1,len))
	
# ## Update canonical vector in mode 1 by SVD 	
# dim.k <- sapply(dim.img, "[", 1)
# csum <- cumsum(c(0,dim.k))
# compk <- matrix(0, csum[m+1], n)			
# idx.rows <- lapply(1:m, 
	# function(i) seq.int(csum[i] + 1, csum[i+1]))
# for (i in 1:m) {
	# compk[idx.rows[[i]],] <- if (ndim.img[i] == 1) {
		# x[[i]] } else {
		# apply(x[[i]], c(1, ndim.img[i] + 1), mean) }
# }
# vk <- if (min(dim(compk)) >= 3) {
	# svds(compk, 1, 1, 0)$u } else { svd(compk, 1, 0)$u }
# for (i in 1:m) 
	# v[[i]][[1]] <- vk[idx.rows[[i]]]	
# v <- scale.v(v, x, type, cnstr, balance)

# ## Update canonical vector in mode 2 by SVD 	
# idx.i <- which(ndim.img >= 2)
# if (length(idx.i) > 0) {
	# dim.k <- sapply(dim.img[idx.i], "[", 2)
	# csum <- cumsum(c(0,dim.k))
	# compk <- matrix(0, tail(csum,1), n)			
	# idx.rows <- lapply(seq_along(idx.i), 
		# function(i) seq.int(csum[i] + 1, csum[i+1]))
	# for (i in seq_along(idx.i)) {
		# ii <- idx.i[i]
		# iprod <- x[[ii]]
		# dims <- dimx[[ii]]
		# if (ndim.img[ii] == 3) 
			# iprod <- apply(iprod, c(1,2,4), mean)
		# dim(iprod) <- c(dims[1], dims[2] * n)
		# iprod <- crossprod(v[[ii]][[1]], iprod)
		# dim(iprod) <- c(dims[2], n)
		# compk[idx.rows[[i]],] <- iprod
	# }
	# vk <- if (min(dim(compk)) >= 3) {
		# svds(compk, 1, 1, 0)$u } else { svd(compk, 1, 0)$u }
	# for (i in seq_along(idx.i)) 
		# v[[idx.i[i]]][[2]] <- vk[idx.rows[[i]]]
	# v <- scale.v(v, x, type, cnstr, balance)
# }
	
# ## Update canonical vector in mode 3 by SVD 	
# idx.i <- which(ndim.img == 3)
# if (length(idx.i) > 0) {
	# dim.k <- sapply(dim.img[idx.i], "[", 3)
	# csum <- cumsum(c(0,dim.k))
	# compk <- matrix(0, tail(csum,1), n)			
	# idx.rows <- lapply(seq_along(idx.i), 
		# function(i) seq.int(csum[i] + 1, csum[i+1]))
	# for (i in seq_along(idx.i)) {
		# ii <- idx.i[i]
		# iprod <- x[[ii]]
		# dims <- dimx[[ii]]	
		# dim(iprod) <- c(dims[1], prod(dims[-1]))
		# iprod <- crossprod(v[[ii]][[1]], iprod)
		# dim(iprod) <- c(dims[2], dims[3] * n)
		# iprod <- crossprod(v[[ii]][[2]], iprod)
		# dim(iprod) <- c(dims[3], n)
		# compk[idx.rows[[i]],] <- iprod
	# }
	# vk <- if (min(dim(compk)) >= 3) {
		# svds(compk, 1, 1, 0)$u } else { svd(compk, 1, 0)$u }
	# for (i in seq_along(idx.i)) 
		# v[[idx.i[i]]][[3]] <- vk[idx.rows[[i]]]
	# v <- scale.v(v, x, type, cnstr, balance)
# }		

# return(v)	
# }



###########################
# SVD-based initialization
# of canonical vectors 
# (Quick Rank 1)
###########################


init.mcca.svd <- function(x, r = 1, type = c("cov", "cor"), 
	cnstr = c("block", "global"), center = TRUE, balance = TRUE)
{
## Check argument x if required
test <- check.arguments(x)
tol <- 1e-14

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
d <- sapply(p, length)
n <- tail(dimx[[1]], 1)
r <- min(as.integer(r), n)

## Scaling constraints
type <- match.arg(type)   # norm or variance constraints
type <- switch(type, cov = "norm", cor = "var")
cnstr <- match.arg(cnstr) # block or global constraints
# ortho <- match.arg(ortho)

## Initialize canonical vectors to 1's
v <- vector("list", m)

## Unfold data along last mode (individuals/objects) and concatenate
for (i in 1:m) 
	dim(x[[i]]) <- c(prod(p[[i]]), n)
x <- do.call(rbind, x)

## Data centering
if (center) 
	x <- x - rowMeans(x)

## Initial SVD
svdfun <- if (max(dim(x)) > 2) {
	function(x) svds(x, k = r) } else {
	function(x) svd(x, nu = r, nv = r) } 
x <- svdfun(x)$u
if (r == 1) dim(x) <- c(length(x), 1)

## Iteratively unfold singular vectors and recalculate SVD
lenx <- sapply(p, prod)
start <- c(0, cumsum(lenx[-m])) + 1
end <- cumsum(lenx)
for (i in 1:m) {
	pi <- p[[i]]
	for (l in 1:r) {
		xi <- x[start[i]:end[i], l]
		for (k in 1:d[i]) {
			if (l == 1) 
				v[[i]][[k]] <- matrix(nrow = pi[k], ncol = r)
			if (k < d[i]) {
				dim(xi) <- c(pi[k], length(xi) / pi[k])
				svdx <- svdfun(xi)
				v[[i]][[k]][, l] <- svdx$u[,1]
				xi <- svdx$v[,1]
			} else {
				v[[i]][[k]][, l] <- xi 
			}
		}
	}
}

## Scale initial canonical vectors as required	
v <- scale.v(v, x, type, cnstr, balance)

return(v)	
}


############################################
# Initialize MCCA by solving pairwise CCA's
############################################



init.mcca.cca <- function(x, r = 1, k = NULL, c = 1, type = c("cov", "cor"), 
	cnstr = c("block", "global"), center = TRUE, balance = TRUE)
{
## Check argument x if required
test <- check.arguments(x)
tol <- 1e-14

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
d <- sapply(p, length)
n <- tail(dimx[[1]], 1)
r <- min(r, n)
if (is.null(k)) k <- n else k <- min(k, n) 
k <- max(k, r)	
	
## Scaling constraints
objective.type <- match.arg(type)   # norm or variance constraints
type <- switch(objective.type, cov = "norm", cor = "var")
cnstr <- match.arg(cnstr) # block or global constraints

## Objective weights
stopifnot(length(c) == 1 || (is.matrix(c) && all(dim(c) == length(x))))
stopifnot(all(c >= 0) && any(c > 0))
if (length(c) == 1) {
	c <- matrix(1/m^2, m, m)
} else {
	c <- (c + t(c)) / (2 * sum(c))
}

## Data centering
if (center) {
	for(i in 1:m) {
	    mu <- rowMeans(x[[i]], dims = ndim.img[i])
	    x[[i]] <- x[[i]] - as.vector(mu)
	}
}

## Calculate SVDs of each unfolded dataset
svdx <- vector("list", m) 
rankx <- integer(m)
for (i in 1:m) {
	svdx[[i]] <- if(k < n) {
		svds(matrix(x[[i]], ncol = n), k = k)
	} else {
		svd(matrix(x[[i]], ncol = n), nv = k)
	}
	pos <- (svdx[[i]][["d"]] > tol)
	rankx[i] <- sum(pos)
	if (all(!pos)) {		
		svdx[[i]][["d"]] <- 0
		svdx[[i]][["u"]] <- 0
		svdx[[i]][["v"]] <- 0
	} else if (!all(pos)) {
		svdx[[i]][["d"]] <- svdx[[i]][["d"]][pos]
		svdx[[i]][["u"]] <- svdx[[i]][["u"]][,pos]
		svdx[[i]][["v"]] <- svdx[[i]][["v"]][,pos]
	}
}

## Reduced data
xred <- if (objective.type == "cov") {
	lapply(svdx, function(x) t(x[["v"]]) * x[["d"]])
} else {
	lapply(svdx, function(x) t(x[["v"]]))	
}

## Find linear combinations of canonical vectors for each dataset
## that maximize objective function
# v <- vector("list", m)
if (objective.type == "cov" && cnstr == "global") {
	xred <- do.call(rbind, xred)
	vred <- if (max(dim(xred)) > 2) {
		svds(xred, k = r)$u
	} else { svd(xred, nu = r)$u }
	if (r == 1) vred <- as.matrix(vred)
	end <- cumsum(rankx)
	vred <- lapply(1:m, function(i) sqrt(m) * 
		vred[seq.int(to = end[i], length.out = rankx[i]),])
	v <- lapply(1:m, function(i) svdx[[i]]$u %*% vred[[i]])

} else if (objective.type == "cov" && cnstr == "block") {
	vred <- vector("list", m) # reduced canonical vectors
	objective0 <- array(0, dim = c(m, m, r))
	for (i in 1:m) {
		vred[[i]] <- array(0, dim = c(rankx[i], r, m))
		for (j in 1:i) {
			mat <- tcrossprod(xred[[i]], xred[[j]])
			svdij <- svd(mat, nu = r, nv = r)
			reff <- length(svdij$d)
			objective0[i,j,1:reff] <- c[i,j] * svdij$d				
			vred[[i]][, 1:reff, j] <-  svdij$u
			if (j < i) vred[[j]][, 1:reff, i] <- svdij$v
		}
	}
	
} else {
	vred <- vector("list", m) # reduced canonical vectors
	objective0 <- array(0, dim = c(m, m, r))
	for (i in 1:m) {
		vred[[i]] <- array(0, dim = c(rankx[i], r, m))
		for (j in 1:i) {
			mat <- tcrossprod(xred[[i]], xred[[j]])
			svdij <- svd(mat, nu = r, nv = r)
			reff <- min(rankx[i], rankx[j], r)
			objective0[i, j, 1:reff] <- c[i,j] * svdij$d			
			vred[[i]][, 1:reff, j] <- svdij$u / svdx[[i]]$d[1:reff]
			if (j < i) 
				vred[[j]][, 1:reff, i] <- svdij$v / svdx[[j]]$d[1:reff] 
		}
	}	
}

#@@@@@@ TO DO Improvement: rank 1 tensor approximation first, LSAP second
	
## Find best overall assignment in pairwise CCA		
## (Linear Sum Assignment Problem)
if (objective.type == "cor" || cnstr == "block") {
 	v <- vector("list", m)
	idxj <- apply(objective0, 3, solve_LSAP, maximum = TRUE)	
	for (i in 1:m) {
		vtmp <- matrix(, rankx[i], r)
		for (k in 1:r) vtmp[,k] <- vred[[i]][, k, idxj[i,k]]
		v[[i]] <- svdx[[i]]$u %*% vtmp
	}
}

## Rank 1 tensor approximation to long canonical vectors
for (i in 1:m) {
	pi <- p[[i]]
	vi <- v[[i]]
	v[[i]] <- lapply(pi, function(p) matrix(nrow = p, ncol = r))
	for (l in 1:r) {
		vil <- vi[, l]
		for (k in 1:d[i]) {
			if (k < d[i]) {
				dim(vil) <- c(pi[k], length(vi) / pi[k])
				svdi <- svd(vil)
				v[[i]][[k]][, l] <- svdx$u
				vil <- svdx$v
			} else {
				v[[i]][[k]][, l] <- vil 
			}
		}
	}
}

## Scale the canonical tensors as needed
## (norm constraints are already enforced)
if (type == "var")
	v <- scale.v(v, x, type = type, balance = balance)

v

}


