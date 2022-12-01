#########################################
# Function to simulate canonical weights 
#########################################

simulate.v <- function(p, r, ortho = TRUE) {
if (any(unlist(p) < 1)) 
	stop("The values in 'dims' must be positive integers.")
stopifnot(r >= 0)
if (!is.list(p)) p <- list(p)
p <- lapply(p, as.integer)
r <- as.integer(r)
if (r == 1) ortho <- FALSE
if (ortho && any(unlist(p) < r))
	warning("The generated canonical weight tensors",
		" are not orthogonal to each other in all modes.")
d <- sapply(p, length)
m <- length(p)
v <- vector("list", m * max(1, r))
dim(v) <- c(m, max(1, r))
for (i in 1:m) {
	if (r == 0) {
		v[[i]] <- lapply(p[[i]], numeric)
		next
	} 
	for (k in 1:d[i]) {
		vik <- matrix(runif(p[[i]][k] * r,-1, 1), p[[i]][k], r)
		vik <- if (ortho) { 
			svd(vik)$u
		} else {
			sweep(vik, 2, sqrt(colSums(vik^2)), "/")
		}
		nik <- ncol(vik)
		for (l in 1:r) 
			v[[i,l]][[k]] <- vik[, min(l, nik)]
	}
}
v
}


################################
# Function to simulate all data
################################



simulate.factor.model <- function(dimx, r, sigma2, xi, 
	scale = c("block", "global"), ortho = FALSE)
{
## Data dimensions
d <- sapply(dimx, length) - 1L
m <- length(dimx)
n <- sapply(dimx, tail, 1)
if (any(n != n[1])) 
	stop(paste("The last values of the components of",
	"'dimx' must all be equal (= number of instances)."))
n <- n[1]
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
pp <- sapply(p, prod)
stopifnot(r >= 0)
r0 <- r <- as.integer(r)
sigma2 <- rep_len(sigma2, r)
stopifnot(length(xi) == 1 || length(xi) == length(dimx))
nxi <- length(xi)
scale <- match.arg(scale)

## Outputs
x <- vector("list", m)
v <- if (scale == "block") {
	simulate.v(p, r, ortho)
} else if (ortho == FALSE) {
	scale.v(simulate.v(p, r, FALSE), "norm", "global", 
		check.args = FALSE)
} else {
	orthogonalize.global(simulate.v(p, r, FALSE))
}
if (r == 0) {
	score <- numeric(n)
} else {
	score <- sapply(sqrt(sigma2), function(sig) rnorm(n, 0, sig))
	score <- sweep(score, 2, colMeans(score))
	score[is.nan(score)] <- 0 
}

for (i in 1:m) {	
	## Simulate signal 
	signal <- 0 
	if (r > 0) {
		vi <- matrix(, pp[i], r)
		for (l in 1:r) 
			vi[,l] <- Reduce(kronecker, rev(v[[i,l]]))
		signal <- tcrossprod(vi, score)
		dim(signal) <- c(p[[i]], n)
	}	
	## Simulate noise
	if (is.numeric(xi)) {
		noise <- rnorm(pp[i] * n, 
			sd = sqrt(xi[min(i, nxi)]))	
	} else {
		if ((nxi == 1 && i == 1) || nxi > 1)
			R <- chol(xi[[i]])
		noise <- crossprod(R, matrix(rnorm(pp[i] * n), pp[i], n))
	}
	dim(noise) <- c(p[[i]], n)		
	noise <- noise - as.vector(rowMeans(noise, dims = d[i]))
	x[[i]] <- signal + noise
}

list(x = x, v = v, score = score)
}


########################################
# Simulate model (2) with block scores
# under independence assumption on block scores
########################################


# X_it = sum_(l=1:r) v_i^l s_it^l + e_it       (2)


# ASSUMPTIONS ON THE SCORES
# with s_it^l, l=1:r, independent for each (i,t)
# but correlation among the s_it, i=1:m (for each (l,t))


# Inputs
# dimx:		dimensions of data tensors (list of length m)
# r:		rank of data decomposition
# sigma:	array of score covariances (m x m x r) 



# simulate.x2 <- function(dimx, r, sigma)
# {
# ## Data dimensions
# m <- length(dimx)
# n <- sapply(dimx, tail, 1)
# if (any(n != n[1])) 
	# stop(paste("Last values in each component of list",
	# "'dimx' must all be equal (number of replications)."))
# n <- n[1]
# p <- sapply(dimx, function(idx) idx[-length(idx)])
# d <- sapply(dimx, length) - 1

# ## Check dimensions of sigma and replicate if needed
# test <- is.matrix(sigma) && dim(sigma) == c(m,m)
# stopifnot(test || all(dim(sigma) == c(m,m,r)))
# if (test) sigma <- replicate(r, sigma)

# ## Outputs
# v <- vector("list", m * r)
# dim(v) <- c(m, r)
# x <- lapply(dimx, function(dims) array(0, dims))
# block.scores <- array(dim = c(n,m,r))

# ## Generate canonical vectors
# rand.vec <- function(pp) {
	# v <- runif(pp)
	# v / sqrt(sum(v^2))
# }
# for (i in 1:m)
# for (l in 1:r)
# v[[i,l]] <- lapply(p[[i]], rand.vec)

# ## Generate data
# for (l in 1:r) {    
	# ## Generate block scores
	# R <- chol(sigma[,,l])
	# tmp <- matrix(rnorm(m*n), m, n)
	# block.scores[,,l] <- crossprod(tmp, R)
            	
	# ## Simulate data 
	# for (i in 1:m) {
		# vil <- Reduce(kronecker, rev(v[[i,l]]))
		# vs <- tcrossprod(vil, block.scores[,i,l])
		# x[[i]] <- x[[i]] + array(vs, dimx[[i]])
	# }
# }

# list(x = x, v = v, block.scores = block.scores)
# }


#####################################
# Function to globally orthogonalize 
# sets of rank-1 tensors
#####################################

# Given rank-1 tensors v(i,l), i = 1, ..., m, and l = 1, ..., r,
# such that for each i, v(i,1), ..., v(i,r) have same dimensions,
# find scaling factors a(i,l) such that the tensors u(i,l) = a(i,l) v(i,l)
# satisfy (1/m) sum(i=1:m) < u(i,k), u(i,l) > = 0 if k != l, = 1 if k = l.   


orthogonalize.global <- function(v)
{
m <- NROW(v)
r <- NCOL(v)
if (r > 1) {
	cp <- aperm(tnsr.rk1.cp(v), c(3, 1, 2))
	a <- matrix(, m, r)
	for (l in 1:r) {
		if (l == 1) {
			a[,1] <- 1 / sqrt(mean(cp[,1,1]))
		} else {
			b <- cp[,1:(l-1),l] * a[,1:(l-1)]
			a[,l] <- qr.Q(qr(b), complete = TRUE)[,l]
			a[,l] <- 1 / sqrt(mean((a[,l] * cp[,l,l])^2))
		}
		for (i in 1:m) 
			v[[i,l]][[1]] <- a[i,l] * v[[i,l]][[1]]
	}
}
scale.v(v, "norm", "global", check.args = FALSE)
}




###########################
# Cosine distance function 
# between rank-1 tensors
###########################

cosine.dist <- function(v1, v2)
{
stopifnot(is.list(v1) && is.list(v2))
stopifnot(length(v1) == length(v2))
len <- length(v1)
out <- numeric(len)
nrmfun <- function(x) sqrt(sum(x^2))
cpfun <- function(x, y) sum(x * y)
for (i in 1:len) {
	nrm1 <- prod(sapply(v1[[i]], nrmfun))
	nrm2 <- prod(sapply(v2[[i]], nrmfun))
	if (nrm1 == 0 || nrm2 == 0) next
	cp12 <- prod(mapply(cpfun, v1[[i]], v2[[i]]))
	out[i] <- 1 - (cp12 / nrm1 / nrm2)
}
dim(out) <- dim(v1)
out
}



###########################
# Cosine distance function 
# between rank-1 tensors 
# taking sign into account
###########################

cosine.dist.mod <- function(v1, v2)
{
distance <- cosine.dist(v1, v2)
idx <- which(colMeans(distance) > 1)
if (length(idx) > 0) 
	distance[,idx] <- 2 - distance[,idx]
distance
}


################################################
# Distance between (projectors associated with)
# spaces spanned by two sets of rank-1 tensors
################################################

space.dist <- function(v1, v2)
{
if (!is.matrix(v1)) dim(v1) <- c(length(v1), 1)
if (!is.matrix(v2)) dim(v2) <- c(length(v2), 1)
v1 <- scale.v(v1, check.args = FALSE)
v2 <- scale.v(v2, check.args = FALSE)
m <- nrow(v1)
r <- ncol(v1)
distance <- numeric(m)
for (i in 1:m) {
	pi <- sapply(v1[[i]], length)
	di <- length(pi)
	v1i <- v2i <- matrix(0, prod(pi), r)
	for (l in 1:r) {
		v1i[,l] <- if (di == 1L) {
			v1[[i,l]][[1]] } else Reduce(kronecker, rev(v1[[i,l]]))
		v2i[,l] <- if (di == 1L) {
			v2[[i,l]][[1]] } else Reduce(kronecker, rev(v2[[i,l]]))
	}
	qr1i <- qr(v1i)
	qr2i <- qr(v2i)
	rk1i <- qr1i$rank
	rk2i <- qr2i$rank
	cp <- crossprod(qr.Q(qr1i), qr.Q(qr2i))
	distance[i] <- sqrt(rk1i + rk2i - 2 * sum(cp^2))
}
distance
}	



#################################
# Function to reshape output of 
# 'mgcca_array' of package RGCCA 	
#################################

reshape.mgcca <- function(v)
{
m <- length(v$a)
r <- ncol(v$a[[1]])
vout <- vector("list", m * r)
dim(vout) <- c(m, r)
for (i in 1:m) {
	is.null.b <- is.null(v$b[[i]])
	is.null.c <- is.null(v$c[[i]])
	for (l in 1:r) {
		vout[[i,l]] <- if (is.null.b) {
			list(v$a[[i]][,l])	
		} else if (is.null(c)) {
			list(v$b[[i]][,l])
		} else {
			list(v$b[[i]][,l], v$c[[i]][,l])
		}	
	}
}
vout	
}



######################################
# Function to reshape output of 'ctd'  	
#####################################

reshape.ctd <- function(v)
{
m <- length(v)
r <- ncol(v[[1]][[1]])
vout <- vector("list", m * r)
dim(vout) <- c(m, r)
for (i in 1:m) {
	di <- length(v[[i]]) - 1L
	for (l in 1:r) {
		vout[[i,l]] <- if (r == 1L) {
			lapply(v[[i]][1:d[i]], as.vector)
		} else {
			lapply(v[[i]][1:d[i]], function(x) x[,l])
		}
	}
}
vout		
}



