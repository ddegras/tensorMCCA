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



simulate.factor.model <- function(dimx, r, sigma2, xi, ortho = FALSE)
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

## Outputs
x <- vector("list", m)
v <- simulate.v(p, r, ortho)
score <- if (r == 0) {
	numeric(n)
} else {
	sapply(sqrt(sigma2), function(sig) rnorm(n, 0, sig))
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
	x[[i]] <- signal + noise
}

list(x = x, v = v, score = score)
}

########################################
# Simulate model (2) with block scores
# under independence assumption on block scores
########################################


# X_it = sum_(l=1:r) v_i^l s_it^l        (2)


# ASSUMPTIONS ON THE SCORES
# with s_it^l, l=1:r, independent for each (i,t)
# but correlation among the s_it, i=1:m (for each (l,t))


# Inputs
# dimx:		dimensions of data tensors (list of length m)
# r:		rank of data decomposition
# sigma:	array of score covariances (m x m x r) 



simulate.x2 <- function(dimx, r, sigma)
{
## Data dimensions
m <- length(dimx)
n <- sapply(dimx, tail, 1)
if (any(n != n[1])) 
	stop(paste("Last values in each component of list",
	"'dimx' must all be equal (number of replications)."))
n <- n[1]
p <- sapply(dimx, function(idx) idx[-length(idx)])
d <- sapply(dimx, length) - 1

## Check dimensions of sigma and replicate if needed
test <- is.matrix(sigma) && dim(sigma) == c(m,m)
stopifnot(test || all(dim(sigma) == c(m,m,r)))
if (test) sigma <- replicate(r, sigma)

## Outputs
v <- vector("list", m * r)
dim(v) <- c(m, r)
x <- lapply(dimx, function(dims) array(0, dims))
block.scores <- array(dim = c(n,m,r))

## Generate canonical vectors
rand.vec <- function(pp) {
	v <- runif(pp)
	v / sqrt(sum(v^2))
}
for (i in 1:m)
for (l in 1:r)
v[[i,l]] <- lapply(p[[i]], rand.vec)

## Generate data
for (l in 1:r) {    
	## Generate block scores
	R <- chol(sigma[,,l])
	tmp <- matrix(rnorm(m*n), m, n)
	block.scores[,,l] <- crossprod(tmp, R)
            	
	## Simulate data 
	for (i in 1:m) {
		vil <- Reduce(kronecker, rev(v[[i,l]]))
		vs <- tcrossprod(vil, block.scores[,i,l])
		x[[i]] <- x[[i]] + array(vs, dimx[[i]])
	}
}

list(x = x, v = v, block.scores = block.scores)
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

