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

simulate.factor.model <- function(dimx, r, score.cov, noise.cov, 
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

## Preprocess score covariance
SIGMA <- score.cov
if (is.vector(SIGMA) && is.numeric(SIGMA)) {
	SIGMA <- rep_len(SIGMA, r)
	stopifnot(all(SIGMA >= 0))
} else if (is.list(SIGMA)) {
	SIGMA <- rep_len(SIGMA, r)
	for (l in 1:r)
		stopifnot(identical(dim(SIGMA[[l]]), c(m,m)))
	SIGMA <- array(unlist(SIGMA), c(m, m, r))
}
stopifnot(is.numeric(SIGMA) || identical(dim(SIGMA), c(m,m,r)))

## Preprocess noise covariance
PSI <- noise.cov
if (is.numeric(PSI)) {
	PSI <- rep_len(PSI, r)
} else if (is.list(PSI)) {
	stopifnot(length(PSI) == r)
	for (l in 1:r) 
		stopifnot(identical(dim(PSI[[l]]), c(pp[l],pp[l])))
}

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
	score <- NULL
} else if (is.vector(SIGMA)) {
	score <- sapply(sqrt(SIGMA), function(sig) rep(rnorm(n, 0, sig), m))
	dim(score) <- c(n, m, r)
} else {
	score <- array(0, dim = c(n, m, r))
	for (l in 1:r) {
		eig <- eigen(SIGMA[,,l], TRUE)
		stopifnot(all(eig$values >= -(1e-14)))
		k <- sum(eig$values > 0)
		if (k == 0) next
		score[,,l] <- matrix(rnorm(n*k), n, k) %*% 
			(t(eig$vectors[,1:k]) * sqrt(eig$values[1:k]))
	}
}

for (i in 1:m) {	
	## Simulate signal 
	signal <- 0 
	if (r > 0) {
		for (l in 1:r) {
			vil <- Reduce(kronecker, rev(v[[i,l]]))
			signal <- signal + tcrossprod(vil, score[,i,l])
		}
		dim(signal) <- c(p[[i]], n)
	}	
	## Simulate noise
	if (is.numeric(PSI)) {
		noise <- rnorm(pp[i] * n, sd = sqrt(PSI[i]))	
	} else {		
		eig <- eigen(PSI[[i]], TRUE)
		stopifnot(all(eig$values >= -(1e-14)))
		k <- sum(eig$values > 0)
		noise <- if (k == 0) {
			numeric(pp[i] * n) } else {
			matrix(rnorm(n*k), n, k) %*% 
				(t(eig$vectors[,1:k]) * sqrt(eig$values[1:k]))
		}
	}
	dim(noise) <- c(p[[i]], n)		
	x[[i]] <- signal + noise
}

list(x = x, v = v, score = score)
}

