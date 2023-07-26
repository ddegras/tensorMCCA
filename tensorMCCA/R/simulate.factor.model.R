#########################################
# Function to simulate canonical weights 
#########################################

simulate.v <- function(p, r, scale = 1,  
	ortho = c("cyclical", "random", "alldim", "maxdim", "none")) 
{
## Check input arguments
stopifnot(all(unlist(p) >= 1))
stopifnot(length(r) == 1 && r >= 0)
stopifnot(length(scale) %in% c(1, length(p)))
stopifnot(all(scale >= 0))

## Preprocess input arguments
if (!is.list(p)) p <- list(p)
p <- lapply(p, as.integer)
d <- sapply(p, length)
m <- length(p)
r <- as.integer(r)
if (r == 0) {
	r <- 1
	scale <- numeric(m)
}
seq2r <- if (r > 1) 2:r else NULL
if (r == 1) ortho == "none"
if (length(scale) < m) scale <- rep(scale, m)
ortho.mode <- NULL
if (is.numeric(ortho)) {
	ortho <- rep_len(as.integer(ortho), m)
	stopifnot(all(ortho > 0))
} else {
	ortho <- match.arg(ortho)
	if (ortho %in% c("cyclical", "random")) {
		x <- lapply(p, function(dims) array(0L, dim = c(dims, 1)))
		ortho.mode <- set.ortho.mode(x, r, method = ortho)
	}
}

## Error message
err <- paste("Orthogonality constraints cannot be accommodated.",
	"Consider using a different type of constraint in 'ortho'",
	"or setting a smaller value for 'r'.")

## Function to generate random vectors with unit norm 
randvec <- function(p, r = 1) {
	if (r == 1) {
		x <- runif(p,-1, 1)
		return(x / sqrt(sum(x^2)))
	}
	x <- matrix(runif(p * r,-1, 1), p, r)
	sweep(x, 2, sqrt(colSums(x^2)), "/")
}

## Output 
v <- vector("list", m * r)
dim(v) <- c(m, r)

## Main loop
for (i in 1:m) {
	if (scale[i] == 0) {
		v[i,] <- replicate(r, lapply(p[[i]], numeric), FALSE)
		next
	}
	if (!is.null(ortho.mode)) {
		v[[i,1]] <- lapply(p[[i]], randvec)
		for (l in seq2r) {
			for (k in 1:d[i]) {
				idx <- which(unlist(ortho.mode[i,1:(l-1),l]) == k)
				if (length(idx) == 0) {
					v[[i,l]][[k]] <- randvec(p[[i]][k])
				} else {
					vecs <- sapply(v[i,idx], "[[", k)
					qrvecs <- qr(vecs)
					if (qrvecs$rank == p[[i]][k]) stop(err)
					v[[i,l]][[k]] <- qr.Q(qrvecs, TRUE)[,qrvecs$rank + 1]
				}
			}	
		}		
	} else {
		kk <- if (is.numeric(ortho)) { 
			ortho[i] 
		} else if (ortho == "alldim") {
			1:d[i]
		} else if (ortho == "maxdim") {
			which.max(p[[i]])
		} else if (ortho == "none") {
			NULL
		}
		for (k in 1:d[i]) {
			vik <- randvec(p[[i]][k], r)
			if (k %in% kk) {
				if (r > p[[i]][k]) stop(err)
				vik <- qr.Q(qr(vik))
			}
			for (l in 1:r) v[[i,l]][[k]] <- vik[,l]			
		}
	} 
	
	## Rescale tensors if needed 
	if (scale[i] != 1) {
		for (l in 1:r) 
			v[[i,l]][[1]] <- v[[i,l]][[1]] * scale[i] 
	}
}
v
}


################################
# Function to simulate all data 
################################

simulate.factor.model <- function(dimx, r, scale.v = 1, score.cov = NULL, 
	noise.cov = NULL, ortho.v = c("cyclical", "random", "alldim", "maxdim", "none"))
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
r <- as.integer(r)
if (is.character(ortho.v)) ortho.v <- match.arg(ortho.v)

## Preprocess score covariance
if (is.null(score.cov)) score.cov <- 1
SIGMA <- score.cov
if (is.vector(SIGMA) && is.numeric(SIGMA)) {
	SIGMA <- rep_len(SIGMA, r)
	stopifnot(all(SIGMA >= 0))
} else if (is.list(SIGMA)) {
	SIGMA <- rep_len(SIGMA, r)
	for (l in 1:r) {
		if (!is.matrix(SIGMA[[l]])) 
			SIGMA[[l]] <- as.matrix(SIGMA[[l]])
		stopifnot(all(dim(SIGMA[[l]]) == m))	
	}	
	SIGMA <- array(unlist(SIGMA), c(m, m, r))
}
# stopifnot(is.numeric(SIGMA) || identical(dim(SIGMA), c(m,m,r)))

## Preprocess noise covariance
if (is.null(noise.cov)) noise.cov <- 1
PSI <- noise.cov
if (is.numeric(PSI)) {
	PSI <- rep_len(PSI, m)
} else if (is.list(PSI)) {
	stopifnot(length(PSI) == m)
	for (i in 1:m) {
		if (!is.matrix(PSI[[i]])) 
			PSI[[i]] <- as.matrix(PSI[[i]])
		stopifnot(all(dim(PSI[[i]]) == pp[i]))
	}	
}

## Outputs
x <- vector("list", m)
if (r == 0) {
	v <- vector("list", m)
	for (i in 1:m) v[[i]] <- lapply(p[[i]], numeric)
} else {
	v <- simulate.v(p, r, scale.v, ortho.v)
}
if (r == 0) {
	score <- array(0, c(n, m, 1))
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

