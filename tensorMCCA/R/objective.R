#####################################
# Function to calculate image scores
# y_it = <v_i, X_it> 
#####################################

image.scores <- function(x, v)
{
m <- length(x) # number of datasets
dimx <- lapply(x, dim)
n <- tail(dimx[[1]], 1) # number of instances per dataset
d <- sapply(dimx, length) - 1 # number of dimensions in images/features 
y <- matrix(0, n, m) # image/individual scores on canonical components
for (i in 1:m) { 		
	iprod <- x[[i]]
	pi <- dimx[[i]]
	if (d[i] == 1) { # 1D case
		y[,i] <- crossprod(iprod, v[[i]][[1]]) 
	} else if (d[i] == 2) { # 2D case
		dim(iprod) <- c(pi[1], pi[2] * n)
		iprod <- crossprod(v[[i]][[1]], iprod)
		dim(iprod) <- c(pi[2], n) 
		y[,i] <- crossprod(iprod, v[[i]][[2]]) 
	} else if (d[i] == 3) { # 3D case
		dim(iprod) <- c(pi[1], prod(pi[2:4]))
		iprod <- crossprod(v[[i]][[1]], iprod)
		dim(iprod) <- c(pi[2], pi[3] * n)
		iprod <- crossprod(v[[i]][[2]], iprod)
		dim(iprod) <- c(pi[3], n) 
		y[,i] <- crossprod(iprod, v[[i]][[3]])
	} else warning("Function not yet supported for tensors of order 4+")
}
y
}



###########################################
# Internal function to calculate objective
###########################################

# This function calculates the -weighted- sum of covariances.
# For speed, it does not check the dimensions of arguments. 
# To calculate the sum of correlations with this function, 
# the canonical vectors should be scaled beforehand to have 
# unit variance, e.g., with function 'scale.v'.

objective.internal <- function(x, v, c = 1)
{
n <- tail(dim(x[[1]]), 1) 
y <- image.scores(x, v)
sum(c * crossprod(y)) / n
}



############################################
# Function to calculate sum of correlations
############################################


objective.cor <- function(x, v, c = 1)
{
# Check arguments
test <- check.arguments(x, v)

## Data dimensions 
m <- length(x)
n <- tail(dim(x[[1]]), 1)

## Number of canonical components 
r <- NCOL(v)

## Objective weights
if (length(c) == 1) {
	c <- matrix(1/m^2, m, m)
} else { 
	c <- (c + t(c)) / (2 * sum(c))
}

## Scale canonical vectors so components have unit variance
if (r == 1) {
	y <- image.scores(x, v)
	y <- y - rowMeans(y)
	nrm <- sqrt(colMeans(y^2))
	nrm[nrm < 1e-15] <- 1 
	y <- scale(y, center = FALSE, scale = nrm)
	out <- sum(c * crossprod(y)) / n
} else {
	out <- numeric(r)
	for (k in 1:r) {
		y <- image.scores(x, v[,k])
		y <- y - rowMeans(y)
		nrm <- sqrt(colMeans(y^2))
		nrm[nrm < 1e-15] <- 1
		y <- scale(y, center = FALSE, scale = nrm)
		out[k] <- sum(c * crossprod(y)) / n
	}
}
sum(out)
}




############################################
# Function to calculate sum of covariances
############################################


objective.cov <- function(x, v, c = 1)
{
## Check arguments
test <- check.arguments(x, v)

## Data dimensions 
m <- length(x)
n <- tail(dim(x[[1]]), 1)

## Number of canonical components 
r <- NCOL(v)

## Objective weights
if (length(c) == 1) {
	c <- matrix(1/m^2, m, m)
} else { 
	c <- (c + t(c)) / (2 * sum(c))
}

## Calculate canonical scores and objective
if (r == 1) {
	y <- image.scores(x, v)
	y <- y - rowMeans(y)
	out <- sum(c * crossprod(y)) / n
} else {
	out <- numeric(r)
	for (k in 1:r) {
		y <- image.scores(x, v[,k])
		y <- y - rowMeans(y)
		out[k] <- sum(c * crossprod(y)) / n
	}
}
sum(out)
}

