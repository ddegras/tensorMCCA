#####################################
# Function to calculate image scores
# y_it = <v_i, X_it> 
#####################################

image.scores <- function(x, v)
{
m <- length(x) # number of datasets
dimx <- lapply(x, dim)
n <- tail(dimx[[1]], 1) # number of instances per dataset
ndim.img <- sapply(dimx, length) - 1 # number of dimensions in images/features 
y <- matrix(0, n, m) # image/individual scores on canonical components
for (i in 1:m) { 		
	iprod <- x[[i]]
	if (ndim.img[i] == 1) { # 1D case
		y[,i] <- crossprod(iprod, v[[i]][[1]]) 
	} else if (ndim.img[i] == 2) { # 2D case
		dim(iprod) <- c(dimx[[i]][1], dimx[[i]][2] * n)
		iprod <- crossprod(v[[i]][[1]], iprod)
		dim(iprod) <- c(dimx[[i]][2], n) 
		y[,i] <- crossprod(iprod, v[[i]][[2]]) 
	} else { # 3D case
		dim(iprod) <- c(dimx[[i]][1], prod(dimx[[i]][2:4]))
		iprod <- crossprod(v[[i]][[1]], iprod)
		dim(iprod) <- c(dimx[[i]][2], dimx[[i]][3] * n)
		iprod <- crossprod(v[[i]][[2]], iprod)
		dim(iprod) <- c(dimx[[i]][3], n) 
		y[,i] <- crossprod(iprod, v[[i]][[3]])
	}
}
return(y)
}




#########################################
# Function for scaling canonical vectors
#########################################

scale.v <- function(v, x = NULL, type = c("norm", "var"), 
	scope = c("block", "global"), balance = TRUE)
{
type <- match.arg(type)   # norm or variance constraints
scope <- match.arg(scope) # block or global constraints
m <- length(v)
ndimv <- sapply(v, length)

## If variance constraints, calculate variances and rescale
if (type == "var") {
	y <- image.scores(x, v)
	sdy <- sqrt(colMeans(y^2))
	for (i in 1:m) {
		if (sdy[i] == 0) next
		for (k in 1:ndimv[i]) 
			v[[i]][[k]] <- v[[i]][[k]] / sdy[i]^(1/ndimv[i])
	}
	if (!balance) return(v)	
}

## Calculate (Euclidean) norms of canonical vectors
nrmv <- vector("list",m)
nrmfun <- function(vv) sqrt(sum(vv^2))
for (i in 1:m) 
	nrmv[[i]] <- sapply(v[[i]], nrmfun)
	
## If variance constraints and balancing requirement,
## make norms of canonical vectors equal in each block
## without changing their product  
if (type == "var") {
	for (i in 1:m) {
		if (any(nrmv[[i]] == 0)) next
		s <- (prod(nrmv[[i]]))^(1/ndimv[i]) / nrmv[[i]]
		for (k in 1:ndimv[i]) 
			v[[i]][[k]] <- s[k] * v[[i]][[k]] 
	}
	return(v)
} 

## If Euclidean norm constraints at the block level,
## rescale each canonical vector by its norm
if (type == "norm" && scope == "block") {
	for (i in 1:m) {
		for (k in 1:ndimv[i]) {
			if (nrmv[[i]][k] > 0) 
				v[[i]][[k]] <- v[[i]][[k]] / nrmv[[i]][k]
		}
	}
	return(v)	
} 

## If Euclidean norm constraints at the global level,
## rescale each canonical vector in a given mode/dimension
## by the global norm of all canonical vectors in this mode
if (type == "norm" && scope == "global") {
	for (k in 1:max(ndimv)) {
		idx <- which(ndimv >= k)
		sum.nrmv <- sqrt(sum((sapply(nrmv[idx], "[", k))^2))
		if (sum.nrmv == 0) next
		for (i in 1:m) 	
			v[[i]][[k]] <- v[[i]][[k]] / sum.nrmv	

	}
	return(v)	
} 

## Other cases (should not happen)
return(NULL)
}



###########################################
# Internal function to calculate objective
###########################################

# This function calculates the -weighted- sum of covariances.
# For speed, it does not check the dimensions of arguments. 
# To calculate the sum of correlations with this funtion, 
# the canonical vectors should be scaled beforehand to have 
# unit variance, .e.g., with function 'scale.v'.

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

## Objective weights
m <- length(v)
if (length(c) == 1) c <- matrix(c,m,m)
c <- c / sum(c)

## Center the data
for(i in 1:m) {
    mu <- rowMeans(x[[i]], dims = ndim.img[i])
    x[[i]] <- x[[i]] - as.vector(mu)
}

## Number of canonical components 
r <- NCOL(v[[1]][[1]])

## Scale canonical vectors so components have unit variance
if (r == 1) {
	y <- image.scores(x, v)
	sdy <- sqrt(colMeans(y^2))
	sdy[sdy < 1e-15] <- 1 
	out <- sum(c * crossprod(y) / tcrossprod(sdy)) / n
} else {
	out <- numeric(r)
	for (k in 1:r) {
		vk <- vector("list", m)
		for (i in 1:m)
			vk[[i]] <- lapply(v[[i]], function(mat) mat[,k])
		y <- image.scores(x, vk)
		sdy <- sqrt(colMeans(y^2))
		sdy[sdy < 1e-15] <- 1 
		out[k] <- sum(c * crossprod(y) / tcrossprod(sdy)) / n
	}
}
out
}

