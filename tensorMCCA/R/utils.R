##############################
# Modified dimension function 
##############################

# Consider vector of length n as 1-by-n matrix

dimfun <- function(x) if (is.vector(x)) length(x) else dim(x)
# {
# dimx <- dim(x)
# if (length(dimx) <= 1) dimx <- c(1,length(x)) 
# dimx
# }


################################################
# Function to calculate norms of rank-1 tensors
################################################


tnsr.rk1.nrm <- function(v, norm = c("block", "global"))
{
stopifnot(is.list(v))
stopifnot(is.vector(v) || is.matrix(v))
norm <- match.arg(norm)
lenv <- length(v) 
out <- numeric(lenv)
sqnrmfun <- function(x) sum(x^2)
for (i in 1:lenv) 
	out[i] <- prod(sapply(v[[i]], sqnrmfun))
dim(out) <- dim(v)
if (norm == "block") return(sqrt(out))
if (is.vector(out)) return(sqrt(mean(out))) 
sqrt(colMeans(out))
}


#######################################
# Function to calculate cross-products 
# between rank-1 tensors
#######################################


tnsr.rk1.cp <- function(v1, v2 = NULL)
{
stopifnot(is.list(v1))
stopifnot(is.null(v2) || is.list(v2))
stopifnot(is.null(v2) || length(v1) == length(v2))

cpfun <- function(x, y) sum(x * y)
if (is.null(v2)) {
	m <- NROW(v1)
	r <- NCOL(v1)
	if (!is.matrix(v1)) dim(v1) <- c(m, r)
	out <- array(dim = c(r, r, m))
	for (i in 1:m) {
		for (k in 1:r) {
			for (l in 1:k) {
				out[k,l,i] <- prod(mapply(cpfun, v1[[i,k]], v1[[i,l]]))
				out[l,k,i] <- out[k,l,i]
			}
		}
	}
} else {
	m <- length(v1)
	out <- numeric(m)
	for (i in 1:m)
		out[i] <- prod(mapply(cpfun, v1[[i]], v2[[i]]))
	dim(out) <- dim(v1)
}
out
}



################################
# Basic functions to calculate 
# outer products
################################

# faster than function(v) Reduce(kronecker, rev(v))

## Outer product of vectors in list with vector(ized) output
outer.prod.nodim <- function(v) {
d <- length(v)
if (d == 1) return(v[[1]])
out <- v[[1]]
for (k in 2:d) {
	out <- tcrossprod(out, v[[k]])
	dim(out) <- NULL
}
out
}

## Regular outer product of vectors
outer.prod <- function(v) {
stopifnot(is.list(v)) 
d <- length(v)
if (d == 1) return(v[[1]])
out <- v[[1]]
for (k in 2:d) {
	out <- tcrossprod(out, v[[k]])
	dim(out) <- NULL
}
dim(out) <- sapply(v, length)
out
}


####################################
# Function to expand rank-1 tensors 
# from lists of vectors to arrays
####################################


tnsr.rk1.expand <- function(v)
{
stopifnot(is.list(v))
if (is.numeric(v[[1]])) 
	return(outer.prod(v))
out <- lapply(v, outer.prod)
dim(out) <- dim(v)
out
}


########################################
# Function to calculate cosine distance 
# between rank-1 tensors
########################################


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


###################################################
# Function to calculate cosine distance 
# between rank-1 tensors taking sign into account
###################################################


# cosine.dist.mod <- function(v1, v2)
# {
# distance <- cosine.dist(v1, v2)
# if (is.matrix(distance)) {
	# idx <- which(colMeans(distance) > 1)
	# if (length(idx) > 0) 
		# distance[,idx] <- 2 - distance[,idx]
# } else {
	# idx <- which(distance > 1)
	# if (length(idx) > 0) 
		# distance[idx] <- 2 - distance[idx]
# }
# distance
# }


#####################################
# Wrapper function to calculate best 
# rank-1 approximation to tensor 
#####################################


# Inputs
# redx: reduced tensor (possibly projected in lower-dimensional space 
#    to handle orthogonality constraints, singleton dimensions dropped
# dimx: original dimensions of x
# 
tnsr.rk1.wrapper <- function(redx, dimx = NULL)
{
## Data dimensions
dimredx <- dim(redx)
if (is.null(dimredx)) { 
	n <- length(redx)
	dimredx <- c(1, n)
	d <- 1
} else {
	d <- length(dimredx)
	n <- dimredx[d]
} 
if (is.null(dimx)) dimx <- dimredx

## Calculate best approximation (canonical weights)
if (d == 1) {
	v <- list(1)
} else if (d == 2) {
	svdx <- tryCatch(svds(redx, 1, 1, 0), error = function(e) NULL)
	if (is.null(svdx)) svdx <- svd(redx, 1, 0)
	v <- list(svdx$u)
} else if (d == 3) {
	v <- tnsr3d.rk1(redx)[1:2]
} else {
	v <- tnsr.rk1(redx)[1:(d-1)]
}

## Reshape tensor as needed
dimv <- dimx[-length(dimx)]
if (any(dimv == 1))	{
	vtmp <- replicate(length(dimv) - 1, 1, FALSE)
	vtmp[dimv > 1] <- v			
	v <- vtmp
}

v	
}


#########################################
# Function to reorient canonical weights 
# so as to maximize objective function
#########################################


reorient <- function(score, w)
{
m <- NROW(w)
flip <- logical(m)
if (m == 1) 
	return(list(flip = flip, objective = sum(w * score^2)))
x <- crossprod(score) * w
if (all(x >= 0)) 
	return(list(flip = flip, objective = sum(x)))
d <- diag(x)
csum <- colSums(x) - d
while(any(csum < 0))
{
	i <- which.min(csum)
	x[,i] <- -x[,i] 
	x[i,] <- x[,i]
	x[i,i] <- -x[i,i] 
	flip[i] <- !(flip[i])
	csum <- colSums(x) - d
}	
list(flip = flip, objective = sum(x))
}

