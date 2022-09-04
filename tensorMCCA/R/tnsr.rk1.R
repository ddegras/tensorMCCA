## INTERNAL FUNCTIONS

##############################################
# Function to approximate a tensor of general 
# order and rank by a rank-1 tensor
##############################################
 

tnsr.rk1 <- function(x, scale = FALSE, maxit = 100, tol = 1e-6)
{
stopifnot(is.numeric(x))
stopifnot(maxit >= 0 && tol >= 0)
eps <- 1e-14
if (is.vector(x)) {
	if (scale) {
		nrm <- sqrt(sum(x^2))
		if (nrm > eps) {
			return(list(x / nrm))
		} else {
			return(list(rep.int(1 / sqrt(length(x)), length(x))))
		}
	} else {
		return(list(x))
	}
}
if (all(x == 0)) {
	if (scale) {
		return(lapply(dim(x), 
			function(xx) rep.int(1/sqrt(length(x)), length(x))))		
	} else {
		return(lapply(dim(x), numeric))
	}
}
d <- length(dim(x)) 
if (d == 3L) { return(tnsr3d.rk1(x, scale, maxit, tol)) }
svdx <- hosvd(x, 1) 
nrmv <- if (scale) 1 else as.numeric(svdx$core)
v <- svdx$vectors
if (d == 2L || maxit == 0) {
	if (scale) {
		return(v)
	} else {
		return(lapply(v, "*", y = nrmv^(1/d)))		
	}
}
if (!scale) v[[d]] <- v[[d]] * nrmv
kronv <- Reduce(kronecker, v)
for (it in 1:maxit) {
	nrmv.old <- nrmv
	kronv.old <- kronv
	for (k in 1:d) {
		v[[k]] <- tvec.prod(x, v[-k], (1:d)[-k])
		if (k < d) v[[k]] <- v[[k]] / sqrt(sum(v[[k]]^2))
	}
	if (scale) {
		v[[d]] <- v[[d]] / sqrt(sum(v[[d]]^2))
	} else {
		nrmv <- sqrt(sum(v[[d]]^2))
	}
	kronv <- Reduce(kronecker, v)
	e <- sqrt(sum((kronv - kronv.old)^2))
	if (e <= tol * max(nrmv, nrmv.old, 1)) break
}

if (scale) {
	return(v)
}
s <- c(rep.int(nrmv^(1/d), d - 1), nrmv^(1/d - 1))
for (k in 1:d) v[[k]] <- v[[k]] * s[k] 
v	
}





##############################################
# Function to approximate a general 3D tensor 
# by a rank-1 tensor
##############################################

tnsr3d.rk1 <- function(x, scale = FALSE, maxit = 100, tol = 1e-6)
{
stopifnot(is.array(x) && length(dim(x)) == 3)
p <- dim(x)
if (all(x == 0)) {
	if (scale) {
		return(lapply(p, 
			function(xx) rep.int(1/sqrt(length(xx)), length(xx))))
	} else {	
		return(lapply(p, numeric))
	}
}
v <- vector("list", 3)
maxit <- as.integer(maxit)
stopifnot(maxit >= 0)

## Rank-1 HOSVD
dim(x) <- c(p[1], p[2] * p[3]) 
svdx <- if (all(dim(x) > 2)) {
	svds(x, k = 1) 
} else { 
	svd(x, nu = 1, nv = 1)
} 
v[[1]] <- svdx$u
s1 <- svdx$d
svdx <- if (all(p[2:3] > 2)) {
	svds(matrix(svdx$v, p[2], p[3]), k = 1) 
} else {
	svd(matrix(svdx$v, p[2], p[3]), nu = 1, nv = 1)
}
v[[2]] <- svdx$u
v[[3]] <- svdx$v
s2 <- svdx$d
nrmv <- if (scale) 1 else s1 * s2
if (maxit == 0) {
	if (scale) { return(v) } 
	return(lapply(v, "*", nrmv^(1/3)))
}
if (scale) v[[3]] <- v[[3]] * nrmv
kronv <- Reduce(kronecker, v)

for (it in 1:maxit) {
	kronv.old <- kronv
	nrmv.old <- nrmv	
	dim(x) <- c(p[1] * p[2], p[3])
	xv <- x %*% v[[3]]
	dim(xv) <- c(p[1], p[2])
	v[[1]] <- xv %*% v[[2]] 
	v[[1]] <- v[[1]] / sqrt(sum(v[[1]]^2))
	v[[2]] <- crossprod(xv, v[[1]]) 
	v[[2]] <- v[[2]] / sqrt(sum(v[[2]]^2))
	dim(x) <- c(p[1], p[2] * p[3])
	xv <- crossprod(x, v[[1]])
	dim(xv) <- p[2:3]
	v[[3]] <- crossprod(xv, v[[2]]) 
	nrmv <- sqrt(sum(v[[3]]^2))	
	if (scale) {
		v[[3]] <- v[[3]] / nrmv	
		nrmv <- 1
	}		
	kronv <- Reduce(kronecker, v)
	e <- sqrt(sum((kronv - kronv.old)^2)) 
	if (e <= tol * max(nrmv, nrmv.old, 1)) break
}

s <- nrmv^(c(1/3, 1/3, -2/3))
for (k in 1:3) v[[k]] <- v[[k]] * s[k] 

v	
	
}



#########################################
# Function to approximate rank-1 tensors
# under orthogonality constraints 
#########################################


## Inputs
# v0:	list of rank-1 tensors (targets) each specified as a list of vectors
# ortho: list of tensors or NULL (orthogonality constraints) 

# Calling m the length of 'v0', 'ortho' should be a matrix (of lists) 
# with m rows. The entry ortho[[i,l]] should be an array with same dimensions
# as those specified by v0[[i]], which is a list of vectors.  
 

tnsr.rk1.ortho <- function(v0, ortho, maxit, tol)
{
m <- length(v0)
d <- sapply(v0, length)
v <- v0
cpfun <- function(x, y) sum(x * y) 

for (i in 1:m) {
objective <- Inf
cp.v0v <- mapply(cpfun, v0[[i]], v[[i]])
cp.vv <- mapply(cpfun, v[[i]], v[[i]])

	for (it in 1:maxit) {
		objective.old <- objective
		## Block coordinate descent

		for (k in 1:d[i]) {
			## Calculate matrix of orthogonal constraints for v(i,k)
			Mik <- if (d[i] == 1) {
				matrix(unlist(ortho[i,]), ncol = ncol(ortho))
			} else { sapply(ortho[i,], tnsr.vec.prod, 
				v = v[[i]][-k], modes = (1:d[i])[-k]) }
			## Apply orthogonal constraints to target vector
			qrM <- qr(Mik)
			rkM <- qrM$rank
			Qik <- qr.Q(qrM)[,1:rkM]
			projv0ik <- if (rkM > 0) {
				v0[[i]][[k]] - Qik %*% crossprod(Qik, v0[[i]][[k]])
			} else { v0[[i]][[k]] }
			## Update v(i,k) and cross-products
			v[[i]][[k]] <- prod(cp.v0v[-k]) / prod(cp.vv[-k]) * projv0ik
			cp.v0v[k] <- cpfun(v0[[i]][[k]], v[[i]][[k]])
			cp.vv[k] <- sum(v[[i]][[k]]^2)
		}	
		objective <- sum((unlist(v0[[i]]) - unlist(v[[i]]))^2)
		delta <- abs(objective - objective.old)
		if (it > 1 && delta <= tol * max(1, objective.old) ) break
	}	
}
v
}







	

