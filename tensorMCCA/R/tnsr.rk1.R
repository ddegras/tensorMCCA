## INTERNAL FUNCTIONS

##############################################
# Function to approximate a tensor of general 
# order and rank by a rank-1 tensor
##############################################
 

tnsr.rk1 <- function(x, scale = FALSE, init = NULL, maxit = 100, 
	tol = 1e-6)
{
stopifnot(is.numeric(x))
d <- max(1L, length(dim(x))) 

## Trivial case: null tensor
if (all(x == 0)) {
	p <- if (is.vector(x)) length(x) else dim(x)
	return(lapply(p, numeric))
}	

## 1D case
if (d == 1L) {
	nrmv <- if (scale) sqrt(sum(x^2)) else 1
	v <- list(x / nrmv)
	return(v)
}

## 2D case
if (d == 2L) {
	svdx <- if (all(dim(x) > 2)) svds(x, 1) else svd(x, 1, 1)
	s <- if (scale) 1 else sqrt(svdx$d[1])
	v <- list(svdx$u * s, svdx$v * s) 
	return(v)
}

## 3D case
if (d == 3L) 
	return(tnsr3d.rk1(x, scale, init, maxit, tol))

## Initialization
if (is.null(init)) {
	svdx <- hosvd(x, 1) 
	v <- svdx$factors
	if (!scale) {
		s <- as.numeric(svdx$core)^(1/d)
		v <- lapply(v, "*", y = s) # vector balancing
	}
} else {
	v <- init
	nrmv <- sapply(v, function(x) sqrt(sum(x^2)))
	p <- dim(x)
	if (any(nrmv == 0)) {
		v <- lapply(p, numeric) 
	} else {
		if (!scale) nrmv <- nrmv / prod(nrmv)^(1/d) # vector balancing
		v <- mapply("/", x = v, y = nrmv, SIMPLIFY = FALSE)
	}
}
objective <- tnsr.vec.prod(x, v, 1:d)
iters <- if (maxit == 0) NULL else (1:maxit)
for (it in iters) {
	objective.prev <- objective
	for (k in 1:d) {
		v[[k]] <- tnsr.vec.prod(x, v[-k], (1:d)[-k])
		nrmv <- sqrt(sum(v[[k]]^2))
		v[[k]] <- v[[k]] / nrmv
	}
	objective <- nrmv^2
	if (objective <= (1+tol) * objective.prev) break
}

if (scale) return(v)
s <- nrmv^(1/d)
v <- lapply(v, "*", y = s)
v	
}




##############################################
# Function to approximate a general tensor  
# by a rank-1 tensor under norm constraints 
# on scores
##############################################
 

tnsr.rk1.score <- function(x, cnstr, maxit = 100, tol = 1e-6)
{
if (all(x == 0)) {
	dimx <- if (is.vector(x)) length(x) else dim(x)
	return(lapply(dimx, numeric))
}

d <- if (is.vector(x)) 1L else length(dim(x))

gradfun <- function(v, x) {
	d <- length(v)
	if (d == 1)	return(list(v[[1]] - x))
	out <- vector("list", d)
	nrmv <- sapply(v, function(vv) sum(vv^2)) 
	for (k in 1:d) 
		out[[k]] <- prod(nrmv[-k]) * v[[k]] - 
			tnsr.vec.prod(x, v[-k], (1:d)[-k])
	out	
}

objfun <- function(a, v, grad, x, cnstr, returnv = FALSE) {
	d <- length(v)
	vnew <- vector("list", d)
	for (k in 1:d) 
		vnew[[k]] <- v[[k]] - a * grad[[k]]
	vnew <- scale.v(vnew, "var", x = cnstr,	check.args = FALSE)
	kronv <- outer.prod.nodim(vnew)
	# kronv <- if (d == 1) { 
		# vnew[[1]]
	# } else {
		# as.vector(Reduce(kronecker, rev(vnew))) 
	# }
	out <- if (returnv) {
		list(objective = sum((x - kronv)^2), v = vnew)
	} else {
		sum((x - kronv)^2)
	}
	out
}

v <- tnsr.rk1(x)
objective <- objfun(0, v, numeric(d), x, cnstr, TRUE)
v <- objective$v
objective <- objective$objective
agrid <- c(10^seq(-4, 1))
for (it in 1:maxit) {
	objective.old <- objective
	grad <- gradfun(v, x)
	objective.best <- objective
	vbest <- v
	for (a in agrid) {
		objective <- objfun(a, v, grad, x, cnstr, TRUE)
		if (objective$objective < objective.best) {
			objective.best <- objective$objective
			vbest <- objective$v
		}
 	}	
 	objective <- objective.best
 	v <- vbest
 	progress <- objective.old - objective 
 	if (progress <= (tol * max(1, objective.old))) break
 }
	
v
}



##############################################
# Function to approximate a general 3D tensor 
# by a rank-1 tensor
##############################################

tnsr3d.rk1 <- function(x, scale = FALSE, init = NULL, maxit = 100, 
	tol = 1e-6)
{
stopifnot(is.array(x) && length(dim(x)) == 3)
p <- dim(x)
if (all(x == 0)) {
	if (scale) {
		return(lapply(p, function(xx) rep(1/sqrt(xx), xx)))
	} else {	
		return(lapply(p, numeric))
	}
}
maxit <- as.integer(maxit)
stopifnot(maxit >= 0)
scalefun <- function(x) x / sqrt(sum(x^2))

## Initialization 
if (is.null(init)) {
	## Rank-1 HOSVD
	v <- vector("list", 3)
	dim(x) <- c(p[1], p[2] * p[3]) 
	svdx <- if (all(dim(x) > 2)) {
		svds(x, k = 1) 
	} else { 
		svd(x, nu = 1, nv = 1)
	}
	dim(x) <- p
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
	if (!scale) 
		v <- lapply(v, "*", y = (s1*s2)^(1/3))
} else {
	stopifnot(length(init) == 3)
	v <- init
	nrmv <- sapply(v, function(x) sqrt(sum(x^2)))
	s <- if (scale) nrmv else nrmv / prod(nrmv)^(1/3)
	s[s == 0] <- 1
	v <- mapply("/", x = v, y = nrmv, SIMPLIFY = FALSE)
}
objective <- tnsr.vec.prod(x, v, 1:3)

iters <- if (maxit == 0) NULL else (1:maxit)
for (it in iters) {
	objective.prev <- objective	
	dim(x) <- c(p[1] * p[2], p[3])
	xv <- x %*% v[[3]]
	dim(xv) <- c(p[1], p[2])
	v[[1]] <- scalefun(xv %*% v[[2]])
	v[[2]] <- scalefun(crossprod(xv, v[[1]]))
	dim(x) <- c(p[1], p[2] * p[3])
	xv <- crossprod(x, v[[1]])
	dim(xv) <- p[2:3]
	v[[3]] <- crossprod(xv, v[[2]]) 
	objective <- sum(v[[3]]^2)
	nrmv <- sqrt(objective)
	v[[3]] <- v[[3]] / nrmv
	if (objective < (1+tol) * objective.prev) break
}
s <- ifelse(scale, 1, nrmv^(1/3)) 
v[[1]] <- as.vector(v[[1]]) * s 
v[[2]] <- as.vector(v[[2]]) * s 
v[[3]] <- as.vector(v[[3]]) * s 

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
if (is.null(ortho)) return(v0)
stopifnot((NROW(ortho) == length(v0)))
if (!is.matrix(ortho)) 
	dim(ortho) <- c(length(ortho), 1)
if (is.list(ortho[[1]]))
	ortho <- tnsr.rk1.expand(ortho)
for (i in 1:m) v0[[i]] <- lapply(v0[[i]], as.matrix)
v <- v0
cpfun <- function(x, y) sum(x * y) 

for (i in 1:m) {
objective <- Inf
cp.v0v <- mapply(cpfun, v0[[i]], v[[i]])
cp.vv <- mapply(cpfun, v[[i]], v[[i]])
pi <- sapply(v0[[i]], length)
vzero <- all(unlist(v) == 0)
if (vzero) next
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
			if (rkM == pi[k]) {
				v[[i]] <- lapply(pi, numeric)
				vzero <- TRUE
				break
			} else if (rkM == 0) {
				projv0ik <- v0[[i]][[k]]
			} else {
				Qik <- qr.Q(qrM)[,1:rkM]
				projv0ik <- v0[[i]][[k]] - 
					Qik %*% crossprod(Qik, v0[[i]][[k]])
			} 
			## Update v(i,k) and cross-products
			v[[i]][[k]] <- prod(cp.v0v[-k]) / prod(cp.vv[-k]) * projv0ik
			cp.v0v[k] <- cpfun(v0[[i]][[k]], v[[i]][[k]])
			cp.vv[k] <- sum(v[[i]][[k]]^2)
		}
		if (vzero) break	
		objective <- sum((unlist(v0[[i]]) - unlist(v[[i]]))^2)
		delta <- abs(objective - objective.old)
		if (it > 1 && delta <= tol * max(1, objective.old) ) break
	}	
}
v
}





# tnsr.rk1.ortho.global <- function(v0, ortho, maxit, tol)
# {
# if (is.null(ortho)) 
	# return(scale.v(v0, type = "norm", scale = "global", check.args = FALSE))
# m <- length(v0)
# d <- sapply(v0, length)
# maxd <- max(d)
# p <- relist(lapply(unlist(v0, FALSE), length), v0)
# pp <- sapply(p, prod)
# if (!is.matrix(ortho)) dim(ortho) <- c(length(ortho), 1)
# northo <- ncol(ortho)
# ortho.arr <- tnsr.rk1.expand(ortho)
# v0.arr <- tnsr.rk1.expand(v0)

# v <- scale.v(v0, check.args = FALSE)


# bdiagfun <- function(x) {
# if (is.matrix(x)) return(x)
# m <- length(x)
# if (m == 1) return(x[[1]])
# nr <- sapply(x, nrow)
# nc <- sapply(x, ncol)
# endr <- cumsum(nr)
# startr <- c(1, endr[-m] + 1)
# endc <- cumsum(nc)
# startc <- c(1, endc[-m] + 1)
# out <- matrix(0, sum(nr), sum(nc))
# for (i in 1:m) 
	# out[startr[i]:endr[i], startc[i]:endc[i]] <- x[[i]]
# out
# }


# for (it in 1:maxit) {
		
# for (k in 1:maxd) {
	
	
	# nrmv <- tnsr.rk1.nrm(v, "block") 
		
	
	# for (i in 1:m) {
# objective <- Inf
# cp.v0v <- mapply(cpfun, v0[[i]], v[[i]])
# cp.vv <- mapply(cpfun, v[[i]], v[[i]])


		# objective.old <- objective
		# ## Block coordinate descent

					# ## Calculate matrix of orthogonal constraints for v(i,k)
			# Mik <- if (d[i] == 1) {
				# matrix(unlist(ortho[i,]), ncol = ncol(ortho))
			# } else { sapply(ortho[i,], tnsr.vec.prod, 
				# v = v[[i]][-k], modes = (1:d[i])[-k]) }
			# ## Apply orthogonal constraints to target vector
			# qrM <- qr(Mik)
			# rkM <- qrM$rank
			# Qik <- qr.Q(qrM)[,1:rkM]
			# projv0ik <- if (rkM > 0) {
				# v0[[i]][[k]] - Qik %*% crossprod(Qik, v0[[i]][[k]])
			# } else { v0[[i]][[k]] }
			# ## Update v(i,k) and cross-products
			# v[[i]][[k]] <- prod(cp.v0v[-k]) / prod(cp.vv[-k]) * projv0ik
			# cp.v0v[k] <- cpfun(v0[[i]][[k]], v[[i]][[k]])
			# cp.vv[k] <- sum(v[[i]][[k]]^2)
		# }	
		# objective <- sum((unlist(v0[[i]]) - unlist(v[[i]]))^2)
		# delta <- abs(objective - objective.old)
		# if (it > 1 && delta <= tol * max(1, objective.old) ) break
	# }	
# }
# v
# }







	

