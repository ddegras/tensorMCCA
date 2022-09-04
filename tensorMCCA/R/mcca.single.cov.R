## INTERNAL FUNCTIONS

####################################################
# Maximize sum of covariances 
# under scaling constraints on canonical weights
# and orthogonality constraints on canonical scores 
# All constraints are at the block level
####################################################

	
mcca.cov.block.score <- function(x, v, w, ortho, sweep, 
	maxit, tol, verbose)
{
	
## Data dimensions
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
m <- length(x)
n <- tail(dimx[[1]], 1)

objective <- numeric(maxit + 1L)
objective[1] <- objective.internal(x, v, w)
if (verbose) 
	cat("\nIteration", 0, "Objective", objective[1])
vbest <- v
objective.best <- objective[1]

score <- matrix(0, n, m) # canonical scores <X_it, v_i>
wiizero <- (diag(w) == 0)
s <- ifelse(wiizero, rep(1,m), diag(w)) # scaling term 
if (sweep == "cyclical") idxi <- 1:m

## Check if some datasets are zero
eps <- 1e-14
xzero <- logical(m)
for (i in 1:m) {
	xzero[i] <- all(abs(range(x[[i]])) <= eps)
	if (xzero[i])
		v[[i]] <- lapply(p[[i]], 
			function(len) rep(1/sqrt(len), len))
}

## Block Coordinate Ascent 
for (it in 1:maxit) {
	if (sweep == "random") idxi <- sample(m)
	for (i in 1:m) { 		
		if (xzero[i]) next
		## Calculate the scores <X_jt, v_j> 
		## After the first algorithm iteration (it = 1), in each 
		## iteration of the i loop, only the inner products associated 
		## with the previous value of i need being updated
		idxj <- if (it == 1 && i == 1) { idxi[-1]  
			} else if (i == 1) { lastidx } else idxi[i-1]
		for (j in idxj) 
			score[, j] <- tnsr.vec.prod(x[[j]], v[[j]], 1:d[j]) 
		lastidx <- idxi[m]
		
		## Set up quadratic program 
		a <- if (wiizero[i]) 0 else x[[i]] # quadratic component
		b <- tnsr.vec.prod(x = x[[i]], 
			v = list(score[, -i] %*% (w[-i, i] / s[i])), 
			modes = d[i] + 1L) # linear component
		
		## Update canonical vectors
		v[[i]] <- optim.block.cov(v[[i]], a, b, maxit, tol) 	
	}								
	
	## Calculate objective value 
	objective[it + 1L] <- objective.internal(x, v, w)
	if (objective[it + 1L] > objective.best) {
		vbest <- v
		objective.best <- objective[it + 1L]
	}
	if (verbose) 
		cat("\nIteration", it, "Objective", objective[it + 1L])
	
	## Check convergence 
	if (it > 1 && abs(objective[it + 1] - objective[it]) <= 
	    	tol * max(1, objective[it])) break
}

list(v = vbest, y = canon.scores(x, vbest), 
	objective = objective.best, iters = it, 
	trace = objective[1:(it+1)])
}







#####################################################
# Maximize sum of covariances under global scaling
# and orthogonality constraints on canonical weights
#####################################################

		
# Note to self: for a given v(i,k) to optimize, 
# cannot rescale the other v(i,kk) to all have norm 1 
# because this will likely break orthogonality constraints :-(
 
	
mcca.cov.global.weight <- function(x, v, w, ortho, sweep, 
	maxit, tol, verbose)
{
	
## Data dimensions
dimx <- lapply(x, dim)
ndimx <- sapply(dimx, length)
d <- ndimx - 1 # number of image dimensions for each dataset
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
m <- length(x) # number of datasets
n <- dimx[[1]][ndimx[1]] # number of individuals/objects

objective <- numeric(maxit+1)
objective[1] <- objective.internal(x, v, w)
if (verbose) 
	cat("\nIteration", 0, "Objective", objective[1])

## Define the blocks of optimization variables 
ngroups <- prod(d)
## Dataset/mode combinations to include in each group
# group <- matrix(0, m, ngroups) 
if (sweep == "cyclical") {
	modes <- mapply(seq.int, 1, d, SIMPLIFY = FALSE)
	group <- expand.grid(modes, KEEP.OUT.ATTRS = FALSE)
	group <- matrix(unlist(group), m, ngroups, byrow = TRUE)
}

## Logical flags for datasets & associated variables to be 
## removed from optimization in case variables become zero
eps <- 1e-14 # numerical tolerance for zero
rmv <- logical(m)
nrm2v <- vector("list", m) 
nrmt.mk <- numeric(m) # products of norms of canonical vectors
# except for the one being optimized

## Test if the objective weights are equal or separable
equal.w <- all(w == w[1])
if (equal.w) {
	wvec <- rep(1/m, m)
	separable.w <- TRUE
} else {
	eigw <- eigen(w, TRUE)
	nonzero <- sum(abs(eigw$values) > (eps * eigw$values[1]))
	separable.w <- (nonzero == 1)
	wvec <- if (separable.w) { 
		eigw$vectors[,1] * sqrt(eigw$values[1]) 
	} else NULL	
}

## Orthogonality constraints
if (is.null(ortho)) {
	northo <- 0L 
} else {
	if (!is.matrix(ortho)) ortho <- as.matrix(ortho)
	northo <- ncol(ortho)
	cpfun <- function(x, y) sum(x * y)
}

for (it in 1:maxit) {
	## Determine blocks of variables to update in each dataset
	## in case of a random sweeping pattern
	if (sweep == "random") {
		group <- replicate(ngroups, sapply(d, sample.int, 
			size = 1, replace = TRUE))
	}	

	## Calculate squared norms of canonical vectors
	for (i in 1:m) 
		nrm2v[[i]] <- sapply(v[[i]], function(x) sum(x^2))

	for (g in 1:ngroups) {		

		## Set to zero canonical tensors that are approximately 
		## zero and remove them from the optimization
		for (i in 1:m) {
			if (any(nrm2v[[i]] < p[[i]] * eps)) { 
				v[[i]] <- lapply(p[[i]], numeric)
				rmv[i] <- TRUE
			} 
		}
		if (all(rmv)) break
		idxg <- which(!rmv)
		groupsize <- length(idxg)

		## Create block matrix of tensor-vector products
		## between X_it and vectors v_il in all modes 
		## l = 1, ..., d(i) except mode k(i,g) 
		## Each block corresponds to one dataset and the 
		## blocks are stacked vertically
		pg <- mapply("[", p[idxg], group[idxg, g])
		start <- cumsum(c(0, pg[-groupsize])) + 1
		end <- cumsum(pg)
		xv <- matrix(0, sum(pg), n)
		vortho <- matrix(0, sum(pg), northo)	
		for (ii in 1:groupsize) {
			i <- idxg[ii]
			k <- group[i,g] # selected mode for optimization
			tvprod <- tnsr.vec.prod(x[[i]], v[[i]][-k], 
				(1:d[i])[-k]) # tensor-vector products	
			nrmt.mk[i] <- sqrt(prod(nrm2v[[i]][-k]))	
			if (separable.w) 
				tvprod <- (wvec[i] / nrmt.mk[i] / sqrt(m*n)) * tvprod	
			xv[start[ii]:end[ii],] <- tvprod

			## Set up orthogonality constraints
			if (northo > 0) {
				for (l in 1:northo) {
					if (d[i] == 1L) {
						vortho[start[ii]:end[ii], l] <- 
						ortho[[i,l]][[k]]	
					} else {
					ip <- prod(mapply(cpfun, v[[i]][-k], 
						ortho[[i,l]][-k]))
					vortho[start[ii]:end[ii], l] <- 
						ortho[[i,l]][[k]] * (ip / nrmt.mk[i])
					}
				}
			}	
		}
		if (northo > 0) 	Q <- qr.Q(qr(vortho))
			
		## Update canonical vectors and objective value
		if (separable.w) { 
			# Case: objective weight matrix separable (SVD)
			if (northo > 0) # orthogonality constraints
				xv <- xv - Q %*% crossprod(Q, xv) 
			svdxv <- if (all(dim(xv) > 2)) {
				svds(xv, k = 1, nu = 1, nv = 0)
			} else { svd(xv, nu = 1, nv = 0) }
			vg <- svdxv$u
		} else {
			# Case: objective weight matrix nonseparable (EVD)
			xv <- tcrossprod(xv)
			s <- m * w / n / tcrossprod(nrmt.mk)
			vv <- vector("list", groupsize)
			for (ii in 1:groupsize) {
				i <- idxg[ii]
				vv[[i]] <- v[[i]][[group[i,g]]] * nrmt.mk[i] / sqrt(m)
				# debugging
				rows <- start[ii]:end[ii] 
				for (jj in 1:groupsize) {
					j <- idxg[jj]
					cols <- start[jj]:end[jj]
					xv[rows, cols] <- s[i,j] * xv[rows, cols]
				}
			}
			## Debugging
			# vv <- unlist(vv)
			# objective.implicit <- crossprod(vv, xv %*% vv)
			# objective.explicit <- objective.internal(x, v, w)
			# objective.test.1 <- (abs(objective.implicit - 
				# objective.explicit) < 1e-6)
			# print(paste("Test: objective (before ortho)", objective.test.1))
			
			if (northo > 0) { # project on orthogonality constraints
				QQtxv <- Q %*% crossprod(Q, xv)
				xv <- xv - QQtxv - t(QQtxv) + 
					tcrossprod(QQtxv %*% Q, Q)
			}				
			eigxv <- if (all(dim(xv) > 2)) {
				eigs_sym(xv, k = 1, which = "LA") 
			} else { eigen(xv, TRUE) }
			vg <- eigxv$vectors[,1]
		}
							
		## Update canonical weight vectors
		for (ii in 1:groupsize) {
			i <- idxg[ii]
			k <- group[i,g]
			idx <- start[ii]:end[ii]
			v[[i]][[k]] <- vg[idx] * (sqrt(m) / nrmt.mk[i])
			nrm2v[[i]][k] <- sum(v[[i]][[k]]^2)
		}
		
	}
									
	## Balance canonical weight vectors  
	v <- scale.v(v, type = "norm", scale = "global", 
		check.args = FALSE)
	objective[it+1] <- objective.internal(x, v, w)
	if (verbose) 
		cat("\nIteration", it, "Objective", objective[it+1])
	
	## Check convergence 
	if (it > 1 && abs(objective[it+1] - objective[it]) <= 
		max(1, tol) * abs(objective[it])) break
}

list(v = v, y = canon.scores(x, v), objective = objective[it+1], 
	iters = it, trace = objective[1:(it+1)])

}






#####################################################
# Maximize sum of covariances 
# under scaling constraints on canonical weights
# and orthogonality constraints on canonical scores
# All constraints are global (across datasets)
#####################################################
	
	
# Note to self: for a given v(i,k) to optimize, 
# cannot rescale the other v(i,kk) to all have norm 1 
# because this will likely break orthogonality constraints :-(
 
	
mcca.cov.global.score <- function(x, v, w, ortho, sweep, 
	maxit, tol, verbose)
{
	
## Data dimensions
dimx <- lapply(x, dim)
ndimx <- sapply(dimx, length)
d <- ndimx - 1 # number of image dimensions for each dataset
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
m <- length(x) # number of datasets
n <- dimx[[1]][ndimx[1]] # number of individuals/objects

objective <- numeric(maxit+1)
objective[1] <- objective.internal(x, v, w)
if (verbose) 
	cat("\nIteration", 0, "Objective", objective[1])

## Define the blocks of optimization variables 
ngroups <- prod(d)
## Dataset/mode combinations to include in each group
# group <- matrix(0, m, ngroups) 
if (sweep == "cyclical") {
	modes <- mapply(seq.int, 1, d, SIMPLIFY = FALSE)
	group <- expand.grid(modes, KEEP.OUT.ATTRS = FALSE)
	group <- matrix(unlist(group), m, ngroups, byrow = TRUE)
}

## Logical flags for datasets & associated variables to be 
## removed from optimization in case variables become zero
eps <- 1e-14 # numerical tolerance for zero
rmv <- logical(m)
nrm2v <- vector("list", m) 
nrmt.mk <- numeric(m) # products of norms of canonical vectors
# except for the one being optimized

## Test if the objective weights are equal or separable
equal.w <- all(w == w[1])
if (equal.w) {
	wvec <- rep(1/m, m)
	separable.w <- TRUE
} else {
	eigw <- eigen(w, TRUE)
	nonzero <- sum(abs(eigw$values) > (eps * eigw$values[1]))
	separable.w <- (nonzero == 1)
	wvec <- if (separable.w) { 
		eigw$vectors[,1] * sqrt(eigw$values[1]) 
	} else NULL	
}

## Orthogonality constraints
if (is.null(ortho)) {
	northo <- 0L 
} else {
	if (!is.matrix(ortho)) ortho <- as.matrix(ortho)
	northo <- ncol(ortho)
}

for (it in 1:maxit) {
	## Determine blocks of variables to update in each dataset
	## in case of a random sweeping pattern
	if (sweep == "random") {
		group <- replicate(ngroups, sapply(d, sample.int, 
			size = 1, replace = TRUE))
	}	

	## Calculate squared norms of canonical vectors
	for (i in 1:m) 
		nrm2v[[i]] <- sapply(v[[i]], function(x) sum(x^2))

	for (g in 1:ngroups) {		

		## Set to zero canonical tensors that are approximately 
		## zero and remove them from the optimization
		for (i in 1:m) {
			if (any(nrm2v[[i]] < p[[i]] * eps)) { 
				v[[i]] <- lapply(p[[i]], numeric)
				rmv[i] <- TRUE
			} 
		}
		if (all(rmv)) break
		idxg <- which(!rmv)
		groupsize <- length(idxg)

		## Create block matrix of tensor-vector products
		## between X_it and vectors v_il in all modes 
		## l = 1, ..., d(i) except mode k(i,g) 
		## Each block corresponds to one dataset and the 
		## blocks are stacked vertically
		pg <- mapply("[", p[idxg], group[idxg, g])
		start <- cumsum(c(0, pg[-groupsize])) + 1
		end <- cumsum(pg)
		xv <- matrix(0, sum(pg), n)
		vortho <- matrix(0, sum(pg), northo)	
		for (ii in 1:groupsize) {
			i <- idxg[ii]
			k <- group[i,g] # selected mode for optimization
			tvprod <- tnsr.vec.prod(x[[i]], v[[i]][-k], 
				(1:d[i])[-k]) # tensor-vector products	
			nrmt.mk[i] <- sqrt(prod(nrm2v[[i]][-k]))	
			if (separable.w) 
				tvprod <- (wvec[i] / nrmt.mk[i] / sqrt(m*n)) * tvprod	
			xv[start[ii]:end[ii],] <- tvprod

			## Set up orthogonality constraints
			if (northo > 0) {
				for (l in 1:northo) {
					if (d[i] == 1L) {
						vortho[start[ii]:end[ii], l] <- 
						ortho[[i,l]][[k]]	
					} else {
					ip <- prod(mapply(cpfun, v[[i]][-k], 
						ortho[[i,l]][-k]))
					vortho[start[ii]:end[ii], l] <- 
						ortho[[i,l]][[k]] * (ip / nrmt.mk[i])
					}
				}
			}	
		}
		if (northo > 0) Q <- qr.Q(qr(vortho))
			
		## Update canonical vectors and objective value
		if (separable.w) { 
			# Case: objective weight matrix separable (SVD)
			if (northo > 0) # orthogonality constraints
				xv <- xv - Q %*% crossprod(Q, xv) 
			svdxv <- if (all(dim(xv) > 2)) {
				svds(xv, k = 1, nu = 1, nv = 0)
			} else { svd(xv, nu = 1, nv = 0) }
			vg <- svdxv$u
		} else {
			# Case: objective weight matrix nonseparable (EVD)
			xv <- tcrossprod(xv)
			s <- m * w / n / tcrossprod(nrmt.mk)
			vv <- vector("list", groupsize)
			for (ii in 1:groupsize) {
				i <- idxg[ii]
				vv[[i]] <- v[[i]][[group[i,g]]] * nrmt.mk[i] / sqrt(m)
				# debugging
				rows <- start[ii]:end[ii] 
				for (jj in 1:groupsize) {
					j <- idxg[jj]
					cols <- start[jj]:end[jj]
					xv[rows, cols] <- s[i,j] * xv[rows, cols]
				}
			}
			## Debugging
			# vv <- unlist(vv)
			# objective.implicit <- crossprod(vv, xv %*% vv)
			# objective.explicit <- objective.internal(x, v, w)
			# objective.test.1 <- (abs(objective.implicit - 
				# objective.explicit) < 1e-6)
			# print(paste("Test: objective (before ortho)", objective.test.1))
			
			if (northo > 0) { # project on orthogonality constraints
				QQtxv <- Q %*% crossprod(Q, xv)
				xv <- xv - QQtxv - t(QQtxv) + 
					tcrossprod(QQtxv %*% Q, Q)
			}				
			eigxv <- if (all(dim(xv) > 2)) {
				eigs_sym(xv, k = 1, which = "LA") 
			} else { eigen(xv, TRUE) }
			vg <- eigxv$vectors[,1]
		}
							
		## Update canonical weight vectors
		for (ii in 1:groupsize) {
			i <- idxg[ii]
			k <- group[i,g]
			idx <- start[ii]:end[ii]
			v[[i]][[k]] <- vg[idx] * (sqrt(m) / nrmt.mk[i])
			nrm2v[[i]][k] <- sum(v[[i]][[k]]^2)
		}
		
	}
									
	## Balance canonical weight vectors  
	v <- scale.v(v, type = "norm", scale = "global", 
		check.args = FALSE)
	objective[it+1] <- objective.internal(x, v, w)
	if (verbose) 
		cat("\nIteration", it, "Objective", objective[it+1])
	
	## Check convergence 
	if (it > 1 && abs(objective[it+1] - objective[it]) <= 
		max(1, tol) * abs(objective[it])) break
}

list(v = v, y = canon.scores(x, v), objective = objective[it+1], 
	iters = it, trace = objective[1:(it+1)])

}



