## INTERNAL FUNCTIONS

####################################################
# Maximize sum of covariances 
# under scaling constraints on canonical weights
# and orthogonality constraints on canonical scores 
# All constraints are at the block level
####################################################

	
mcca.cov.bca.block <- function(x, v, w, ortho, sweep, 
	maxit, tol, verbose)
{
	
## Determine data dimensions
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
m <- length(x)
n <- tail(dimx[[1]], 1)

## Set up objective values
objective <- numeric(maxit + 1L)
objective[1] <- objective.internal(x, v, w)
if (verbose) 
	cat("\nIteration", 0, "Objective", objective[1])
vbest <- v
objective.best <- objective[1]

score <- matrix(0, n, m) # canonical scores <X_it, v_i>
s <- diag(w) # scaling term 
s[s == 0] <- 1
if (sweep == "cyclical") idxi <- 1:m
if (!is.null(ortho) && !is.matrix(ortho)) 
	dim(ortho) <- c(m, 1)

## Check if some datasets are zero
eps <- 1e-14
xzero <- logical(m)
for (i in 1:m) {
	xzero[i] <- all(abs(range(x[[i]])) <= eps)
	if (xzero[i])
		v[[i]] <- lapply(p[[i]], 
			function(len) rep(1/sqrt(len), len))
}
if (all(xzero))
	return(list(v = v, score = score, objective = 0, iters = 1, 
		trace = 0))
	
## Block coordinate ascent 
for (it in 1:maxit) {
	if (sweep == "random") idxi <- sample(m)
	for (ii in 1:m) { 	
		i <- idxi[ii]	
		if (xzero[i]) next
		## Calculate the scores <X_jt, v_j> 
		## After the first algorithm iteration (it = 1), in each 
		## iteration of the i loop, only the inner products associated 
		## with the previous value of i need being updated
		idxj <- if (it == 1 && ii == 1) { idxi[-1]  
			} else if (ii == 1) { lastidx } else idxi[ii-1]
		for (j in idxj) 
			score[, j] <- tnsr.vec.prod(x[[j]], v[[j]], 1:d[j]) 
		
		## Set up quadratic program 
		a <- if (w[i,i] == 0) 0 else x[[i]] # quadratic component
		b <- tnsr.vec.prod(x[[i]], 
			if (m == 2) { list(score[, -i] * (w[-i, i] / s[i]))
			} else list(score[, -i] %*% (w[-i, i] / s[i])), 
			d[i] + 1L) # linear component
		
		## Update canonical vectors
		v[[i]] <- optim.block.cov(v[[i]], a, b, ortho[i,], maxit, tol)
	}								
	lastidx <- idxi[m]
	
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

list(v = vbest, score = canon.scores(x, vbest), 
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
 
	
mcca.cov.bca.global <- function(x, v, w, ortho, sweep, 
	maxit, tol, verbose)
{
	
## Determine data dimensions
dimx <- lapply(x, dim)
ndimx <- sapply(dimx, length)
d <- ndimx - 1 # number of image dimensions for each dataset
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
m <- length(x) # number of datasets
n <- dimx[[1]][ndimx[1]] # number of individuals/objects

## Set up objective values
objective <- numeric(maxit+1)
objective[1] <- objective.internal(x, v, w)
if (verbose && is.null(ortho)) 
	cat("\nIteration", 0, "Objective", objective[1])
# Starting value typically not feasible so don't show its objective value

## Catch trivial case
if (all(unlist(x) == 0)) {
	p <- unlist(p)
	return(list(v = relist(rep(1/sqrt(p), p), v), 
		score = matrix(0, n, m), objective = 0, iters = 1, trace = 0))
}

## Identify type of orthogonality constraints
ortho.type <- if (is.numeric(ortho[[1]])) { "score" 
	} else if (is.list(ortho[[1]])) { "weight"
	} else { "none" }

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

## Test if the objective weights are separable
## Recall: weights are nonnegative and add up to 1
csumw <- colSums(w)
separable.w <- isTRUE(all.equal(w, tcrossprod(csumw)))

## Ensure that orthogonality constraints are in matrix format
if (!is.null(ortho)) {
	if (!is.matrix(ortho)) dim(ortho) <- c(m,1)
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
		if (!is.null(ortho)) 
			vortho <- matrix(0, sum(pg), northo)	
		for (ii in 1:groupsize) {
			i <- idxg[ii]
			k <- group[i,g] # selected mode for optimization
			tvprod <- tnsr.vec.prod(x[[i]], v[[i]][-k], 
				(1:d[i])[-k]) # tensor-vector products	
			nrmt.mk[i] <- sqrt(prod(nrm2v[[i]][-k]))	
			if (separable.w) 
				tvprod <- (csumw[i] / nrmt.mk[i] / sqrt(m*n)) * tvprod
			rows <- start[ii]:end[ii]	
			xv[rows,] <- tvprod

			## Set up orthogonality constraints
			if (ortho.type == "weight") {			
				for (l in 1:northo) {
					if (d[i] == 1L) {
						vortho[rows, l] <- ortho[[i,l]][[k]]	
					} else {
					iprod <- prod(mapply(cpfun, v[[i]][-k], 
						ortho[[i,l]][-k]))
					vortho[rows, l] <- ortho[[i,l]][[k]] * 
						(iprod / nrmt.mk[i])
					}
				}
			} else if (ortho.type == "score") {
				for (l in 1:northo) {
					if (d[i] == 1L) {
						vortho[rows, l] <- ortho[[i,l]]	
					} else {
					tvprod <- tnsr.vec.prod(ortho[[i,l]], 
						v[[i]][-k], (1:d[i])[-k])
					vortho[rows, l] <- tvprod / nrmt.mk[i]
					}
				}				
			}	
		}
		if (!is.null(ortho)) Q <- qr.Q(qr(vortho))
			
		## Update canonical vectors and objective value
		if (separable.w) { 
			# Case: objective weight matrix separable (SVD)
			if (!is.null(ortho)) # orthogonality constraints
				xv <- xv - Q %*% crossprod(Q, xv) 
			svdxv <- if (all(dim(xv) > 2)) {
				svds(xv, k = 1, nu = 1, nv = 0)
			} else { svd(xv, nu = 1, nv = 0) }
			vg <- svdxv$u
			if (all(is.nan(vg))) 
				vg <- rep(1/sqrt(length(vg)), length(vg))
		} else {
			# Case: objective weight matrix nonseparable (EVD)
			xv <- tcrossprod(xv)
			s <- m * w / n / tcrossprod(nrmt.mk)
			# vv <- vector("list", groupsize)
			for (ii in 1:groupsize) {
				i <- idxg[ii]
				# vv[[i]] <- v[[i]][[group[i,g]]] * nrmt.mk[i] / sqrt(m)
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
			
			if (!is.null(ortho)) { # project on orthogonality constraints
				QQtxv <- Q %*% crossprod(Q, xv)
				xv <- xv - QQtxv - t(QQtxv) + 
					tcrossprod(QQtxv %*% Q, Q)
			}				
			eigxv <- if (all(dim(xv) > 2)) {
				eigs_sym(xv, k = 1, which = "LA") 
			} else { eigen(xv, TRUE) }
			vg <- eigxv$vectors[,1]
			if (all(is.nan(vg)))
				vg <- rep(1/sqrt(length(vg)), length(vg))
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

list(v = v, score = canon.scores(x, v), objective = objective[it+1], 
	iters = it, trace = objective[1:(it+1)])

}




