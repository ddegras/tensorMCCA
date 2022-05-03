

###############################
# Internal function for MCCA 
# Maximize sum of covariances 
# under block norm constraints 
###############################

# This function is not meant to be called by the user
# It calculates a single set of canonical vectors
# The optimization is conducted one data block at a time
	
mCCA.single.block.cov <- function(x, v, c, sweep, maxit, 
	tol, verbose)
{
	
## Data dimensions
dimx <- lapply(x, dim)
ndimx <- sapply(dimx, length)
d <- ndimx - 1
m <- length(x)
n <- dimx[[1]][ndimx[1]]

objective <- numeric(maxit+1)
objective[1] <- objective.internal(x, v, c)
if (verbose) 
	cat("\nIteration",0,"Objective",objective[1])

iprod <- matrix(0, m, n) # inner tensor products <X_it, v_i>
ciizero <- diag(c) == 0
s <- ifelse(ciizero, rep(1,m), diag(c)) # scaling term 
if (sweep == "cyclical") idxi <- 1:m

for (it in 1:maxit) {
	if (sweep == "random") idxi <- sample(m)
	for (i in 1:m) { 		
		## Calculate the inner products <X_jt, v_j> 
		## After the first algorithm iteration (it = 1), in each 
		## iteration of the i loop, only the inner products associated 
		## with the previous value of i need being updated
		idxj <- if (it == 1 && i == 1) { idxi[-1]  
			} else if (i == 1) { lastidx } else idxi[i-1]
		for (j in idxj) {
			## Calculate tensor-vector products between X_j and 
			## all vectors v_(jk) for k != i
			tvprod <- x[[j]] 
			pj <- head(dimx[[j]], d[j])
			vj <- v[[j]]
			for (k in 1:d[j]) {
				if (k < d[j]) {
					dim(tvprod) <- c(pj[k], prod(pj[(k+1):d[j]], n))
					tvprod <- crossprod(vj[[k]], tvprod)
					dim(tvprod) <- c(pj[k+1], length(tvprod) / pj[k+1])	
				} else {
					iprod[j,] <- crossprod(vj[[k]], tvprod)
				}				
			}
			# if (d[j] == 1) { # 1D case
				# iprod[j,] <- crossprod(vj[[1]], tvprod) 
			# } else if (d[j] == 2) { # 2D case
				# dim(tvprod) <- c(pj[1], prod(pj[-1], n))
				# tvprod <- crossprod(vj[[1]], tvprod)
				# dim(tvprod) <- c(pj[2], n) 
				# iprod[j,] <- crossprod(vj[[2]], tvprod) 
			# } else { # 3D case
				# dim(tvprod) <- c(pj[1], prod(pj[-1], n))
				# tvprod <- crossprod(vj[[1]], tvprod)
				# dim(tvprod) <- c(pj[2], pj[3] * n)
				# tvprod <- crossprod(vj[[2]], tvprod)
				# dim(tvprod) <- c(pj[3], n) 
				# iprod[j,] <- crossprod(vj[[3]], tvprod)
			# }
		}
		lastidx <- idxi[m]
		
		## Set up quadratic program 
		a <- if (ciizero[i]) 0 else x[[i]] # quadratic component
		w <- crossprod(iprod[-i,, drop = FALSE], c[-i, i] / s[i])
		b <- x[[i]]	# linear component
		dim(b) <- c(prod(dimx[[i]][-ndimx[i]]), n)
		b <- b %*% w
		if (d[i] > 1) dim(b) <- dimx[[i]][-ndimx[i]]
		
		## Update canonical vectors
		v[[i]] <- optim.block.cov(v[[i]], a, b, maxit, tol) 	
					
	}								
	
	## Calculate objective value (sum of correlations)
	objective[it+1] <- objective.internal(x, v, c)	
	if (verbose) 
		cat("\nIteration",it,"Objective",objective[it+1])
	
	## Check convergence 
	if (it > 1 && abs(objective[it+1]-objective[it]) <= 
	    	tol * max(1,objective[it])) break
}

list(v = v, y = image.scores(x, v), objective = objective[it+1], 
	iters = it, trace = objective[1:(it+1)])

}







################################
# Internal function for MCCA 
# Maximize sum of covariances 
# under global norm constraints 
################################

# This function is not meant to be called by the user
# It calculates a single set of canonical vectors
# The optimization is conducted one data block at a time
	
mCCA.single.global.cov <- function(x, v, c, sweep, 
	maxit, tol, verbose)
{
	
## Data dimensions
dimx <- lapply(x, dim)
ndimx <- sapply(dimx, length)
d <- ndimx - 1 # number of image dimensions for each dataset
m <- length(x) # number of datasets
n <- dimx[[1]][ndimx[1]] # number of individuals/objects

objective <- numeric(maxit+1)
objective[1] <- objective.internal(x, v, c)
if (verbose) 
	cat("\nIteration",0,"Objective",objective[1])

## Define the blocks of optimization variables 
ngroups <- max(d)
## Dataset/mode combinations to include in each group
group <- matrix(0, m, ngroups) 
if (sweep == "cyclical") {
	for (i in 1:m) 
		group[i,1:d[i]] <- 1:d[i]
}

## Logical flags for datasets & associated variables to be 
## removed from optimization in case variables become zero
eps <- 1e-14 # numerical tolerance for zero
rmv <- logical(m)

## Test if the objective weights are equal or separable
cequal <- all(c == c[1])
if (cequal) csep <- TRUE
if (!cequal) {
	eigc <- eigen(c, TRUE)
	csep <- (sum(eigc$values > eps) == 1)
	cvec <- if (csep) { 
		eigc$vectors[,1] * sqrt(eigc$values[1]) } else NULL
}

for (it in 1:maxit) {
	## Determine blocks of variables to update in each dataset
	## in case of a random sweeping pattern
	if (sweep == "random") {
		group <- matrix(0, m, ngroups)
		for (i in 1:m)
			group[i,sample.int(ngroups, d[i])] <- 1:d[i] 
	}
	
	for (g in 1:ngroups) {
		## Calculate norms of canonical vectors and rescale 
		## vectors that are not part of the current update 
		## so their product is one
		nrmt <- numeric(m)
		for (i in 1:m) {
			nrmv <- sapply(v[[i]], function(x) sqrt(sum(x^2)))
			nrmt[i] <- prod(nrmv)
			if (nrmt[i] < eps) { 
				nrmt[i] <- 0
				v[[i]] <- lapply(dimx[[i]][-d[i]], numeric)
				rmv[i] <- TRUE
			} else if (group[i,g] > 0) {
				for (k in 1:d[i]) {
					if (k == group[i,g]) next
					v[[i]][[k]] <- v[[i]][[k]] / nrmv[k] 
				}
			}
		}

		## Create block matrix of tensor-vector products
		## between X_it and vectors v_il in all modes 
		## l = 1, ..., d(i) except mode k(i,g) 
		## Each block corresponds to one dataset and the 
		## blocks are stacked vertically
		idxg <- which(group[,g] > 0 & (!rmv))
		if (length(idxg) == 0) next
		p <- mapply("[", dimx[idxg], group[idxg,g])
		groupsize <- length(idxg)
		start <- cumsum(c(0, p[-groupsize])) + 1
		end <- cumsum(p)
		xv <- matrix(0, sum(p), n)		
		for (ii in 1:groupsize) {
			i <- idxg[ii]
			tvprod <- x[[i]] 
			pi <- head(dimx[[i]], d[i])
			vi <- v[[i]]
			kopt <- group[i,g] # selected mode for optimization
			## If the selected mode is not next to last in X_i
			## permute the dimensions of X_i to enforce it
			if (kopt < d[i]) {
				perm <- c((1:d[i])[-kopt], kopt, d[i]+1)
				tvprod <- aperm(tvprod, perm)
			}
			## Calculate successive tensor-vector products
			for (k in 1:d[i]) {
				if (k == kopt) next
				dim(tvprod) <- c(pi[k], length(tvprod)/pi[k])
				tvprod <- crossprod(vi[[k]], tvprod)
			}
			dim(tvprod) <- c(pi[kopt], n)
			if (csep && (!cequal)) 
				tvprod <- cvec[i] * tvprod
			
			## Fill matrix block			
			xv[start[ii]:end[ii],] <- tvprod
		}
				
		## Update canonical vectors and objective value
		if (csep) { 
			# Case: objective weight matrix separable (SVD)
			svdxv <- if (max(dim(xv)) > 2) {
				svds(xv, k = 1, nu = 1, nv = 0)
			} else { svd(xv, nu = 1, nv = 1) }
			vg <- svdxv$u
			if (ii == groupsize && all(group[,g] > 0)) {
				objective[it+1] <- ifelse(cequal,
					(svdxv$d)^2 * c[1] * m / n,
					(svdxv$d)^2 * m / n) }
		} else {
			# Case: objective weight matrix nonseparable (EVD)
			xv <- tcrossprod(xv)
			for (ii in 1:groupsize) {
				i <- idxg[ii]
				rows <- start[ii]:end[ii] 
				for (jj in 1:groupsize) {
					j <- idxg[jj]
					cols <- start[ij]:end[jj]
					xv[rows, cols] <- c[i, j] * xv[rows, cols]
				}			
			}
			eigxv <- if (max(dim(xv)) > 2) {
				eigs_sym(xv, k = 1) } else { eigen(xv, TRUE) }
			vg <- eigxv$vectors[,1]
			if (ii == groupsize && all(group[,g] > 0)) {
				objective[it+1] <- ifelse(cequal,
					eigxv$values[1] * c[1] * m / n,
					eigxv$values[1] * m / n) }
		}
		## Rescale solution to meet global norm constraint 
		s <- sqrt(m - sum(nrmt[group[,g] == 0]^2))
		vg <-  s * vg
		for (ii in 1:groupsize) {
			i <- idxg[ii]
			k <- group[i,g]
			v[[i]][[k]] <- vg[start[ii]:end[ii]]
		}
	
					
	}								
	
	## Balance canonical vectors  
	v <- scale.v(v, type = "norm", cnstr = "global")
	
	## Debugging: Compare objective value calculated in loop 
	## to full calculation
	# test <- objective.internal(x, v, c)
	if (any(group[,ngroups] == 0))
		objective[it+1] <- objective.internal(x, v, c)
	if (verbose) 
		cat("\nIteration",it,"Objective",objective[it+1])
	# browser()
	
	## Check convergence 
	if (it > 1 && abs(objective[it+1]-objective[it]) <= 
	    	tol * max(1,objective[it])) break
}

list(v = v, y = image.scores(x, v), objective = objective[it+1], 
	iters = it, trace = objective[1:(it+1)])

}



