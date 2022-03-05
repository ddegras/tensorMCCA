

###################################################
# Internal function for MCCA 
# Maximize sum of covariances one block at a time
###################################################

# This function is not meant to be called by the user
# It calculates a single set of canonical vectors
# The optimization is conducted one data block at a time
	
mCCA.single.cov <- function(x, v, c, maxit = 1000, tol = 1e-6, verbose)
{
	
## Data dimensions
dimx <- lapply(x, dim)
ndimx <- sapply(dimx, length)
ndim.img <- ndimx - 1
m <- length(x)
n <- dimx[[1]][ndimx[1]]

objective <- numeric(maxit+1)
objective[1] <- objective.internal(x, v, c)
if (verbose) 
	cat("\nIteration",0,"Objective",objective[1])

iprod <- matrix(0, m, n) # inner tensor products <X_it, v_i>
ciizero <- diag(c) == 0
s <- ifelse(ciizero, rep(1,m), diag(c)) # scaling term 

for (it in 1:maxit) {
	for (i in 1:m) { 		
		## Calculate/update the inner products <X_jt, v_j> 
		## After the first algorithm iteration (it = 1), in each 
		## iteration of the i loop, only the <X_(i-1)t, v_(i-1)>, t=1:n, 
		## (resp. the <X_mt, v_m), t=1:n) need updating if i > 1 
		## (resp. if i = 1)
		idxj <- if (it == 1 && i == 1) { 2:m  
			} else if (i == 1) { m } else (i-1)
		for (j in idxj) {
			## Calculate tensor-vector products between X_j and 
			## all vectors v_(jk) for k != i
			tvprod <- x[[j]] 
			pj <- head(dimx[[j]], ndim.img[j])
			vj <- v[[j]]
			if (ndim.img[j] == 1) { # 1D case
				iprod[j,] <- crossprod(vj[[1]], tvprod) 
			} else if (ndim.img[j] == 2) { # 2D case
				dim(tvprod) <- c(pj[1], prod(pj[-1], n))
				tvprod <- crossprod(vj[[1]], tvprod)
				dim(tvprod) <- c(pj[2], n) 
				iprod[j,] <- crossprod(vj[[2]], tvprod) 
			} else { # 3D case
				dim(tvprod) <- c(pj[1], prod(pj[-1], n))
				tvprod <- crossprod(vj[[1]], tvprod)
				dim(tvprod) <- c(pj[2], pj[3] * n)
				tvprod <- crossprod(vj[[2]], tvprod)
				dim(tvprod) <- c(pj[3], n) 
				iprod[j,] <- crossprod(vj[[3]], tvprod)
			}
		}
		
		## Set up quadratic and linear components of objective to maximize 
		a <- if (ciizero[i]) 0 else x[[i]] # quadratic component
		w <- crossprod(iprod[-i,, drop = FALSE], c[-i, i] / s[i])
		b <- x[[i]]	# linear component
		dim(b) <- c(prod(dimx[[i]][-ndimx[i]]), n)
		b <- b %*% w
		if (ndim.img[i] > 1) dim(b) <- dimx[[i]][-ndimx[i]]
		
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



