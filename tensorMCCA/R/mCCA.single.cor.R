

###################################################
# Internal function for mCCA with mixed 2D/3D data
# (maximize sum of correlations one block at a time)
###################################################

# This function is not meant to be called by the user
# It calculates a single set of canonical vectors
# The optimization is conducted one data block at a time
	
mCCA.single.cor <- function(x, v, c, maxit = 1000, tol = 1e-6, 
	balance, verbose)
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

w <- matrix(0, m, n)
for (it in 1:maxit) {
	for (i in 1:m) { 		
		## Calculate/update objective weights w_jt = <X_jt, v_j> 
		## After the first algorithm iteration (it = 1), in each 
		## iteration of the i loop, only w_{i-1,t}, t=1:n, 
		## (resp. w_mt, t=1:n) need updating if i > 1 (resp. i = 1)
		idxj <- if (it == 1 && i == 1) {
			2:m } else if (i == 1) { m 
			} else (i-1)
		for (j in idxj) {
			tvprod <- x[[j]] # tvprod stands for tensor-vector product
			if (ndim.img[j] == 1) { # 1D case
				w[j,] <- c[i,j] * crossprod(v[[j]][[1]], tvprod) 
			} else if (ndim.img[j] == 2) { # 2D case
				dim(tvprod) <- c(dimx[[j]][1], prod(dimx[[j]][-1]))
				tvprod <- crossprod(v[[j]][[1]], tvprod)
				dim(tvprod) <- c(dimx[[j]][2], n) 
				w[j,] <- c[i,j] * crossprod(v[[j]][[2]], tvprod) 
			} else { # 3D case
				dim(tvprod) <- c(dimx[[j]][1], prod(dimx[[j]][-1]))
				tvprod <- crossprod(v[[j]][[1]], tvprod)
				dim(tvprod) <- c(dimx[[j]][2], dimx[[j]][3] * n)
				tvprod <- crossprod(v[[j]][[2]], tvprod)
				dim(tvprod) <- c(dimx[[j]][3], n) 
				w[j,] <- c[i,j] * crossprod(v[[j]][[3]], tvprod)
			}
		}
		colsumw <- colSums(w[-i,, drop = FALSE]) / n
		tvprod <- x[[i]]		
		dim(tvprod) <- c(prod(dimx[[i]][-ndimx[i]]), n)
		tvprod <- tvprod %*% colsumw
		dim(tvprod) <- dimx[[i]][-ndimx[i]]
		
		## Update canonical vectors
		v[[i]] <- optim.cor(v[[i]], tvprod, x[[i]], maxit, tol)				
	}								
	
	## Calculate objective value (sum of correlations)
	objective[it+1] <- objective.internal(x, v, c)	

	if (verbose) 
		cat("\nIteration",it,"Objective",objective[it+1])

	## Optionally balance the canonical vectors of each dataset
	## so they have equal norm
	if (balance) {
		for (i in 1:m) {
			nrm <- sapply(v[[i]], function(x) sqrt(sum(x^2)))
			if (any(nrm < 1e-15)) next
			d <- length(nrm)
			s <- prod(nrm)^(1/d) / nrm
			for (k in 1:d)
				v[[i]][[k]] <- s[k] * v[[i]][[k]]			
		}
	}
	
	## Check convergence 
	if (it > 1 && abs(objective[it+1]-objective[it]) <= 
	    	tol * max(1,objective[it])) break
}

list(v = v, y = image.scores(x, v), objective = objective[it+1], 
	iters = it, trace = objective[1:(it+1)])

}



