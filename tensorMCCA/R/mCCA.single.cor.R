

###################################################
# Internal function for mCCA with mixed 2D/3D data
# (maximize sum of correlations one block at a time)
###################################################

# This function is not meant to be called by the user
# It calculates a single set of canonical vectors
# The optimization is conducted one data block at a time
	
mCCA.single.cor <- function(x, v, c, sweep, maxit, tol, balance, verbose)
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

iprod <- matrix(0, m, n)
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
					tvprod <- drop(crossprod(vj[[k]], tvprod))
					dim(tvprod) <- c(pj[k+1], length(tvprod) / pj[k+1])	
				} else {
					iprod[j,] <- crossprod(vj[[k]], tvprod)
				}				
			}
		}
		lastidx <- idxi[m]
		
		## Set up linear program
		w <- crossprod(iprod[-i,, drop = FALSE], c[-i, i] / n)
		if (d[i] == 1) {
			a <- x[[i]] %*% w	
		} else {
			a <- x[[i]]
			dim(a) <- c(prod(dimx[[i]][-ndimx[i]]), n)
			a <- a %*% w
			dim(a) <- dimx[[i]][-ndimx[i]]
		}
		
		## Update canonical vectors
		v[[i]] <- optim.block.cor(v[[i]], a, x[[i]], maxit, tol)				
	}								
	
	## Calculate objective value (sum of correlations)
	objective[it+1] <- objective.internal(x, v, c)	

	if (verbose) 
		cat("\nIteration",it,"Objective",objective[it+1])

	## Optionally balance the canonical vectors of each dataset
	## so they have equal norm
	if (balance) {
		for (i in 1:m) {
			nrmv <- sapply(v[[i]], function(x) sqrt(sum(x^2)))
			nrmt <- prod(nrmv)
			s <- if (nrmt < 1e-15) { 
				numeric(d[i]) } else nrmt^(1/d[i]) / nrmv
			for (k in 1:d[i])
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



