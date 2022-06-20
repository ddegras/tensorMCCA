

################################################
# Calculate gradient of mCCA objective function
################################################

	
mCCA.gradient <- function(x, v, c)
{
	
## Data dimensions
dimx <- lapply(x, dim)
lenx <- sapply(dimx, prod)
ndimx <- sapply(dimx, length)
d <- ndimx - 1
m <- length(x)
n <- dimx[[1]][ndimx[1]]

## Calculate tensor-vector products X_it x v_i^(-k) 
## and scores  X_it, v_i > 
tvprod <- vector("list", m)
score <- matrix(0, n, m)
for (i in 1:m) { 	
	xi <- x[[i]]
	pi <- p[[i]]
	for (k in 1:d[i]) {
		sigma <- c(1:ndimx[i])[-k], k)
		xi <- aperm(x[[i]], sigma)
		dim(xi) <- c(lenx[i] / n / pi[k], n * pi[k])
		vi <- Reduce(rev(v[[i]][-k]), kronecker)
		vtxi <- crossprod(vi, xi)
		dim(vtxi) <- c(n, pi[k])
		tvprod[[i]][[k]] <- vtxi
		}
	score[,i] <- vtxi %*% v[[i]][k]
}

## Calculate gradient
for (i in 1:m) 
	grad[[i]] <- lapply(tvprod[[i]], crossprod, 
		y = score %*% ((2/n) * c[,i]))

grad
		
}



