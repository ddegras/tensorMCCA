

############################################
# Initialize MCCA by solving pairwise CCA's
############################################



init.mcca.cca <- function(x, k = NULL, c = 1, objective = c("cov", "cor"), 
	cnstr = c("block", "global"), center = TRUE)
{
## Check argument x if required
test <- check.arguments(x)
tol <- 1e-14

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
pp <- sapply(p, prod)
d <- sapply(p, length)
n <- sapply(dimx, tail, 1)
if (!all(n == n[1])) 
	stop(paste("The last modes of the components of 'x'",
		"must have the same dimension"))
n <- n[[1]]
if (is.null(k)) k <- n else k <- min(k, n) 
	
## Scaling constraints
objective.type <- match.arg(objective)   # norm or variance constraints
# type <- switch(objective.type, cov = "norm", cor = "var")
if (objective.type == "cor") cnstr <- "block"
cnstr <- match.arg(cnstr) # block or global constraints

## Objective weights
stopifnot(length(c) == 1 || (is.matrix(c) && all(dim(c) == length(x))))
stopifnot(all(c >= 0) && any(c > 0))
if (length(c) == 1) {
	c <- matrix(1/m^2, m, m)
} else {
	c <- (c + t(c)) / (2 * sum(c))
}

## Data centering
if (center) {
	xbar <- vector("list", m)
	for(i in 1:m) 
	    xbar[[i]] <- rowMeans(x[[i]], dims = d[i])
}

## Find the ranks of unfolded datasets
rankx <- numeric(m)
for (i in 1:m) {
	rankx[i] <- if (center) {
		qr(matrix(x[[i]], pp[i], n) - xbar[[i]])$rank	
	} else { qr(matrix(x[[i]], pp[i], n) - xbar[[i]])$rank }
}

## Effective rank of SVD 	
k <- pmin(rankx, k)

## Calculate SVDs of each unfolded dataset
reducex <- switch(objective.type, 
	cov = (k <= pp / 2), cor = rep(TRUE, m))
u <- dd <- xred <- vector("list", m) 
for (i in which(reducex)) {
	xmat <- if (center) { matrix(x[[i]] - xbar[[i]], pp[i], n)  
		} else { matrix(x[[i]], pp[i], n) }
	svdx <- if(max(n, pp[i]) > 2) {
		svds(xmat, k[i])
	} else {
		svd(xmat, nu = k[i], nv = k[i])
	}
	
	## Reduce data
	u[[i]] <- svdx$u  # left singular vectors
	dd[[i]] <- svdx$d # singular values
	xred[[i]] <- if (objective.type == "cov") { # reduced data
		t(svdx$v) * svdx$d }  else t(svdx$v) 
}
if (any(reducex)) rm(xmat, svdx)


## CASE: mCCA with global constraints
if (objective.type == "cov" && cnstr == "global") { 
	## Calculate rank of objective weight matrix
	const <- all(c == c[1])
	if (const) { 
		rankc <- 1
	 } else {
		eigc <- eigen(c, TRUE)
		rankc <- sum(eigc$values > tol)
		if (rankc == 1) 
			s <- sqrt(eigc$values[1]) * eigv$vectors[,1] # scaling factors
	}
	
	## Concatenate data
	len <- pp
	len[reducex] <- k[reducex]
	cumlen <- c(0, cumsum(len))
	block <- lapply(1:m, function(i) (cumlen[i]+1):cumlen[i+1])
	xmat <- matrix(nrow = sum(pp), ncol = n)
	noscale <- (const || rankc > 1)
	for (i in 1:m) {
		xmat[block[[i]],] <- if (reducex[i] && noscale) { 
			xred[[i]] 
			} else if (reducex[i]) { 
				s[i] * xred[[i]] 
			} else if (center) { 
				matrix(x[[i]] - xbar[[i]], pp[[i]], n)
		} else matrix(x[[i]], pp[[i]], n)
	}

	## Perform SVD of concatenated data if objective weight matrix has rank 1
	## or EVD of associated covariance matrix otherwise
	if (rankc == 1) {
		v <- if (max(dim(xmat)) > 2) {
			svds(x, k = 1)$u
		} else { svd(x, nu = 1, nv = 0)$u }
	} else {
		xmat <- tcrossprod(xmat)
		for (i in 1:m)
		for (j in 1:m)
			xmat[block[[i]], block[[j]]] <- c[i,j] * x[block[[i]], block[[j]]]
		v <- if (ncol(xmat) > 2) {
			eigs_sym(x, k = 1)$vectors
		} else { eigen(x, TRUE)$vectors[,1] }
	}
	
	## Recover (long) block canonical vectors
	vcopy <- sqrt(m) * v
	v <- vector("list", m)
	for (i in 1:m) 
		v[[i]] <- if (reducex[i]) { u[[i]] %*% vcopy[block[[i]]]
			} else { vcopy[block[[i]]] }
	rm(vcopy)		

} else { 
## CASE: mCCA with block constraints

	## Calculate CCA for each pair of datasets 	 
	v <- lapply(pp, function(nr) matrix(nrow = nr, ncol = n)) # canonical vectors
	objective <- matrix(0, m, m)
	for (i in 1:m) {
		xi <- if (reducex[i]) { 
			xred[[i]]
		} else if (center) {
			matrix(x[[i]] - xbar[[i]], pp[i], n)
		} else {
			matrix(x[[i]], pp[i], n)
		}		
		for (j in 1:i) {
			xj <- if (i == j) { xi 
			} else if (reducex[j]) { 
				xred[[j]]
			} else if (center) {
				matrix(x[[j]] - xbar[[j]], pp[j], n)
			} else {
				matrix(x[[j]], pp[j], n)
			}
			RSpectra.flag <- (max(nrow(xi), nrow(xj)) > 2)
			svdij <- if (i == j && RSpectra.flag) {
				svds(xi, k = 1, nv = 0)	
			} else if (i == j) {
				svd(xi, nu = 1, nv = 0)			
			} else if (RSpectra.flag) {
				svds(tcrossprod(xi, xj), k = 1)
			} else {
				svd(tcrossprod(xi, xj), nu = 1, nv = 1)
			} 
			objective[i,j] <- c[i,j] * ifelse(i == j, 
				svdij$d[1]^2, svdij$d[1])			
			v[[i]][,j] <- if (objective.type == "cor") {
				u[[i]] %*% (svdij$u / dd[[i]])
			} else if (reducex[i]) {
				u[[i]] %*% svdij$u 
			} else {
				svdij$u		
			}				
			if (j < i) {
				objective[j,i] <- objective[i,j]
				v[[j]][, i] <- if (objective.type == "cor") {
					u[[j]] %*% (svdij$v / dd[[j]])
				} else if (reducex[i]) { 
					u[[j]] %*% svdij$v
				} else {
					svdij$v
				}
			}
		}
	}	
}
 

## Find best overall assignment in pairwise CCA		
## (Linear Sum Assignment Problem)
if (cnstr == "block") {
	idxj <- solve_LSAP(objective, maximum = TRUE)
	for (i in 1:m) 
		v[[i]] <- v[[i]][, idxj[i]]
}

## Approximate long canonical vectors by rank-1 tensors 
for (i in 1:m) {
	pi <- p[[i]]
	vi <- v[[i]]
	v[[i]] <- vector("list", d[i])
	for (k in 1:d[i]) {
		if (k < d[i]) {
			dim(vi) <- c(pi[k], length(vi) / pi[k])
			svdi <- if (min(dim(vi)) > 2) svds(vi, 1) else svd(vi, 1, 1)
			v[[i]][[k]] <- svdi$u
			vi <- svdi$v
		} else {
			v[[i]][[k]] <- vi 
		}
	}
}

## Scale the canonical tensors as needed
## (norm constraints are already enforced)
if (objective.type == "cor")
	v <- scale.v(v, x, type = "var")

v

}


