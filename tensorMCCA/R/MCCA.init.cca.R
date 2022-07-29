

############################################
# Initialize MCCA by solving pairwise CCA's
############################################



mcca.init.cca <- function(x, k = NULL, w = 1, 
	objective = c("cov", "cor"), norm = c("block", "global"), 
	center = TRUE, search = c("exhaustive", "approximate"))
{
	
################
# Preprocessing
################


## Check argument x
test <- check.arguments(x, w = w)
eps <- 1e-14

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (omit last dimension = instances/cases)
cnst.set <- sapply(x, function(xx) all(abs(xx - xx[1]) < eps))
if (all(cnst.set)) {
	v <- vector("list", m)
	for (i in 1:m) 
		v[[i]] <- lapply(p[[i]], numeric)
	return(v)
} else if (any(cnst.set)) {
	m0 <- m
	dimx0 <- dimx
	p0 <- p
	x <- x[!cnst.set]
	m <- sum(!cnst.set)
	p <- p[!cnst.set]
}
pp <- sapply(p, prod)
d <- sapply(p, length)
n <- tail(dimx[[1]], 1)
if (is.null(k)) { k <- n } else { k <- min(k, n) }



## Scaling constraints
objective.type <- match.arg(objective)   
if (objective.type == "cor" && norm == "global")
	stop(paste("Argument values 'objective = cor'", 
		"and 'norm = global' are incompatible."))
norm <- match.arg(norm) # block or global constraints

## Search method in optimization
search <- if (identical(search, 
	c("exhaustive", "approximate"))) {
	ifelse(m <= 5, "exhaustive", "approximate")
} else { match.arg(search) }

## Objective weights
if (is.matrix(w)) 
	w <- w[!cnst.set, !cnst.set]
stopifnot(length(w) == 1 || 
	(is.matrix(w) && all(dim(w) == length(x))))
stopifnot(all(w >= 0) && any(w > 0))
if (length(w) == 1) {
	w <- matrix(1/m^2, m, m)
} else {
	w <- (w + t(w)) / (2 * sum(w))
}

## Data centering
if (center) {
	xbar <- vector("list", m)
	for(i in 1:m) 
	    xbar[[i]] <- as.vector(rowMeans(x[[i]], dims = d[i]))
}



###########################
# SVD of unfolded datasets
###########################


## Find the ranks of unfolded datasets
rankx <- numeric(m)
for (i in 1:m) {
	rankx[i] <- if (center) {
		qr(matrix(x[[i]], pp[i], n) - xbar[[i]])$rank	
	} else { qr(matrix(x[[i]], pp[i], n))$rank }
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
	svdx <- if (min(n, pp[i]) > max(2, k[i])) {
		svds(xmat, k[i])
	} else {
		svd(xmat, nu = k[i], nv = k[i])
	}
	
	## Reduce data
	u[[i]] <- svdx$u  # left singular vectors
	dd[[i]] <- svdx$d # singular values
	xred[[i]] <- if (objective.type == "cov") { 
		# reduced data
		t(svdx$v) * svdx$d }  else t(svdx$v) 
}
if (any(reducex)) rm(xmat, svdx)



#####################################
# Case: mCCA with global constraints
#####################################


if (objective.type == "cov" && norm == "global") { 
	## Calculate rank of objective weight matrix
	const <- all(w == w[1])
	if (const) { 
		rankw <- 1
	 } else {
		eigw <- eigen(w, TRUE)
		rankw <- sum(eigw$values > eps)
		if (rankw == 1) # scaling factors
			s <- sqrt(eigw$values[1]) * eigv$vectors[,1] 			
	}	
	## Concatenate data
	len <- pp
	len[reducex] <- k[reducex]
	cumlen <- c(0, cumsum(len))
	block <- lapply(1:m, function(i) (cumlen[i]+1):cumlen[i+1])
	xmat <- matrix(nrow = sum(pp), ncol = n)
	noscale <- (const || rankw > 1)
	for (i in 1:m) {
		xmat[block[[i]],] <- if (reducex[i] && noscale) { 
			xred[[i]] 
			} else if (reducex[i]) { 
				s[i] * xred[[i]] 
			} else if (center) { 
				matrix(x[[i]] - xbar[[i]], pp[[i]], n)
		} else matrix(x[[i]], pp[[i]], n)
	}
	## Perform SVD of concatenated data if objective weight 
	## matrix has rank 1 or EVD of associated covariance 
	## matrix otherwise
	if (rankw == 1) {
		v <- if (max(dim(xmat)) > 2) {
			svds(x, k = 1)$u
		} else { svd(x, nu = 1, nv = 0)$u }
	} else {
		xmat <- tcrossprod(xmat)
		for (i in 1:m)
		for (j in 1:m)
			xmat[block[[i]], block[[j]]] <- 
				w[i,j] * x[block[[i]], block[[j]]]
		v <- if (ncol(xmat) > 2) {
			eigs_sym(xmat, k = 1)$vectors
		} else { eigen(xmat, TRUE)$vectors[,1] }
	}
	
	## Recover (long) block canonical vectors
	vcopy <- sqrt(m) * v
	v <- vector("list", m)
	for (i in 1:m) 
		v[[i]] <- if (reducex[i]) { u[[i]] %*% vcopy[block[[i]]]
			} else { vcopy[block[[i]]] }
	rm(vcopy)		

}



#################################### 
# Case: mCCA with block constraints
#################################### 


## Calculate CCA for each pair of datasets 	
# Note: for the problem of maximizing the sum of correlations, 
# the canonical vectors are expressed up to scale factors. 
# Only their directions matter at this stage, not their scale.   

if (norm == "block") {
	v <- lapply(pp, function(nr) matrix(0, nr, m)) 
	# canonical vectors
	for (i in 1:m) {
		xi <- if (reducex[i]) { 
			xred[[i]]
		} else if (center) {
			matrix(x[[i]] - xbar[[i]], pp[i], n)
		} else {
			matrix(x[[i]], pp[i], n)
		}		
		for (j in 1:i) {
			if (i == j && objective.type == "cor") {
				v[[i]][,i] <- u[[i]][,1] 
				next
			}
			xj <- if (i == j) { xi 
			} else if (reducex[j]) { 
				xred[[j]]
			} else if (center) {
				matrix(x[[j]] - xbar[[j]], pp[j], n)
			} else {
				matrix(x[[j]], pp[j], n)
			}
			RSpectra.flag <- (min(nrow(xi), nrow(xj)) > 2)
			svdij <- if (i == j && RSpectra.flag) {
				svds(xi, k = 1, nv = 0)	
			} else if (i == j) {
				svd(xi, nu = 1, nv = 0)			
			} else if (RSpectra.flag) {
				svds(tcrossprod(xi, xj), k = 1)
			} else {
				svd(tcrossprod(xi, xj), nu = 1, nv = 1)
			} 		
			v[[i]][,j] <- if (objective.type == "cor") {
				u[[i]] %*% (svdij$u / dd[[i]])
			} else if (reducex[i]) {
				u[[i]] %*% svdij$u 
			} else {
				svdij$u		
			}				
			if (i == j) next
			v[[j]][,i] <- if (objective.type == "cor") {
				u[[j]] %*% (svdij$v / dd[[j]])
			} else if (reducex[j]) { 
				u[[j]] %*% svdij$v
			} else {
				svdij$v
			}
		}
	}	
}



#####################################
# Approximate long canonical vectors 
# by rank-1 tensors 
#####################################

vcopy <- v
v <- vector("list", m^2)
dim(v) <- c(m, m)
for (i in 1:m) {
	pi <- p[[i]]
	test <- all(pi > 2)
	for (j in 1:m) {
		if (d[i] == 1) {
			vij <- vcopy[[i]][,j]
			v[[i,j]] <- list(vij / sqrt(sum(vij^2)))
		} else if (d[i] == 2) {
			vij <- matrix(vcopy[[i]][,j], pi[1], pi[2])
			svdij <- if (test) { svds(vij, 1) 
				} else svd(vij, 1, 1)
			v[[i,j]] <- list(svdij$u, svdij$v)
		} else {
			vij <- array(vcopy[[i]][,j], pi)
			v[[i,j]] <- tnsr.rk1(vij, scale = TRUE)
		}
	}
}
rm(vcopy)

# After this stage, all canonical vectors have unit norm



####################################
# Calculate image scores associated 
# with each canonical tensor
####################################


score <- canon.scores(x, v)
# score <- array(dim = c(n, m, m))
# for (i in 1:m)
# for (j in 1:m)
# {
	# vij <- if (d[i] == 1) { unlist(v[[i,j]]) 
		# } else { Reduce(kronecker, rev(v[[i,j]])) }
	# dim(vij) <- NULL
	# score[,i,j] <- colSums(vij * x[[i]], dims = d[i])
	# if (center) score[,i,j] <- score[,i,j] - sum(vij * xbar[[i]])
# }



####################################
# Rescale canonical vectors 
# if maximizing sum of correlations
####################################


if (objective.type == "cor") {
	for (i in 1:m) {
		s <- sqrt(colMeans(score[,i,]^2))
		s[s <= eps] <- Inf
		score[,i,] <- sweep(score[,i,], 2, s, "/")
		for (j in 1:m) 
			v[[i,j]] <- lapply(v[[i,j]], "/", y = s[j]^(1/d[i]))
	}
}




#######################################
# Case: approximate search in 
# multidimensional assignment problem
#######################################


# Dimensionwise heuristic: find best canonical 
# tensor one dataset at a time while keeping the 
# other canonical tensors fixed
 
if (norm == "block" && search == "approximate") {
	## Initialization @@@@ does not do what it's supposed to do. FIX IT
	objective <- 0
	part <- matrix(0, m, m) # partial objective values
	for (i in 1:m)
	for (j in 1:i)
		part[i,j] <- part[j,i] <- 
			w[i,j] * mean(score[,i,j] * score[,j,i])
	diag(part) <- -Inf
	assignment <- integer(m)
	for (i in 1:m) {
		idx <- arrayInd(which.max(part), c(m,m))
		assignment[idx[1]] <- idx[2] 
		part[idx[1],] <- part[,idx[2]] <- -Inf 
	}	
	## Dimensionwise heuristic
	count <- 0
	repeat{
		count <- count + 1
		assignment.old <- assignment
		for (i in 1:m) { ## dimension to update
			part <- matrix(0, m, m)
			for (j in 1:m) { # candidate assignment
				assignment[i] <- j
				for (k in 1:m) # go through all datasets 
					part[j,k] <- w[i,k] * 
					mean(score[,i,j] * score[,k,assignment[k]])	
			}
			assignment[i] <- which.max(colSums(part))
		}
		if (count == 10 || 
			identical(assignment, assignment.old)) break
	}
	v <- v[cbind(1:m, assignment)]
}



######################################
# Case: exhaustive search in
# multidimensional assignment problem
######################################


if (norm == "block" && search == "exhaustive") {
	part <- array(dim = rep(m, 4))
	for (i in 1:m) 
	for (j in 1:i) 
	{			
		part[i, j, , ] <- (w[i,j] / n) * 
			crossprod(score[,i,], score[,j,])
		part[j, i, , ] <- t(part[i, j, , ])
	}
	objective <- array(0, dim = rep(m, m))
	for (i in 1:m)
	for (j in 1:m)
	{
		if (i == j) { 
			# self-contribution c(i,i) sum_t <v_i, x_it >^2
			perm <- 1:m
			perm[c(1,i)] <- c(i,1)
			objective <- aperm(objective, perm)
			dim(objective) <- c(m, m^(m-1))
			objective <- objective + diag(part[i,i,,])
			dim(objective) <- rep(m, m)
			objective <- aperm(objective, perm)
		} else { 
			# c(i,j) sum_t <v_i, x_it> <v_j, x_jt>  for i != j 
			perm <- 1:m
			perm[1:2] <- c(i,j)
			perm[-(1:2)] <- setdiff(1:m, c(i,j))
			objective <- aperm(objective, perm)
			dim(objective) <- c(m^2, m^(m-2))
			objective <- objective + as.vector(part[i,j,,])
			dim(objective) <- rep(m, m)
			iperm <- integer(m)
			iperm[perm] <- 1:m
			objective <- aperm(objective, iperm)
		}
	}
	assignment <- arrayInd(which.max(objective), rep(m,m))
	dim(assignment) <- NULL
	v <- v[cbind(1:m, assignment)]
}



#################
# Postprocessing
#################


## Drop singleton dimensions
for (i in 1:m)
	v[[i]] <- lapply(v[[i]], drop)

## Put back canonical tensors for constant datasets
if (any(cnst.set)) {
	vv <- vector("list", m0)
	vv[!cnst.set] <- v
	for (i in which(cnst.set))
		vv[[i]] <- lapply(p0[[i]], numeric)
	v <- vv
}

v

}





