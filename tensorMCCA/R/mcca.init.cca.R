mcca.init.cca <- function(x, k = NULL, w = 1, 
	objective = c("cov", "cor"), scale = c("block", "global"), 
	center = TRUE, search = c("exhaustive", "approximate"),
	maxit = 100, tol = .0001)
{
	
################
# Preprocessing
################

## Check argument x
test <- check.arguments(x, w = w)
eps <- 1e-14

## Data dimensions
m <- length(x) 
dimx <- lapply(x, dim) 
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
pp <- sapply(p, prod)
n <- tail(dimx[[1]], 1)

## Check for constant datasets
cnst.set <- sapply(x, function(xx) all(abs(xx - xx[1]) < eps))
if (sum(cnst.set | pp == 1) >= (m-1)) {
	v <- vector("list", m)
	for (i in 1:m) 
		v[[i]] <- lapply(p[[i]], function(len) rep(1/sqrt(len), len))
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
k <- if (is.null(k)) pmin(pp, n) else rep_len(k, m)

## Scaling constraints
objective.type <- match.arg(objective)   
scale <- match.arg(scale) # block or global constraints
if (objective.type == "cor" && scale == "global")
	stop(paste("Argument values 'objective = cor'", 
		"and 'scale = global' are incompatible."))

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
rankx <- integer(m)
for (i in 1:m) {
	rankx[i] <- if (center) {
		qr(matrix(x[[i]], pp[i], n) - xbar[[i]])$rank			
	} else { qr(matrix(x[[i]], pp[i], n))$rank }
}

## Effective rank of SVD 	
k <- pmin(k, rankx)

## Calculate compact/truncated SVD of each unfolded dataset
u <- v <- vector("list", m)
for (i in 1:m) {
	xmat <- if (center) { matrix(x[[i]] - xbar[[i]], pp[i], n)  
		} else { matrix(x[[i]], pp[i], n) } 
	test <- (max(3, 2 * k[i]) <= min(pp[i], n))
	if (test) 
		svdx <- tryCatch(svds(xmat, k[i]), error = function(e) NULL)
	if (!test || is.null(svdx))
		svdx <- svd(xmat, nu = k[i], nv = k[i])
	if (objective.type == "cov") {
		u[[i]] <- svdx$u
		# v[[i]] <- sweep(svdx$v, 2, svdx$d[1:rankx[i]], "/")
		v[[i]] <- sweep(svdx$v, 2, svdx$d[1:rankx[i]], "*")
	} else if (objective.type == "cor") {
		# u[[i]] <- sweep(svdx$u, 2, svdx$d[1:rankx[i]], "*")
		u[[i]] <- sweep(svdx$u, 2, svdx$d[1:rankx[i]], "/")
		v[[i]] <- svdx$v		
	}
}
rm(xmat, svdx)



##########################################
# Calculate CCA for each pair of datasets 	
##########################################

## Canonical vectors
a <- lapply(pp, function(nr) matrix(0, nr, m)) 
# a[[i]][,j] is the first left canonical vector associated to Xi Xj'
# For j > i, a[[j]][,i] is calculated as the first right canonical vector 
# associated to Xi Xj' 
 
for (i in 1:m) {	
	for (j in 1:i) {
		if (j == i) {
			a[[i]][,i] <- u[[i]][,1] 
			next
		}
		test <- all(k[c(i,j)] > 2)
		if (test) 
			svdij <- tryCatch(svds(crossprod(v[[i]], v[[j]]), k = 1), 
				error = function(e) NULL)
		if (!test || is.null(svdij)|| is.null(svdij$u))
			svdij <- svd(crossprod(v[[i]], v[[j]]), nu = 1, nv = 1)		
		a[[i]][,j] <- u[[i]] %*% svdij$u
		a[[j]][,i] <- u[[j]] %*% svdij$v
	}
}



#####################################
# Approximate long canonical vectors 
# by rank-1 tensors 
#####################################



v <- vector("list", m^2)
dim(v) <- c(m, m)
for (i in 1:m) {
	if (objective.type == "cor")
		xbar <- as.vector(rowMeans(x[[i]], dims = d[i]))
	for (j in 1:m) {
		aij <- array(a[[i]][,j], dim = p[[i]])
		v[[i,j]] <- if (objective.type == "cov") {
			tnsr.rk1(aij, scale = TRUE, 	maxit = maxit, tol = tol)
		} else {
			tnsr.rk1.score(aij, cnstr = x[[i]] - xbar, maxit = maxit, 
				tol = tol)
		}
	}
}


# After this stage, all canonical tensor weights 
# are scaled to have unit norm or unit variance 
# as required



#############################
# Calculate canonical scores 
#############################


score <- canon.scores(x, v)
score <- score - rep(colMeans(score), each = n)




#######################################
# Case: approximate search in 
# multidimensional assignment problem
#######################################


# Dimensionwise heuristic: find best canonical 
# tensor one dataset at a time while keeping the 
# other canonical tensors fixed
 
if (search == "approximate") {
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


if (search == "exhaustive") {
	part <- array(dim = rep(m, 4))
	# Indices: 1 = dataset 1, 2 = dataset 2, 
	# 3 = candidate for dataset 1, 4 = candidate for dataset 2
	for (i in 1:m) 
	for (j in 1:i) 
	{			
		part[i, j, , ] <- (w[i,j] / n) * 
			crossprod(score[,i,], score[,j,])
		part[j, i, , ] <- t(part[i, j, , ])
	}
	
	objective <- array(0, dim = rep(m, m + 2))	
	for (i in 1:m)
	for (j in 1:m)
	{
		if (i == j) { 
			# self-contribution c(i,i) sum_t <v_i, x_it >^2
			perm <- 1:(m+2)
			perm[c(3,2+i)] <- c(2+i,3)
			objective <- aperm(objective, perm)
			dim(objective) <- c(m, m, m^m)
			objective[i,i,] <- diag(part[i,i,,])
			dim(objective) <- rep(m, m + 2)
			objective <- aperm(objective, perm)
		} else { 
			# c(i,j) sum_t <v_i, x_it> <v_j, x_jt>  for i != j 
			perm <- 1:(m+2)
			perm[3:4] <- c(2+i,2+j)
			perm[-(1:4)] <- setdiff(3:(m+2), c(2+i,2+j))
			objective <- aperm(objective, perm)
			dim(objective) <- c(m, m, m^m)
			objective[i,j,] <- as.vector(part[i,j,,])
			dim(objective) <- rep(m, m + 2)
			iperm <- integer(m)
			iperm[perm] <- 1:(m+2)
			objective <- aperm(objective, iperm)
		}
	}
	
	flip <- lapply(0:floor(m/2), combn, x = m, simplify = FALSE)
	flip <- unlist(flip, FALSE)
	nflip <- length(flip)
	dim(objective) <- c(m^2, m^m)
	objective.best <- numeric(nflip)
	idx.best <- integer(nflip)	
	for (k in 1:nflip) {
		vec <- rep(1, m)
		vec[flip[[k]]] <- -1
		mask <- as.numeric(tcrossprod(vec))
		vals <- crossprod(mask, objective)
		idx.best[k] <- which.max(vals)
		objective.best[k] <- vals[idx.best[k]]
	}
	k <- which.max(objective.best)
	assignment <- arrayInd(idx.best[k], rep(m,m))
	dim(assignment) <- NULL
	v <- v[cbind(1:m, assignment)]
	for (i in flip[[k]]) 
		v[[i]][[1]] <- - v[[i]][[1]]
}



####################################
# Case: maximize sum of covariances
# with global norm constraint
####################################


# >>> TO BE COMPLETED (OR NOT) <<<



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
		vv[[i]] <- lapply(p0[[i]], function(len) rep(1/sqrt(len), len))
	v <- vv
}

v

}



