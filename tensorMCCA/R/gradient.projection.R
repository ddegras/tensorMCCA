#######################################
# FUNCTIONS FOR GRADIENT ASCENT METHOD
#######################################


# Functions to calculate the gradient of the objective function
# and to project vector solutions onto constraint spaces


objective.gradient <- function(x, v, c)
{
## Data dimensions
m <- length(x) # number of datasets
dimx <- lapply(x, dim) # tensor dimensions for each dataset
p <- lapply(dimx, function(x) x[-length(x)]) 
# tensor dimensions for each dataset, last mode (instances) omitted
d <- sapply(p, length) 
# numbers of image dimensions for each dataset
n <- tail(dimx[[1]], 1) # numbers of instances per dataset

## Calculate gradient components 
## 1) Vectors g_ikt = v_i1 x ... x v_i(k-1) x v_i(k+1) x ... x v_id(i) x X_it
## for i = 1,...,m (dataset), k = 1, ... , d_i (mode), t = 1, ..., n (replication) 	
## 2) Inner products < v_i, X_it > = <v_ik, g_ikt> for all i, t and for any k 
g <- vector("list",m)
scores <- matrix(,n,m)
for (i in 1:m) {
	pi <- p[[i]]
	if (d[i] == 1) { # 1D case
		g[[i]][[1]] <- x[[i]]
		scores[,i] <- crossprod(x[[i]],v[[i]][[1]])
	} else if (d[i] == 2) { # 2D case
		xmat <- aperm(x[[i]], c(1,3,2))
		dim(xmat) <- c(pi[1] * n, pi[2])
		xv <- xmat %*% v[[i]][[2]]
		dim(xv) <- c(pi[1], n)
		g[[i]][[1]] <- xv
		xmat <- x[[i]]
		dim(xmat) <- c(pi[1], pi[2] * n)
		xv <- crossprod(xmat, v[[i]][[1]])
		dim(xv) <- c(pi[2], n)
		g[[i]][[2]] <- xv
		scores[,i] <- crossprod(xv, v[[i]][[2]])
	} else { # 3D case
		xmat <- x[[i]]
		dim(xmat) <- c(pi[1], prod(pi[2:3], n))
		xv <- crossprod(v[[i]][[1]], xmat) # v_i1 x X_i
		dim(xv) <- c(pi[2], pi[3] * n)
		xvv <- crossprod(v[[i]][[2]], xv) # v_i1 x v_i2 x X_i
		dim(xvv) <- c(pi[3], n) 
		g[[i]][[3]] <- xvv
		dim(xv) <- c(pi[2:3], n)
		xv <- aperm(xv, c(2,1,3))
		dim(xv) <- c(pi[3], pi[2] * n)
		xvv <- crossprod(v[[i]][[3]], xv) # v_i1 x v_i3 x X_i
		dim(xvv) <- c(pi[2], n)
		g[[i]][[2]] <- xvv
		xmat <- aperm(x[[i]], c(2,3,1,4))
		dim(xmat) <- c(pi[2], pi[1] * pi[3] * n)
		xv <- crossprod(v[[i]][[2]], xmat)
		dim(xv) <- c(pi[3], pi[1] * n)
		xvv <- crossprod(v[[i]][[3]], xv)
		dim(xvv) <- c(pi[1], n)
		g[[i]][[1]] <- xvv		
		scores[,i] <- crossprod(xvv, v[[i]][[1]])
	}	
}

## Combine the components to form gradient
grad <- vector("list",m)
w <- tcrossprod(scores, c/n) # w_ti = sum_j c_ij < v_j, X_jt > / n
for (i in 1:m) {
	for (k in 1:d[i]) {
		grad[[i]][[k]] <- g[[i]][[k]] %*% w[,i]
	}
}
list(grad = grad, scores = scores)	
}



##########
# Wrapper
##########


# mCCA.gradient <- function(x, v, c, maxit = 1000, tol = 1e-6, 
	# type = c("norm", "var"), scope = c("block", "global"), 
	# balance = TRUE, verbose = TRUE)
# {
	
# ## Data dimensions
# dimx <- lapply(x, dim)
# ndimx <- sapply(dimx, length)
# d <- ndimx - 1
# m <- length(x)
# n <- dimx[[1]][ndimx[1]]
# type <- match.arg(type)   
# scope <- match.arg(scope) 

# objective <- numeric(maxit+1)
# objective[1] <- objective.internal(x, v, c)
# if (verbose) 
	# cat("\nIteration",0,"Objective",objective[1])

# for (it in 1:maxit) {

	# ## Calculate gradient
	# grad <- objective.gradient(x, v, c)
	# y <- grad$scores # y_it = < v_i, X_it > 
	# grad <- grad$grad

	# ## Calculate inner products <g_i, X_it> 
	# ## with g_i partial gradient for v_i	
	# gx <- matrix(, n, m) 
	# for (i in 1:m) {
		# dims <- dimx[[i]]
		# if (d[i] == 1) {
			# gx[,i] <- crossprod(x[[i]], grad[[i]][[1]])
		# } else if (d[i] == 2) {
			# xmat <- x[[i]]
			# dim(xmat) <- c(dims[1], dims[2] * n)
			# xmat <- crossprod(grad[[i]][[1]], xmat)
			# dim(xmat) <- c(dims[2], n)
			# gx[,i] <- crossprod(xmat, grad[[i]][[2]])
		# } else {
			# xmat <- x[[i]]
			# dim(xmat) <- c(dims[1], dims[2] * dims[3] * n)
			# xmat <- crossprod(grad[[i]][[1]], xmat)
			# dim(xmat) <- c(dims[2], dims[3] * n)
			# xmat <- crossprod(grad[[i]][[2]], xmat)
			# dim(xmat) <- c(dims[3], n)
			# gx[,i] <- crossprod(xmat, grad[[i]][[3]])
		# }
	# }
	
	# ## Determine step size(s)
	# ## (use as many step sizes as there are constraints)
	# ## Calculate step sizes that respect variance constraints	
	# if (type == "var") {
			# alpha <- vector("list",m)
			# for (i in 1:m) {
				# if (d[i] == 1) {
					# alpha[[i]] <- max(0, -2 * sum(y[,i]) / sum(gx[,i]^2))
				# } else if (d[i] == 2) {
				# cubic.coef <- 
			# }
			# vv <- v # candidate new solutions
			# for (i in 1:m)
			# for (k in 1:d[i])
				# vv[[i]][[k]] <- vv[[i]][[k]] + alpha[i] * grad[[i]][[k]]
		# } else if (type == "norm" && scope == "global") {
			# alpha <- numeric(max(d))
			# for (k in seq_along(alpha)) {
				# idx <- which(d >= k)		
				# vk <- unlist(lapply(v[idx], "[[", k))
				# gradk <- unlist(lapply(grad[idx], "[[", k))
				# alpha[k] <- -2 * sum(vk * gradk) / sum(vk^2)
			# }			
		
				
			
		# ## Check that these step sizes increase objective
		# yy <- image.scores(x, vv) # candidate new scores
		# cyty <- (c/n) * crossprod(y)
		# cyytyy <- (c/n) * crossprod(yy)
		# cytyy <- (c/n) * crossprod(y, yy)
		# power.set <- unlist(lapply(1:n, 
			# function(k) data.frame(combn(n,k))), 
			# recursive = FALSE, use.names = FALSE)
		# vals <- sapply(power.set, function(idx) sum(cyytyy[idx,idx], 
			# cyty[-idx,-idx], cytyy[-idx,idx], cytyy[idx,-idx]))
		# if (any(vals > 0)) {
			# idx <- power.set[[which.max(vals)]]
			# v[idx] <- vv[idx]
		# }
		# objective[it+1] <- max(objective[it], vals) 	
	
	# } else {
		
		
	# }
	
	# quad.coef.obj <- c(objective[it], # degree 0
		# 2 * sum(colMeans(gx) * colSums(c), colMeans(scores) * rowSums(c)), # degree 1
		# sum(crossprod(gx) * c / n)) # degree 2
	
	# ## Calculate 
	
	
	# if (verbose) 
		# cat("\nIteration",it,"Objective",objective[it+1])

	# ## Optionally balance the canonical vectors of each dataset
	# ## so they have equal norm
	# if (balance) {
		# for (i in 1:m) {
			# nrm <- sapply(v[[i]], function(x) sqrt(sum(x^2)))
			# if (any(nrm == 0)) next
			# d <- length(nrm)
			# s <- prod(nrm)^(1/d) / nrm
			# for (k in 1:d)
				# v[[i]][[k]] <- s[k] * v[[i]][[k]]			
		# }
	# }
	
	# ## Check convergence 
	# if (it > 1 && abs(objective[it+1]-objective[it]) <= 
	    	# tol * max(1,objective[it])) break
# }
# list(v = v, y = y, objective = objective[it+1], iters = it, 
	# trace = objective[1:(it+1)])
# }
