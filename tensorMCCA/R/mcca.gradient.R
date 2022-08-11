

	
mcca.gradient.retraction.scale <- function(x, v, w, 
	scale, cnstr, maxit, tol, verbose)
{
	
## Data dimensions
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
m <- length(x)
n <- tail(dimx[[1]], 1)

## Set up objective values and objective function
objective <- numeric(maxit + 1L)
objective[1] <- objective.internal(x, v, w)
if (verbose) 
	cat("\nIteration", 0, "Objective", objective[1])
objective.fn <- function(alpha, v, grad, scale, cnstr, x, w)
{
	vnew <- relist(unlist(v) + alpha * unlist(grad), v)
	vnew <- scale.v(v = vnew, scale = scale, cnstr = cnstr, 
		x = x, check.args = FALSE)
	objective.internal(x, vnew, w)
}

## Grid of step sizes
alpha.grid <- c(0, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)
nalpha <- length(alpha.grid)

for (it in 1:maxit) {

	## Calculate gradient 
	grad <- objective.gradient(x, v, w)
	
	## Grid search
	vals <- sapply(alpha.grid, objective.fn, v = v, grad = grad, 
		scale = scale, cnstr = cnstr, x = x, w = w)
	
	## Line search 
	idx <- which.max(vals)
	lb <- alpha.grid[max(1L, idx - 1L)]
	ub <- alpha.grid[min(nalpha, idx + 1L)]
	optim.alpha <- optimize(objective.fn, c(lb, ub), v = v, grad = grad, 
		scale = scale, cnstr = cnstr, x = x, maximum = TRUE, tol = 1e-7)

	## Update solution and objective
	if (optim.alpha$objective > vals[idx]) {
		alpha <- optim.alpha$maximum
		objective[it + 1L] <- optim.alpha$objective
	} else {
		alpha <- alpha.grid[idx]
		objective[it + 1L] <- vals[idx]
	}
	v <- relist(unlist(v) + alpha * unlist(grad), v)
	v <- scale.v(v = vnew, scale = scale, cnstr = cnstr, 
		x = x, check.args = FALSE)
	if (verbose) 
		cat("\nIteration", it, "Objective", objective[it + 1L])

	## Check progress
	increase <- objective[it + 1L] - objective[it]
	if (it > 1 && increase <= tol * max(1, abs(objective[it]))) break
	
}								

list(v = v, y = canon.scores(x, v), objective = objective[it + 1L], 
	iters = it, trace = objective[1:(it+1)])
}






mcca.gradient.retraction.rotation <- function(x, v, w, 
	scale, cnstr, maxit, tol, verbose)
{
	
## Data dimensions
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
m <- length(x)
n <- tail(dimx[[1]], 1)

## Set up objective values and objective function
objective <- numeric(maxit + 1L)
objective[1] <- objective.internal(x, v, w)
if (verbose) 
	cat("\nIteration", 0, "Objective", objective[1])
objective.fn <- function(alpha, v, h, scale, cnstr, x, w)
{
	vnew <- relist(unlist(v) + alpha * unlist(grad), v)
	vnew <- scale.v(v = vnew, scale = scale, cnstr = cnstr, 
		x = x, check.args = FALSE)
	objective.internal(x, vnew, w)
}

## Grid of step sizes
alpha.grid <- c(0, 1e-5, 5e-5, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)
nalpha <- length(alpha.grid)

for (it in 1:maxit) {

	## Calculate gradient 
	grad <- objective.gradient(x, v, w)
	
	## Grid search
	vals <- sapply(alpha.grid, objective.fn, v = v, grad = grad, 
		scale = scale, cnstr = cnstr, x = x, w = w)
	
	## Line search 
	idx <- which.max(vals)
	lb <- alpha.grid[max(1L, idx - 1L)]
	ub <- alpha.grid[min(nalpha, idx + 1L)]
	optim.alpha <- optimize(objective.fn, c(lb, ub), v = v, grad = grad, 
		scale = scale, cnstr = cnstr, x = x, maximum = TRUE, tol = 1e-7)

	## Update solution and objective
	if (optim.alpha$objective > vals[idx]) {
		alpha <- optim.alpha$maximum
		objective[it + 1L] <- optim.alpha$objective
	} else {
		alpha <- alpha.grid[idx]
		objective[it + 1L] <- vals[idx]
	}
	v <- relist(unlist(v) + alpha * unlist(grad), v)
	v <- scale.v(v = vnew, scale = scale, cnstr = cnstr, 
		x = x, check.args = FALSE)
	if (verbose) 
		cat("\nIteration", it, "Objective", objective[it + 1L])

	## Check progress
	increase <- objective[it + 1L] - objective[it]
	if (it > 1 && increase <= tol * max(1, abs(objective[it]))) break
	
}								

list(v = v, y = canon.scores(x, v), objective = objective[it + 1L], 
	iters = it, trace = objective[1:(it+1)])
}


