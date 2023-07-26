

	
mcca.gradient.scale <- function(x, v, w, type = c("norm", "var"), 
	scale = c("block", "global"), maxit, tol, verbose)
{

## Match arguments
type <- match.arg(type)
scale <- match.arg(scale)
	
## Data dimensions
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
m <- length(x)
n <- tail(dimx[[1]], 1)

## Set up objective values and objective function
objective <- numeric(maxit + 1L)
v <- scale.v(v = v, type = type, scale = scale,
	x = x, check.args = FALSE)
objective[1] <- objective.internal(x, v, w)
if (verbose) 
	cat("\nIteration", 0, "Objective", objective[1])
objective.fn <- function(alpha, v, grad, type, scale, x, w)
{
	vnew <- relist(unlist(v) + alpha * unlist(grad), v)
	vnew <- scale.v(v = vnew, scale = scale, type = type, 
		x = x, check.args = FALSE)
	objective.internal(x, vnew, w)
}

## Grid of step sizes
alpha.grid <- c(0, 10^seq.int(-5, 3)) 
nalpha <- length(alpha.grid)

for (it in 1:maxit) {

	## Calculate gradient 
	grad <- objective.gradient(x, v, w)
	
	## Grid search
	vals <- sapply(alpha.grid, objective.fn, v = v, grad = grad, 
		scale = scale, type = type, x = x, w = w)
	
	## Line search 
	idx <- which.max(vals)
	lb <- alpha.grid[max(1L, idx - 1L)]
	ub <- alpha.grid[min(nalpha, idx + 1L)]
	optim.alpha <- optimize(objective.fn, c(lb, ub), v = v, 
		grad = grad, scale = scale, type = type, x = x, w = w, 
		maximum = TRUE, tol = 1e-7)

	## Update solution and objective
	if (optim.alpha$objective > vals[idx]) {
		alpha <- optim.alpha$maximum
		objective[it + 1L] <- optim.alpha$objective
	} else {
		alpha <- alpha.grid[idx]
		objective[it + 1L] <- vals[idx]
	}
	v <- relist(unlist(v) + alpha * unlist(grad), v)
	v <- scale.v(v = v, scale = scale, type = type, 
		x = x, check.args = FALSE)
	if (verbose) 
		cat("\nIteration", it, "Objective", objective[it + 1L])

	## Check progress
	increase <- objective[it + 1L] - objective[it]
	if (it > 1 && increase <= tol * max(1, abs(objective[it]))) break	
}								

list(v = v, score = canon.scores(x, v), objective = objective[it + 1L], 
	iters = it, trace = objective[1:(it+1)])
}






mcca.gradient.rotate <- function(x, v, w, maxit, tol, verbose)
{

## Data dimensions
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
m <- length(x)
n <- tail(dimx[[1]], 1)

## Set up objective values and objective function
objective <- numeric(maxit + 1L)
v <- scale.v(v = v, type = "norm", scale = "block", 
	check.args = FALSE)
objective[1] <- objective.internal(x, v, w)
if (verbose) 
	cat("\nIteration", 0, "Objective", objective[1])
objective.fn <- function(theta, v, h, x, w, fixed)
{
	vnew <- relist(cos(theta) * unlist(v) + 
		sin(theta) * unlist(h), v)
	if (length(fixed) > 0) {
		for (idx in 1:nrow(fixed)) {
			i <- fixed[idx,1]
			k <- fixed[idx,2]
			vnew[[i]][[k]] <- v[[i]][[k]]
		}
	}
	objective.internal(x, vnew, w)
}

## Grid of rotation parameters
ntheta <- 10
theta.grid <- seq(-pi, pi, length.out = ntheta + 1L)[-1L]

## Initialize other quantities
eps <- 1e-14 # numeric tolerance for zero
h <- vector("list", m) # projected gradient
hzero <- matrix(FALSE, m, max(d))


for (it in 1:maxit) {

	## Calculate gradient 
	grad <- objective.gradient(x, v, w)
	
	## Project each partial gradient on unit sphere 
	## orthogonally to corresponding component of current solution
	for (i in 1:m) {
		for (k in 1:d[i]) {
			gik <- grad[[i]][[k]]
			vik <- v[[i]][[k]]
			nrmvik2 <- sum(vik^2)
			vikzero <- (nrmvik2 <= (p[[i]][[k]] * eps^2)) 
			hik <- if (vikzero) gik else {
				gik - (sum(gik * vik) / nrmvik2) * vik }
			nrmhik <- sqrt(sum(hik^2))
			hzero[i,k] <- (nrmhik <= (p[[i]][[k]] * eps))
			h[[i]][[k]] <- if (hzero[i,k]) vik else (hik / nrmhik)
		}
	}
	fixed <- NULL
	if (any(hzero)) {
		vold <- v
		fixed <- which(hzero, TRUE)
	}
	
	## Grid search
	vals <- sapply(theta.grid, objective.fn, v = v, h = h, 
		x = x, w = w, fixed = fixed)
	
	## Line search 
	idx <- which.max(vals)
	lb <- theta.grid[max(1L, idx - 1L)]
	ub <- theta.grid[min(ntheta, idx + 1L)]
	optim.theta <- optimize(objective.fn, c(lb, ub), v = v, 
		h = h, x = x, w = w, fixed = fixed, maximum = TRUE, tol = 1e-7)

	## Update solution and objective
	if (optim.theta $objective > vals[idx]) {
		theta <- optim.theta$maximum
		objective[it + 1L] <- optim.theta$objective
	} else {
		theta <- theta.grid[idx]
		objective[it + 1L] <- vals[idx]
	}

	v <- relist(cos(theta) * unlist(v) + sin(theta) * unlist(h), v)
	if (length(fixed) > 0) {
		for (idx in 1:nrow(fixed)) {
			i <- fixed[idx,1]
			k <- fixed[idx,2]
			v[[i]][[k]] <- vold[[i]][[k]]
		}
	}
	if (verbose) 
		cat("\nIteration", it, "Objective", objective[it + 1L])

	## Check progress
	increase <- objective[it + 1L] - objective[it]
	if (it > 1 && increase <= tol * max(1, abs(objective[it]))) break
	
}								

list(v = v, score = canon.scores(x, v), objective = objective[it + 1L], 
	iters = it, trace = objective[1:(it+1)])
}


