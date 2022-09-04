mcca.cov <- function(x, r = 1, w = 1, scale = c("block", "global"), 
	ortho = c("score", "weight"), optim = c("bca", "grad.scale", 
	"grad.rotate"), init = c("cca", "svd", "random"), maxit = 1000, 
	tol = 1e-6, sweep = c("cyclical", "random"), control = list(), 
	verbose = FALSE)
{

## Check arguments
test <- if (is.list(init)) { 
	check.arguments(x, init, w)
} else check.arguments(x, NULL, w)
eps <- 1e-14
r0 <- r 
r <- as.integer(r)

## Data dimensions
m <- length(x) 
dimx <- lapply(x, dim) 
d <- sapply(dimx, length) - 1L 
p <- mapply(head, dimx, d, SIMPLIFY = FALSE) 
n <- tail(dimx[[1]], 1) 

## Objective weights
w <- if (length(w) == 1) {
	matrix(1 / (m^2), m, m)
} else {
	(w + t(w)) / (2 * sum(w)) 
}

## Match other arguments
scale <- match.arg(scale)
ortho <- match.arg(ortho)
optim <- match.arg(optim)
sweep <- match.arg(sweep)
if (is.character(init)) init <- match.arg(init)
if (optim == "grad.rotate" && scale == "global")
	stop(paste("If argument 'optim' is equal to 'grad.rotate',",
	"'scale' must be equal to 'block'."))

## Set initialization method
if (is.character(init)) 
	mcca.init <- switch(init, cca = mcca.init.cca, 
		svd = mcca.init.svd, random = mcca.init.random)

## Center data
for(i in 1:m) {
    xbar <- rowMeans(x[[i]], dims = d[i])
    x[[i]] <- x[[i]] - as.vector(xbar)
}
if (ortho == "weight") x0 <- x

## Adjust number of canonical components 
r0 <- r
pp <- sapply(p, prod)
r <- if (ortho == "score" && scale == "block") { 
	min(n - 1, max(pp), r0)
} else if (ortho == "score" && scale == "global") {
	min(n - 1, sum(pp), r0)
} else if (ortho == "weight" && scale == "block") {
	min(max(pp), r0)
} else {
	min(sum(pp), r0)
}

## Set orthogonality constraints
ortho.cnstr <- NULL
if (r > 1) {
	if (ortho == "weight" && scale == "block") {
		ortho.mode <- set.ortho.mode(x, r, 
			cnstr = control$ortho$cnstr, 
			method = control$ortho$method)
	} else if (ortho == "score" && scale == "block") {
		ortho.cnstr <- vector("list", m * (r-1))
		dim(ortho.cnstr) <- c(m, r-1)
	} else if (ortho == "score" && scale == "global") {
		ortho.cnstr <- vector("list", r-1) 
	}
}

## Set up initialization parameters
test <- (is.list(control) && !is.null(control$init))
if (identical(init, "cca")) {
	init.args <- list(k = NULL, w = w, 
		objective = "cov", scale = "block", center = FALSE, 
		search = ifelse(m <= 5, "exhaustive", "approximate"))
	if (test) {
		names. <- intersect(names(control$init), 
			c("k", "scale", "search"))
		init.args[names.] <- control$init[names.]
	}
} else if (identical(init, "svd")) {
	init.args <- list(objective = "cov", center = FALSE)
} else if (identical(init, "random")) {
	init.args <- list(objective = "cov")
}

## Create output objects 
v <- vector("list", m * r) # canonical vectors
dim(v) <- c(m, r)
block.score <- array(0, c(n, m, r)) # canonical scores
global.score <- matrix(0, n, r) 
objective <- iters <- numeric(r)
input <- list(objective = "covariance", r = r0, w = w, 
	scale = scale, ortho = ortho, init = init, maxit = maxit, 
	tol = tol, sweep = sweep, control = control) 



## MAIN LOOP
for (l in 1:r) {	
	
	## Deflate data as required
	if (l > 1) { 	
		if (ortho == "score" && scale == "block") { 
			score <- block.score[,, 1:(l-1), drop = FALSE]
			## Calculate orthogonality constraints explicitly	
			for (i in 1:m)
				ortho.cnstr[[i,l-1]] <- tnsr.vec.prod(x[[i]], 
					score[,i,l-1], d[i]+1)
		} else if (ortho == "score" && scale == "global") { 
			score <- global.score[,1:(l-1)]
			for (i in 1:m)
				ortho.cnstr[[i,l-1]] <- tnsr.vec.prod(x[[i]], 
					score[,l-1], d[i]+1)
		} else if (ortho == "weight" && scale == "block") {
			ortho.cnstr <- set.ortho.mat(v = v[,1:(l-1)], 
				modes = ortho.mode[, 1:(l-1), l])
			x <- deflate.x(x0, v = v[, 1:(l-1)], 
				ortho.mode = ortho.mode[, 1:(l-1), l], 
				check.args = FALSE)	
		} 
	}

	## Initialize canonical weights
	if (is.character(init)) {
		init.args$x <- if (l == 1 || ortho == "weight") {
			x } else { 
			deflate.x(x, score = score, scope = scale, 
				check.args = FALSE)
		} 
		v0 <- do.call(mcca.init, init.args)
	} else if (is.list(init)) {
		v0 <- if (is.vector(init)) {
			init } else { init[, min(l, ncol(init))] }
		if (l > 1 && ortho == "weight" && scale == "block") {
			v0 <- tnsr.rk1.mat.prod(v0, mat = ortho.cnstr$mat, 
				modes = ortho.cnstr$modes, transpose.mat = FALSE)
		}
		v0 <- scale.v(v0, type = "norm", scale = scale, 
			check.args = FALSE)
	}
	
	## Run MCCA and store results
	if (verbose) cat("\n\nMCCA: Component", l, "\n")
	out <- if (optim == "bca" && scale == "block") {
		mcca.cov.bca.block(x = x, v = v0, w = w, 
			ortho = if (ortho == "score" && l > 1) {
			ortho.cnstr[,1:(l-1)] } else NULL, 
			sweep = sweep, maxit = maxit, tol = tol, 
			verbose = verbose)
	} else if (optim == "bca" && scale == "global") { 
		mcca.cov.bca.global(x = x, v = v0, w = w, 
			ortho = if (l == 1) { NULL 
			} else if (ortho == "weight") { v[, 1:(l-1)] 
			} else ortho.cnstr[,1:(l-1)], 
			sweep = sweep, maxit = maxit, tol = tol, 
			verbose = verbose)
	} else if (optim == "grad.scale") {
		mcca.gradient.scale(x = x, v = v0, w = w, 
			scale = scale, type = "norm", maxit = maxit, 
			tol = tol, verbose = verbose)
	} else { # optim == "grad.rotate"
		mcca.gradient.rotate(x = x, v = v0, w = w, 
			maxit = maxit, tol = tol, verbose = verbose)
	}	
	v[,l] <- out$v 
	objective[l] <- out$objective
	block.score[,,l] <- out$y 
	global.score[,l] <- rowMeans(block.score[,,l]) 
	iters[l] <- out$iters

	## Transform back canonical weights to original search space 
	## if needed
	if (l > 1 && ortho == "weight" && scale == "block") {
		v[,l] <- tnsr.rk1.mat.prod(v = v[,l], 
			mat = ortho.cnstr$mat, modes = ortho.cnstr$modes, 
			transpose.mat = TRUE)
	}

	## Monitor objective value
	if (objective[l] <= eps) break 
		
} 

## Re-order results according to objective values if needed
o <- order(objective, decreasing = TRUE)
o <- o[objective[o] > eps]
if (length(o) == 0) o <- 1
if (!identical(o, 1:r)) {
	v <- v[,o]
	block.score <- block.score[,,o]
	global.score <- global.score[,o]
	objective <- objective[o]
	iters <- iters[o]
}

r <- length(o)
if (r == 1) {
	dim(block.score) <- c(n, m)
	dim(global.score) <- NULL
} 

## Clean up scores
block.score[abs(block.score) < eps] <- 0
global.score[abs(global.score) < eps] <- 0

list(v = v, block.score = block.score, global.score = global.score,
	objective = objective, iters = iters, input = input)
}


