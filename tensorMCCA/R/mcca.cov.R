mcca.cov <- function(x, r = 1L, w = NULL, scale = c("block", "global"), 
	ortho = c("score", "weight"), optim = c("bca", "grad.scale", 
	"grad.rotate"), init = c("cca", "svd", "random"), maxit = 1000, 
	tol = 1e-6, sweep = c("cyclical", "random"), control = list(), 
	verbose = FALSE)
{

## Check arguments
test <- check.arguments(x = x, w = w, 
	v = if (is.list(init)) init else NULL)
eps <- 1e-14
r <- as.integer(r)

## Data dimensions
m <- length(x) 
dimx <- lapply(x, dim) 
d <- sapply(dimx, length) - 1L 
p <- mapply(head, dimx, d, SIMPLIFY = FALSE) 
n <- tail(dimx[[1]], 1) 

## Objective weights
w <- if (is.null(w)) {
	(1 - diag(m)) / m / (m-1)
} else if (length(w) == 1) {
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
if (optim %in% c("grad.scale", "grad.rotate") && 
	ortho == "score" && r > 1)
	stop("The gradient-based methods do not support orthogonality ",
	"constraints on canonical scores. If 'optim' is set to 'grad.scale' ",
	"or 'grad.rotate', 'ortho' must be set to 'weight' or 'r' to 1.")

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
pp <- sapply(p, prod)
r <- if (ortho == "score" && scale == "block") { 
	min(n - 1, max(pp), r)
} else if (ortho == "score" && scale == "global") {
	min(n - 1, sum(pp), r)
} else if (ortho == "weight" && scale == "block") {
	min(max(pp), r)
} else {
	min(sum(pp), r)
}

## Set orthogonality constraints
if (r > 1) {
	if (ortho == "weight" && scale == "block") {
		ortho.mode <- set.ortho.mode(x, r, 
			cnstr = control$ortho$cnstr, 
			method = control$ortho$method)
	} else if (ortho == "score") {
		ortho.cnstr <- vector("list", m * (r-1))
		dim(ortho.cnstr) <- c(m, r-1)
	} 
}

## Set up initialization parameters
test <- (is.list(control) && !is.null(control$init))
init.args <- NULL
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
	init.args <- list(objective = "cov", scale = scale, 
		center = FALSE)
} else if (identical(init, "random")) {
	init.args <- list(objective = "cov")
}

## Create output objects 
v <- vector("list", m * r) # canonical weights
dim(v) <- c(m, r)
block.score <- array(0, c(n, m, r)) # canonical scores
global.score <- matrix(0, n, r) 
objective <- rep(NA, r)
iters <- numeric(r)
call.args <- list(objective.type = "cov", r = NULL, w = w, 
	scale = scale, ortho.type = ortho, ortho.cnstr = NULL, 
	optim = optim, init.method = NULL, init.args = init.args, 
	init.val = NULL, maxit = maxit, tol = tol, sweep = sweep, 
	control = control) 
if (ortho == "weight" && scale == "block" && r > 1)
	call.args$ortho.cnstr <- ortho.mode
if (is.character(init)) { 
	call.args$init.method <- init
} else {
	call.args$init.val <- as.matrix(init)
}

ones <- function(len) rep(1/sqrt(len), len) 

## MAIN LOOP
for (l in 1:r) {	
	
	if (verbose) cat("\n\nMCCA: Component", l, "\n")

	## Reduce data by orthogonal projections if required
	if (l > 1) { 	
		if (ortho == "score" && scale == "block") { 
			score <- block.score[,, 1:(l-1), drop = FALSE]
			## Calculate orthogonality constraints explicitly	
			for (i in 1:m)
				ortho.cnstr[[i,l-1]] <- tnsr.vec.prod(x[[i]], 
					score[,i,l-1], d[i]+1)
		} else if (ortho == "score" && scale == "global") { 
			score <- global.score[,1:(l-1), drop = FALSE]
			for (i in 1:m)
				ortho.cnstr[[i,l-1]] <- tnsr.vec.prod(x[[i]], 
					score[,l-1], d[i]+1)
		} else if (ortho == "weight" && scale == "block") {
			ortho.cnstr <- set.ortho.mat(v = v[,1:(l-1)], 
				modes = ortho.mode[, 1:(l-1), l])
			if (any(ortho.cnstr$vzero)) {
				v[,l] <- relist(lapply(unlist(p), numeric), p)
				next
			}
			x <- mapply(tnsr.mat.prod, x = x0, 
				mat = ortho.cnstr$mat, modes = ortho.cnstr$modes, 
				SIMPLIFY = FALSE)			
			# x <- deflate.x(x0, v = v[, 1:(l-1)], 
				# ortho.mode = ortho.mode[, 1:(l-1), l], 
				# check.args = FALSE)	
		} 
	}
	
	## More data preprocessing
	if (l == 1 || (ortho == "weight" && scale == "block")) {
		## Remove datasets that are zero 
		xzero <- sapply(x, function(a) all(a == 0))
		xnzero <- !xzero
		if (any(xzero)) x <- x[xnzero]
		## Drop singleton dimensions
		if (any(xnzero)) {	
			dimxx <- lapply(x, 
				function(xx) if (is.array(xx)) dim(xx) else c(1,n))
			dimv <- lapply(dimxx, function(x) x[-length(x)])
			for (i in 1:length(x)) {
				if (all(dimv[[i]] > 1)) next
				x[[i]] <- drop(x[[i]])
				if (!is.array(x[[i]]))
					dim(x[[i]]) <- c(1, n)
			}
		}
	}

	## Trivial case: all datasets equal to zero
	if (all(xzero)) {
		v[,l] <- relist(lapply(unlist(p), ones), p)
		objective[l] <- 0
		block.score[,,l] <- 0
		global.score[,l] <- 0
		if (l > 1 && ortho == "weight" && scale == "block") 
			next else break
	}
	
	## Specify canonical weights for datasets equal to zero
	for (i in which(xzero)) {
		v[[i,l]] <- if (scale == "block") {
			lapply(p[[i]], ones)	
		} else lapply(p[[i]], numeric)
	}

	## Trivial case: all but one dataset are zero
	if (sum(xnzero) == 1) {
		i <- which(xnzero)
		v[[i,l]] <- tnsr.rk1.wrapper(x, dimx[[i]])	
		if (l > 1 && ortho == "weight" && scale == "block") {
		v[[i,l]] <- tnsr.rk1.mat.prod(v = v[[i,l]], 
			mat = ortho.cnstr$mat[[i]], modes = ortho.cnstr$modes[[i]], 
			transpose.mat = TRUE)
		}
		block.score[,,l] <- canon.scores(x0, v[,l], FALSE)
		global.score[,l] <- rowMeans(block.score[,,l])
		objective[l] <- sum(crossprod(block.score[,,l]) * w) / n
		next	
	}

	## Initialize canonical weights
	if (is.character(init)) {
		init.args$x <- if (l == 1) {
			x
		} else if (ortho == "weight" && scale == "block") {
			x 	
		} else if (ortho == "score" && scale == "block") {
			deflate.x(x, score = score[,xnzero,], 
				scope = "block", check.args = FALSE)
		} else if (ortho == "score" && scale == "global") {
			deflate.x(x, score = score, scope = "global", 
				check.args = FALSE)
		} else if (ortho == "weight" && scale == "global") {
			deflate.x(x, v = v[,1:(l-1)], check.args = FALSE)
		}
		if (init == "cca") init.args$w <- w[xnzero, xnzero]	
		v0 <- do.call(mcca.init, init.args)
	} else if (is.list(init)) {
		v0 <- if (is.vector(init)) {
			init[xnzero] } else { init[xnzero, min(l, ncol(init))] }
		if (l > 1 && ortho == "weight" && scale == "block") 
			v0 <- tnsr.rk1.mat.prod(v0, mat = ortho.cnstr$mat[xnzero], 
				modes = ortho.cnstr$modes[xnzero], transpose.mat = FALSE)
		v0 <- scale.v(v0, type = "norm", scale = scale, check.args = FALSE)
	}
	if (l > 1 && ortho == "score" && scale == "block") {
		v0 <- tnsr.rk1.ortho(v0 = v0, maxit = 100L, tol = 1e-6, 
			ortho = ortho.cnstr[xnzero, 1:(l-1), drop = FALSE])
		v0 <- scale.v(v0, check.args = FALSE)
	}	
	if (all(unlist(v0) == 0)) {
		v[,l] <- relist(lapply(unlist(p), numeric), p)
		next
	}

	## Run MCCA
	out <- if (optim == "bca" && scale == "block") {
		mcca.cov.bca.block(x = x, v = v0, w = w[xnzero,xnzero], 
			ortho = if (ortho == "score" && l > 1) {
			ortho.cnstr[xnzero,1:(l-1)] } else NULL, 
			sweep = sweep, maxit = maxit, tol = tol, 
			verbose = verbose)
	} else if (optim == "bca" && scale == "global") { 
		mcca.cov.bca.global(x = x, v = v0, w = w[xnzero,xnzero], 
			ortho = if (l == 1) { NULL 
			} else if (ortho == "weight") { v[xnzero, 1:(l-1)] 
			} else ortho.cnstr[,1:(l-1)], sweep = sweep, 
			maxit = maxit, tol = tol, verbose = verbose)
	} else if (optim == "grad.scale") {
		mcca.gradient.scale(x = x, v = v0, w = w[xnzero,xnzero], 
			scale = scale, type = "norm", maxit = maxit, 
			tol = tol, verbose = verbose)
	} else if (optim == "grad.rotate") {
		mcca.gradient.rotate(x = x, v = v0, w = w[xnzero,xnzero], 
			maxit = maxit, tol = tol, verbose = verbose)
	}	
	objective[l] <- out$objective
	block.score[,,l] <- out$score
	global.score[,l] <- rowMeans(block.score[,,l]) 
	iters[l] <- out$iters
	
	## Post-process canonical weights
	for (i in 1:length(x)) {
		singleton <- (dimv[[i]] == 1)
		if (!any(singleton)) next			
		vil <- replicate(length(dimv[[i]]), 1, FALSE)
		if (all(singleton)) {
			vil[[1]] <- out$v[[i]][[1]]
		} else {
			vil[!singleton] <- out$v[[i]]			
		}
		out$v[[i]] <- vil		
	}
	
	## Transform back canonical weights to original search space 
	if (l > 1 && ortho == "weight" && scale == "block") {
		out$v <- tnsr.rk1.mat.prod(v = out$v, 
			mat = ortho.cnstr$mat[xnzero], modes = ortho.cnstr$modes[xnzero], 
			transpose.mat = TRUE)
	}	
	v[xnzero,l] <- out$v
	
	## Check orientation of canonical weights (can the objective be
	## increased by flipping the orientation of some canonical tensors?)
	test <- reorient(out$score, w[xnzero,xnzero])
	if (any(test$flip)) {
		objective[l] <- test$objective
		for (i in which(test$flip)) {
			v[[i,l]][[1]] <- (-v[[i,l]][[1]])
			block.score[,i,l] <- -block.score[,i,l]
		}
		global.score[,l] <- rowMeans(block.score[,,l])
	}
	
} 

## Retain feasible solutions
feasible <- !is.na(objective)
if (!all(feasible)) {
	objective <- objective[feasible]
	v <- v[, feasible, drop = FALSE]
	block.score <- block.score[,, feasible, drop = FALSE]
	global.score <- global.score[, feasible, drop = FALSE]
	iters <- iters[feasible]
	r <- sum(feasible)
}

## Re-order results according to objective values
o <- order(objective, decreasing = TRUE)
if (!identical(o, 1:r)) {
	v <- v[,o, drop = FALSE]
	block.score <- block.score[,,o, drop = FALSE]
	global.score <- global.score[,o, drop = FALSE]
	objective <- objective[o]
	iters <- iters[o]
}

call.args$r <- r
if (ortho == "score" && r > 1) 
	call.args$ortho.cnstr <- ortho.cnstr

list(v = v, block.score = block.score, global.score = global.score,
	objective = objective, iters = iters, call.args = call.args)
}


