mcca.cov <- function(x, r = 1L, w = NULL, scope = c("block", "global"), 
	ortho = c("score", "weight"), optim = c("bca", "grad.scale", 
	"grad.rotate"), init = c("cca", "svd", "random"), maxit = 1000, 
	tol = 1e-6, sweep = c("cyclical", "random"), control = list(), 
	verbose = FALSE)
{

## Check arguments
test <- check.arguments(x = x, w = w, 
	v = if (is.list(init)) init else NULL)
r <- as.integer(r)
stopifnot(r > 0)

## Data dimensions
m <- length(x)
dimx <- lapply(x, dimfun)
d <- sapply(dimx, length) - 1L 
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
p[d == 0] <- 1
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
scope <- match.arg(scope)
ortho <- match.arg(ortho)
optim <- match.arg(optim)
if (is.character(init)) 
	init <- match.arg(init)
sweep <- match.arg(sweep)
if (optim %in% c("grad.scale", "grad.rotate") && r > 1 &&
	(ortho == "score" || scope == "global"))
	stop("The gradient-based methods do not support orthogonality ",
	"constraints on canonical scores nor do they support global ",
	"orthogonality constraints. If 'optim' is set to 'grad.scale' ",
	"or 'grad.rotate', please either set 'r = 1' or set ",
	"'ortho = weight' and 'scope = block'.")

## Set initialization method
if (is.character(init)) 
	mcca.init <- switch(init, cca = mcca.init.cca, 
		svd = mcca.init.svd, random = mcca.init.random)


## Calculate tensor means across last dimension
xbar <- vector("list", m)
uncentered <- logical(m)
for(i in 1:m) {
    xbar[[i]] <- if (d[i] == 0) {
    	mean(x[[i]])
    } else {
	    as.vector(rowMeans(x[[i]], dims = d[i]))
	}
	uncentered[i] <- any(abs(xbar[[i]]) > 1e-16)
}
uncentered <- which(uncentered)

## Adjust number of canonical components 
pp <- sapply(p, prod)
r <- if (ortho == "score" && scope == "block") { 
	min(n - 1, max(pp), r)
} else if (ortho == "score" && scope == "global") {
	min(n - 1, sum(pp), r)
} else if (ortho == "weight" && scope == "block") {
	min(max(pp), r)
} else {
	min(sum(pp), r)
}

## Set orthogonality constraints
if (r > 1) {
	if (ortho == "weight" && scope == "block") {
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
	init.args <- list(k = NULL, w = w, objective = "cov", 
		scope = scope, optim = NULL)
	if (test) {
		names. <- intersect(names(control$init), 
			c("k", "optim", "maxit", "tol"))
		init.args[names.] <- control$init[names.]
	}
} else if (identical(init, "svd")) {
	init.args <- list(objective = "cov", scope = scope)
	if (test) {
		names. <- intersect(names(control$init), 
			c("k", "optim", "maxit", "tol"))
		init.args[names.] <- control$init[names.]
	}
} else if (identical(init, "random")) {
	init.args <- list(objective = "cov", scope = scope)
}

## Create output objects 
v <- vector("list", m * r) # canonical weights
dim(v) <- c(m, r)
block.score <- array(0, c(n, m, r)) # canonical scores
global.score <- matrix(0, n, r) 
objective <- rep(NA, r)
iters <- numeric(r)
call.args <- list(objective.type = "cov", r = NULL, w = w, 
	scope = scope, ortho.type = ortho, ortho.cnstr = NULL, 
	optim = optim, init.method = NULL, init.args = init.args, 
	init.val = NULL, maxit = maxit, tol = tol, sweep = sweep, 
	control = control) 
if (ortho == "weight" && scope == "block" && r > 1)
	call.args$ortho.cnstr <- ortho.mode
if (is.character(init)) { 
	call.args$init.method <- init
} else {
	call.args$init.val <- as.matrix(init)
}

## MAIN LOOP
for (l in 1:r) {	
	
	if (verbose) cat("\n\nMCCA: Component", l, "\n")
	
	## Create copy of original data for centering and deflating 
	xl <- x
	for (i in uncentered) 
		xl[[i]] <- xl[[i]] - xbar[[i]]
	
	## Prepare orthogonality constraints and deflate data
	if (l > 1) {
		# browser()
		if (ortho == "score" && scope == "block") {
			for (i in 1:m)
				ortho.cnstr[[i,l-1]] <- tnsr.vec.prod(xl[[i]], 
					block.score[,i,l-1], d[i]+1)
			xl <- deflate.x(xl, ortho.cnstr[,1:(l-1)])
		} else if (ortho == "score" && scope == "global") { 
			for (i in 1:m)
				ortho.cnstr[[i,l-1]] <- tnsr.vec.prod(xl[[i]], 
					global.score[,l-1], d[i]+1)
			xl <- deflate.x(xl, ortho.cnstr[,1:(l-1)], "global")
		} else if (ortho == "weight" && scope == "block") {
			ortho.cnstr <- set.ortho.mat(v = v[,1:(l-1), drop=FALSE], 
				modes = matrix(ortho.mode[, 1:(l-1), l], m, l-1))
			for (i in 1:m)
				xl[[i]] <- tnsr.mat.prod(xl[[i]], ortho.cnstr$mat[[i]],
					ortho.cnstr$modes[[i]])	
		} 
	}
	
	## Remove datasets that are zero 
	xzero <- sapply(xl, function(a) all(a == 0))
	xnzero <- !xzero
	if (any(xzero)) xl <- xl[xnzero]
	wl <- w[xnzero, xnzero]
	
	## Trivial case: all datasets equal to zero
	if (all(xzero)) {
		v[,l] <- relist(lapply(unlist(p), numeric), p)
		objective[l] <- 0
		if (l > 1 && ortho == "weight" && scope == "block") 
			next else break
	}

	# ## Drop singleton dimensions
	# dimxl <- vector("list", sum(xnzero))
	# for (i in seq_along(xl)) {
		# dimxl[[i]] <- dim(xl[[i]])
		# if (is.null(dimxl[[i]]) || all(dimxl[[i]] > 1)) next
		# xl[[i]] <- drop(xl[[i]])
	# }
	
	## Set canonical weights for datasets equal to zero
	for (i in which(xzero)) 
		v[[i,l]] <- lapply(p[[i]], numeric)

	## Specify orthogonality constraints in optimization
	cnstr <- if (l == 1) {
		NULL
	} else if (scope == "block" && ortho == "weight") { 
		NULL 
	} else if (scope == "global" && ortho == "weight") {
		v[xnzero,1:(l-1)]
	} else if (ortho == "score") {
		ortho.cnstr[xnzero,1:(l-1)] 
  }
	
	## Initialize canonical weights
	if (is.character(init)) {
		init.args$x <- xl
		if (init %in% c("cca", "svd")) 
			init.args$w <- wl	
		# if (l == r) debug(mcca.init)
		v0 <- do.call(mcca.init, init.args)
	} else if (is.list(init)) {
		v0 <- if (is.vector(init)) {
			init[xnzero] 
		} else { 
			init[xnzero, min(l, ncol(init))] 
		}
		if (l > 1 && ortho == "weight" && scope == "block") 
			v0 <- tnsr.rk1.mat.prod(v0, mat = ortho.cnstr$mat[xnzero], 
				modes = ortho.cnstr$modes[xnzero], transpose.mat = FALSE)
	}

	## Enforce scaling and orthogonality constraints
	## @@@@ FINISH FUNCTION 'make.feasible' so that it handle 
	## ALL combination of 'ortho' and 'scope' (= 'scale')
	if (l > 1) {
	  v0 <- tnsr.rk1.ortho(v0, cnstr, maxit = 100L, tol = 1e-6)
	}
	v0 <- scale.v(v0, type = "norm", scope = scope, check.args = FALSE)
	
	

	## Run TMCCA
	out <- if (optim == "bca" && scope == "block") {
		mcca.cov.bca.block(x = xl, v = v0, w = wl, 
			ortho = cnstr, sweep = sweep, maxit = maxit, 
			tol = tol, verbose = verbose)
	} else if (optim == "bca" && scope == "global") { 
		mcca.cov.bca.global(x = xl, v = v0, w = wl, 
			ortho = cnstr, sweep = sweep, maxit = maxit, 
			tol = tol, verbose = verbose)
	} else if (optim == "grad.scale") {
		mcca.gradient.scale(x = xl, v = v0, w = wl, 
			scope = scope, type = "norm", maxit = maxit, 
			tol = tol, verbose = verbose)
	} else if (optim == "grad.rotate") {
		mcca.gradient.rotate(x = xl, v = v0, w = wl, 
			maxit = maxit, tol = tol, verbose = verbose)
	}	
	objective[l] <- out$objective
	block.score[,,l] <- out$score
	global.score[,l] <- rowMeans(block.score[,,l]) 
	iters[l] <- out$iters

	# ## Post-process canonical weights
	# for (i in seq_along(xl)) {
		# dimxli <- dimxl[[i]]
		# if (length(dimxli) <= 2 || all(dimxli > 1)) next
		# vil <- vector("list", length(dimxli))
		# vil[dimxli == 1] <- list(1)
		# vil[dimxli > 1] <- out$v[[i]]		
		# out$v[[i]] <- vil		
	# }	
	
	## Transform back canonical weights to original search space 
	if (l > 1 && ortho == "weight" && scope == "block") {
		out$v <- tnsr.rk1.mat.prod(v = out$v, 
			mat = ortho.cnstr$mat[xnzero], 
			modes = ortho.cnstr$modes[xnzero], 
			transpose.mat = TRUE)
	}	
	
	
	v[xnzero,l] <- out$v
	
	## Check orientation of canonical weights (can the objective be
	## increased by flipping the orientation of some canonical tensors?)
	test <- reorient(out$score, wl)
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


