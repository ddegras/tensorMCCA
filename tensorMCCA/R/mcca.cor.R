mcca.cor <- function(x, r = 1L, w = 1, ortho = c("score", "weight"), 
	optim = c("bca", "grad.scale"), init = c("cca", "svd", "random"), 
	maxit = 1000, tol = 1e-6, sweep = c("cyclical", "random"), 
	control = list(), verbose = FALSE)
{

## Check arguments
test <- check.arguments(x = x, w = w, 
	v = if (is.list(init)) init else NULL)
eps <- 1e-14
# r0 <- r 
r <- as.integer(r)
stopifnot(r > 0)

## Data dimensions
m <- length(x) 
dimx <- lapply(x, dim) 
d <- sapply(dimx, length) - 1L 
p <- mapply(head, dimx, d, SIMPLIFY = FALSE) 
pp <- sapply(p, prod)
n <- tail(dimx[[1]], 1) 

## Objective weights
w <- if (length(w) == 1) { 
	matrix(1 / (m^2), m, m) 
} else { 
	(w + t(w)) / (2 * sum(w))
}

## Match other arguments
ortho <- match.arg(ortho)
sweep <- match.arg(sweep)
optim <- match.arg(optim)
if (optim == "gradient.scale" && ortho == "score" && r > 1)
	stop("The gradient-based methods do not support orthogonality ",
	"constraints on canonical scores. If 'optim' is set to 'grad.scale', ",
	"'ortho' must be set to 'weight' or 'r' to 1.")
if (is.character(init)) init <- match.arg(init)

## Set initialization method
if (is.character(init)) 
	mcca.init <- switch(init, cca = mcca.init.cca, 
		svd = mcca.init.svd, random = mcca.init.random)

## Center data
for(i in 1:m) {
    mu <- rowMeans(x[[i]], dims = d[i])
    x[[i]] <- x[[i]] - as.vector(mu)
}
if (ortho == "weight") x0 <- x

## Adjust number of canonical components as needed
r <- switch(ortho, score = min(n - 1, max(pp), r),
	weight = min(max(pp), r)) 

## Set orthogonality constraints
if (r > 1) {
	if (ortho == "weight") {
		## Set tensor modes on which constraints apply
		ortho.mode <- set.ortho.mode(x, r, 
			cnstr = control$ortho$cnstr, 
			method = control$ortho$method)
	} else {
		ortho.cnstr <- vector("list", m * (r-1))
		dim(ortho.cnstr) <- c(m, r-1)
	}
}

## Set up initialization parameters
test <- (is.list(control) && !is.null(control$init))
init.args <- NULL
if (identical(init, "cca")) {
	init.args <- list(k = NULL, w = w, 
		objective = "cor", scope = "block", optim = NULL)
	if (test) {
		names. <- intersect(names(control), 
			c("k", "optim"))
		init.args[names.] <- control[names.]
	}
} else if (identical(init, "svd")) {
	init.args <- list(objective = "cor")
} else if (identical(init, "random")) { 
	init.args <- list(objective = "cor")
}

## Create output objects 
feasible <- logical(r)
v <- vector("list", m * r) # canonical weights
dim(v) <- c(m, r)
block.score <- array(0, c(n, m, r)) 
global.score <- matrix(0, n, r) 
objective <- rep(NA, r)
iters <- numeric(r)
call.args <- list(objective.type = "cor", r = NULL, w = w, 
	scope = "block", ortho.type = ortho, ortho.cnstr = NULL, 
	optim = optim, init.method = NULL, init.args = init.args, 
	init.val = NULL, maxit = maxit, tol = tol, sweep = sweep, 
	control = control) 
if (ortho == "weight" && r > 1)
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
	
	## Prepare orthogonality constraints
	if (l > 1) {
		if (ortho == "weight") {
			## Deflate + reduce data
			ortho.cnstr <- set.ortho.mat(v = v[, 1:(l-1)], 
				modes = ortho.mode[, 1:(l-1), l])
			if (any(ortho.cnstr$vzero)) {
				v[,l] <- relist(lapply(unlist(p), numeric), p)
				next
			}
			x <- mapply(tnsr.mat.prod, x = x0, 
				mat = ortho.cnstr$mat, modes = ortho.cnstr$modes, 
				SIMPLIFY = FALSE)
		} else if (ortho == "score") {
			## Calculate orthogonality constraints explicitly	
			for (i in 1:m)
				ortho.cnstr[[i,l-1]] <- tnsr.vec.prod(x[[i]], 
					block.score[,i,l-1], d[i]+1)
		}
	}

	## More data preprocessing
	if (l == 1 || ortho == "weight") {
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

	## Trivial case: all or all but one datasets equal to zero
	if (sum(xzero) >= (m-1)) {
		v[,l] <- relist(lapply(unlist(p), numeric), p)
		if (l == 1) objective[1] <- 0
		if (l > 1 && ortho == "weight") next else break
	}
	
	## Specify canonical weights for datasets equal to zero
	for (i in which(xzero)) 
		v[[i,l]] <- lapply(p[[i]], numeric)
		
	## Initialize canonical weights
	if (is.character(init)) {
		init.args$x <- if (ortho == "score" && l > 1) {
			deflate.x(x, score = block.score[,,1:(l-1)], 
				scope = "block", check.args = FALSE)
			} else x
		if (init == "cca") init.args$w <- w[xnzero, xnzero]			
		v0 <- do.call(mcca.init, init.args)
	} else { # case: 'init' is a list
		v0 <- if (is.vector(init)) {
			init[xnzero] } else { init[xnzero, min(l, ncol(init))] }
		if (ortho == "weight" && l > 1) 
			v0 <- tnsr.rk1.mat.prod(v0, mat = ortho.cnstr$mat[xnzero], 
				modes = ortho.cnstr$modes[xnzero], transpose.mat = FALSE)
		v0 <- scale.v(v0, type = "var", x = x, check.args = FALSE)
	}
	if (ortho == "score" && l > 1) {
		v0 <- tnsr.rk1.ortho(v0, ortho.cnstr[xnzero, 1:(l-1), drop = FALSE], 
			maxit = 100L, tol = 1e-6)
		v0 <- scale.v(v0, type = "var", x = x, check.args = FALSE)
	}
	if (all(unlist(v0) == 0)) {
		v[,l] <- relist(lapply(unlist(p), numeric), p)
		next
	}
	
	## Run MCCA
	out <- if (optim == "bca") {
		mcca.cor.bca(x = x, v = v0, w = w[xnzero,xnzero], 
		ortho = if (ortho == "score" && l > 1) {
			ortho.cnstr[xnzero, 1:(l-1), drop = FALSE] } else NULL,
		sweep = sweep, maxit = maxit, tol = tol, verbose = verbose)
	} else { 
		mcca.gradient.scale(x = x, v = v0, w = w[xnzero,xnzero], 
			scope = "block", type = "var", maxit = maxit, 
			tol = tol, verbose = verbose)
	}
	objective[l] <- out$objective
	block.score[,xnzero,l] <- out$score
	global.score[,l] <- rowMeans(block.score[,xnzero,l])	
	iters[l] <- out$iters	
	
	## Reinsert singleton dimensions in canonical weights if needed
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
	if (ortho == "weight" && l > 1) {
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
		idx <- which(xnzero)[test$flip]
		for (i in idx) {
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


