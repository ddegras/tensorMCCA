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
dimx <- lapply(x, dimfun) 
d <- sapply(dimx, length) - 1L 
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
p[d == 0] <- 1
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

## Calculate tensor means across last dimension
xbar <- vector("list", m)
uncentered <- logical(m)
for(i in 1:m) {
    xbar[[i]] <- if (d[i] == 1) {
    	mean(x[[i]])
    } else {
	    as.vector(rowMeans(x[[i]], dims = d[i]))
	}
	uncentered[i] <- any(abs(xbar[[i]]) > 1e-16)
}
uncentered <- which(uncentered)

## Adjust number of canonical components
pp <- sapply(p, prod)
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
	init.args <- list(k = NULL, w = w, objective = "cor", 
		scope = "block", optim = NULL)
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

## MAIN LOOP

for (l in 1:r) {	

	if (verbose) cat("\n\nMCCA: Component", l, "\n")

	## Create copy of original data before centering and deflating 
	xl <- x
	for (i in uncentered) 
		xl[[i]] <- xl[[i]] - xbar[[i]]
	
	## Prepare orthogonality constraints and deflate data
	if (l > 1) {
		if (ortho == "score") {
			for (i in 1:m)
				ortho.cnstr[[i,l-1]] <- tnsr.vec.prod(x[[i]], 
					block.score[,i,l-1], d[i]+1)
			xl <- deflate.x(xl, ortho.cnstr[,1:(l-1)])
		} else if (ortho == "weight") {
			ortho.cnstr <- set.ortho.mat(v = v[, 1:(l-1)], 
				modes = ortho.mode[, 1:(l-1), l])
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

	## Drop singleton dimensions
	dimxl <- vector("list", sum(xnzero))
	for (i in seq_along(xl)) {
		dimxl[[i]] <- dim(xl[[i]])
		if (is.null(dimxl[[i]]) || all(dimxl[[i]] > 1)) next
		dim(xl[[i]]) <- drop(dimxl[[i]])
	}

	## Trivial case: all datasets equal to zero
	if (all(xzero)) {
		v[,l] <- relist(lapply(unlist(p), numeric), p)
		objective[l] <- 0
		if (l > 1 && ortho == "weight") next else break
	}
	
	## Set canonical weights for datasets equal to zero
	for (i in which(xzero)) 
		v[[i,l]] <- lapply(p[[i]], numeric)
		
	## Specify orthogonality constraints in optimization
	cnstr <- if (l == 1 || ortho == "weight") {
		NULL
	} else if (ortho == "score") { 
		ortho.cnstr[xnzero, 1:(l-1), drop = FALSE]
	}

	## Initialize canonical weights
	if (is.character(init)) {
		init.args$x <- xl
		if (init %in% c("cca","svd")) 
			init.args$w <- wl			
		v0 <- do.call(mcca.init, init.args)
	} else if (is.list(init)) { 
		v0 <- if (is.vector(init)) {
			init[xnzero] 
		} else { 
			init[xnzero, min(l, ncol(init))]
		}
		if (l > 1 && ortho == "weight") 
			v0 <- tnsr.rk1.mat.prod(v0, mat = ortho.cnstr$mat[xnzero], 
				modes = ortho.cnstr$modes[xnzero], transpose.mat = FALSE)
	}
		
	## Enforce scaling and orthogonality constraints
	if (l > 1) {
		v0 <- tnsr.rk1.ortho(v0, cnstr, maxit = 100L, tol = 1e-6)
	}
	v0 <- scale.v(v0, type = "var", x = x, check.args = FALSE)
	
	## Run TMCCA
	out <- if (optim == "bca") {
		mcca.cor.bca(x = x, v = v0, w = wl, ortho = cnstr,
		sweep = sweep, maxit = maxit, tol = tol, verbose = verbose)
	} else { 
		mcca.gradient.scale(x = x, v = v0, w = wl, 
			scope = "block", type = "var", maxit = maxit, 
			tol = tol, verbose = verbose)
	}
	objective[l] <- out$objective
	block.score[,xnzero,l] <- out$score
	global.score[,l] <- rowMeans(block.score[,xnzero,l])	
	iters[l] <- out$iters	
	
	## Post-process canonical weights
	for (i in seq_along(xl)) {
		dimxli <- dimxl[[i]]
		if (is.null(dimxli) || all(dimxli > 1)) next
		vil <- vector("list", length(dimxli))
		vil[dimxli == 1] <- list(1)
		vil[dimxli > 1] <- out$v[[i]]		
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


